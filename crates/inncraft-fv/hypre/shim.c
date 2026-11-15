/*
 * Inncraft: a thin C shim over the public API of hypre (D20), compiled only
 * with the cargo feature `hypre`. It solves one symmetric system A x = b,
 * given in CSR form on the host, with hypre's PCG preconditioned by one
 * BoomerAMG V-cycle, and reports iterations, residual and timings. It is the
 * reference against which the in-house multigrid of `inncraft-fv::amg` is
 * checked; it is never needed to build or run Inncraft.
 *
 * Only hypre's documented user interface is called (HYPRE_IJ*, HYPRE_ParCSR*,
 * HYPRE_BoomerAMG*, HYPRE_PCG*); no hypre code is copied. hypre is licensed
 * under MIT or Apache-2.0.
 *
 * hypre must be built without MPI (cmake -DHYPRE_ENABLE_MPI=OFF); with
 * -DHYPRE_ENABLE_CUDA=ON the setup and the solve run on the GPU.
 */
#include <stdlib.h>
#include <time.h>
#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_krylov.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_utilities.h"

#ifndef HYPRE_SEQUENTIAL
#error "Inncraft's hypre reference needs a sequential hypre build: cmake -DHYPRE_ENABLE_MPI=OFF"
#endif

typedef struct {
    double strong_threshold; /* BoomerAMG strength threshold */
    int coarsen_type;        /* e.g. 8 PMIS, 10 HMIS */
    int interp_type;         /* e.g. 6 extended+i, 14 extended */
    int relax_type;          /* e.g. 6 symmetric hybrid Gauss-Seidel, 18 l1-Jacobi */
    int num_sweeps;
    int agg_num_levels;      /* levels of aggressive coarsening */
    int pmax;                /* largest number of interpolation entries per row */
    int max_coarse;          /* largest coarsest level */
    int print_level;
} inncraft_hypre_options;

typedef struct {
    int iterations;
    int status;              /* hypre's error flag after the solve */
    double final_rel_res;    /* ||r|| / ||b|| reported by PCG */
    double setup_seconds;
    double solve_seconds;
} inncraft_hypre_result;

static double now_seconds(void)
{
    struct timespec t;
    timespec_get(&t, TIME_UTC);
    return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}

/* 1 when the hypre library runs on a GPU. */
int inncraft_hypre_uses_gpu(void)
{
#if defined(HYPRE_USING_GPU)
    return 1;
#else
    return 0;
#endif
}

/* The hypre version string. */
const char *inncraft_hypre_version(void)
{
    return HYPRE_RELEASE_VERSION;
}

int inncraft_hypre_solve(int n, const int *row_ptr, const int *col, const double *val,
                         const double *rhs, double *x, double tol, int max_iter,
                         const inncraft_hypre_options *o, inncraft_hypre_result *res)
{
    MPI_Comm comm = 0;
    HYPRE_IJMatrix A;
    HYPRE_IJVector b, u;
    HYPRE_ParCSRMatrix pA;
    HYPRE_ParVector pb, pu;
    HYPRE_Solver pcg, amg;
    HYPRE_Int i, ilower = 0, iupper = n - 1;
    HYPRE_Int *ncols, *rows, *cols;
    double t0;

    if (!HYPRE_Initialized()) {
        HYPRE_Initialize();
    }
#if defined(HYPRE_USING_GPU)
    HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE);
    HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE);
#endif

    ncols = (HYPRE_Int *)malloc(sizeof(HYPRE_Int) * (size_t)n);
    rows = (HYPRE_Int *)malloc(sizeof(HYPRE_Int) * (size_t)n);
    cols = (HYPRE_Int *)malloc(sizeof(HYPRE_Int) * (size_t)row_ptr[n]);
    for (i = 0; i < n; i++) {
        rows[i] = i;
        ncols[i] = row_ptr[i + 1] - row_ptr[i];
    }
    for (i = 0; i < row_ptr[n]; i++) {
        cols[i] = col[i];
    }

    /* matrix and vectors are filled on the host and moved to the device */
    HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &A);
    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize_v2(A, HYPRE_MEMORY_HOST);
    HYPRE_IJMatrixSetValues(A, n, ncols, rows, cols, val);
    HYPRE_IJMatrixAssemble(A);
    HYPRE_IJVectorCreate(comm, ilower, iupper, &b);
    HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize_v2(b, HYPRE_MEMORY_HOST);
    HYPRE_IJVectorSetValues(b, n, rows, rhs);
    HYPRE_IJVectorAssemble(b);
    HYPRE_IJVectorCreate(comm, ilower, iupper, &u);
    HYPRE_IJVectorSetObjectType(u, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize_v2(u, HYPRE_MEMORY_HOST);
    for (i = 0; i < n; i++) {
        x[i] = 0.0;
    }
    HYPRE_IJVectorSetValues(u, n, rows, x);
    HYPRE_IJVectorAssemble(u);
#if defined(HYPRE_USING_GPU)
    HYPRE_IJMatrixMigrate(A, HYPRE_MEMORY_DEVICE);
    HYPRE_IJVectorMigrate(b, HYPRE_MEMORY_DEVICE);
    HYPRE_IJVectorMigrate(u, HYPRE_MEMORY_DEVICE);
#endif
    HYPRE_IJMatrixGetObject(A, (void **)&pA);
    HYPRE_IJVectorGetObject(b, (void **)&pb);
    HYPRE_IJVectorGetObject(u, (void **)&pu);

    HYPRE_ParCSRPCGCreate(comm, &pcg);
    HYPRE_PCGSetMaxIter(pcg, max_iter);
    HYPRE_PCGSetTol(pcg, tol);
    HYPRE_PCGSetTwoNorm(pcg, 1);
    HYPRE_PCGSetPrintLevel(pcg, o->print_level);
    HYPRE_PCGSetLogging(pcg, 1);

    HYPRE_BoomerAMGCreate(&amg);
    HYPRE_BoomerAMGSetStrongThreshold(amg, o->strong_threshold);
    HYPRE_BoomerAMGSetCoarsenType(amg, o->coarsen_type);
    HYPRE_BoomerAMGSetInterpType(amg, o->interp_type);
    HYPRE_BoomerAMGSetRelaxType(amg, o->relax_type);
    HYPRE_BoomerAMGSetNumSweeps(amg, o->num_sweeps);
    HYPRE_BoomerAMGSetAggNumLevels(amg, o->agg_num_levels);
    HYPRE_BoomerAMGSetPMaxElmts(amg, o->pmax);
    HYPRE_BoomerAMGSetMaxCoarseSize(amg, o->max_coarse);
    HYPRE_BoomerAMGSetPrintLevel(amg, o->print_level);
    HYPRE_BoomerAMGSetTol(amg, 0.0);
    HYPRE_BoomerAMGSetMaxIter(amg, 1);
#if defined(HYPRE_USING_GPU)
    HYPRE_BoomerAMGSetKeepTranspose(amg, 1);
#endif
    HYPRE_PCGSetPrecond(pcg, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
                        (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, amg);

    t0 = now_seconds();
    HYPRE_ParCSRPCGSetup(pcg, pA, pb, pu);
    res->setup_seconds = now_seconds() - t0;
    t0 = now_seconds();
    HYPRE_ParCSRPCGSolve(pcg, pA, pb, pu);
    res->solve_seconds = now_seconds() - t0;
    res->status = (int)HYPRE_GetError();
    HYPRE_ClearAllErrors();

    {
        HYPRE_Int its;
        HYPRE_Real rr;
        HYPRE_PCGGetNumIterations(pcg, &its);
        HYPRE_PCGGetFinalRelativeResidualNorm(pcg, &rr);
        res->iterations = (int)its;
        res->final_rel_res = (double)rr;
    }

#if defined(HYPRE_USING_GPU)
    HYPRE_IJVectorMigrate(u, HYPRE_MEMORY_HOST);
#endif
    HYPRE_IJVectorGetValues(u, n, rows, x);

    HYPRE_BoomerAMGDestroy(amg);
    HYPRE_ParCSRPCGDestroy(pcg);
    HYPRE_IJMatrixDestroy(A);
    HYPRE_IJVectorDestroy(b);
    HYPRE_IJVectorDestroy(u);
    free(ncols);
    free(rows);
    free(cols);
    return 0;
}
