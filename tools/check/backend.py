#!/usr/bin/env python3
"""D60: no CUDA type, call, header or source outside the CUDA backend crate.

    python3 -m tools.check.backend

The solver crates use only the vendor-neutral interfaces (the backend and
executor traits of foehncraft-fv, foehncraft-flow and foehncraft-bench); the
NVIDIA code lives in crates/foehncraft-cuda (Rust host code and the CUDA C
reference of G1) and the PTX entry points in gpu/ptx-kernels. A backend for
another vendor goes into a crate of its own and is added to ALLOWED.

Fails (exit code 1) when a file outside the allowed directories
- uses the CUDA bindings (`cudarc`), CUDA driver or runtime types and calls
  (`CUdeviceptr`, `CUstream`, `cuLaunchKernel`, `cudaMalloc`, ...), NVRTC,
  CUDA C headers or kernel qualifiers (`__global__`, `__device__`);
- is a CUDA source (`.cu`, `.cuh`);
- is a Cargo manifest of a crate that depends on `cudarc` (the workspace
  manifest only declares the version).
Documentation that names CUDA in prose is not a finding.
"""
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
ALLOWED = ("crates/foehncraft-cuda/", "gpu/")
SKIP_DIRS = {"target", ".git", "node_modules", "__pycache__"}
CODE = [
    re.compile(r"\bcudarc\b"),
    re.compile(r"\bCU(deviceptr|stream|context|module|function|event|graph|result)\b"),
    re.compile(r"\bcu(LaunchKernel|MemAlloc|MemcpyHtoD|MemcpyDtoH|CtxCreate|StreamCreate)\w*\b"),
    re.compile(r"\bcuda(Malloc|Memcpy|Free|Stream\w*|DeviceSynchronize|LaunchKernel|Graph\w*)\b"),
    re.compile(r"\bnvrtc\w*\b"),
    re.compile(r"#\s*include\s*[<\"](cuda|cuda_runtime|nvrtc|cublas|cusparse)\w*\.h"),
    re.compile(r"\b__(global|device|shared|constant)__\b"),
]


def code_lines(path, text):
    """Lines of code: Rust and Python comments and doc comments are prose."""
    for no, line in enumerate(text.splitlines(), 1):
        s = line.strip()
        if path.suffix == ".rs" and s.startswith("//"):
            continue
        if path.suffix in (".py", ".toml", ".sh") and s.startswith("#"):
            continue
        yield no, line


def main():
    findings = []
    for p in ROOT.rglob("*"):
        rel = p.relative_to(ROOT).as_posix()
        if any(part in SKIP_DIRS for part in p.relative_to(ROOT).parts) or not p.is_file():
            continue
        if rel.startswith(ALLOWED) or rel == "tools/check/backend.py":
            continue
        if p.suffix in (".cu", ".cuh"):
            findings.append(f"{rel}: CUDA source outside the backend crate")
            continue
        if p.name == "Cargo.toml":
            text = p.read_text()
            if rel != "Cargo.toml" and re.search(r"^\s*cudarc\s*=", text, re.M):
                findings.append(f"{rel}: depends on cudarc")
            continue
        if p.suffix not in (".rs", ".py", ".sh", ".c", ".h", ".cpp"):
            continue
        text = p.read_text(errors="replace")
        for no, line in code_lines(p, text):
            for pat in CODE:
                if pat.search(line):
                    findings.append(f"{rel}:{no}: {line.strip()}")
                    break
    for f in findings:
        print(f, file=sys.stderr)
    if findings:
        print(f"D60: {len(findings)} use(s) of CUDA outside {', '.join(ALLOWED)} (above)", file=sys.stderr)
        return 1
    print("D60 backend check passed: CUDA only in " + ", ".join(ALLOWED))
    return 0


if __name__ == "__main__":
    sys.exit(main())
