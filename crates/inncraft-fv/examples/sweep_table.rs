//! Prints the order-of-accuracy table of docs/verification/stage2.md:
//! `cargo run --release -p inncraft-fv --example sweep_table`.
use inncraft_fv::solver::Settings;
use inncraft_fv::verify::{observed_orders, order_sweep};
use inncraft_mesh::generate::Kind;

fn main() -> anyhow::Result<()> {
    let s = Settings::default();
    println!(
        "| Mesh | Cells (levels) | L2 error (finest) | Observed order | Outer it. (finest) | ns/cell/it. (finest) |"
    );
    println!("|---|---|---:|---|---:|---:|");
    for kind in [Kind::Hex, Kind::Tet, Kind::Poly, Kind::Mixed] {
        for (dim, levels) in [(3u8, vec![1, 2, 3]), (2, vec![2, 3, 4])] {
            let r = order_sweep(
                kind,
                dim,
                4,
                &levels,
                0.15,
                &s,
                "cpu",
                std::path::Path::new(""),
            )?;
            let o = observed_orders(&r, dim);
            let f = r.last().unwrap();
            println!(
                "| {} {dim}D | {} | {:.3e} | {} | {} | {:.0} |",
                kind.name(),
                r.iter()
                    .map(|x| x.cells.to_string())
                    .collect::<Vec<_>>()
                    .join(", "),
                f.l2,
                o.iter()
                    .map(|x| format!("{x:.2}"))
                    .collect::<Vec<_>>()
                    .join(", "),
                f.outer,
                f.ns_per_cell_iteration
            );
        }
    }
    Ok(())
}
