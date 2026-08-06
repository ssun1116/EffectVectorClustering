import csv
import math
import pathlib

import clustering_helpers as unbiased
import module_browser


HERE = pathlib.Path(__file__).resolve().parents[1]
OUT_DIR = HERE / "analysis_outputs"
LOWRES_HTML = OUT_DIR / "html" / "low_resolution_module_browser.html"
LOWRES_SUMMARY = OUT_DIR / "tables" / "low_resolution_module_summary.tsv"
LOWRES_MEMBERS = OUT_DIR / "tables" / "low_resolution_module_members.tsv"


def write_tsv(path, rows, columns):
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            clean = {}
            for key in columns:
                val = row.get(key, "")
                if isinstance(val, float):
                    clean[key] = "" if math.isnan(val) else f"{val:.6g}"
                else:
                    clean[key] = val
            writer.writerow(clean)


def main():
    (OUT_DIR / "tables").mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "html").mkdir(parents=True, exist_ok=True)
    labels, matrix = unbiased.load_plotly_heatmap(unbiased.HTML_PATH)
    segments, prefix = unbiased.recursive_segments(
        matrix,
        min_size=6,
        max_size=180,
        min_gain=0.055,
        min_child_within=0.02,
    )
    summaries, member_rows = unbiased.summarize_segments(labels, matrix, segments, prefix)
    for row in summaries:
        row["auto_theme"] = ""
        row["theme_hits"] = ""
    summary_cols = [
        "module_id",
        "manual_label",
        "auto_theme",
        "theme_hits",
        "start",
        "end",
        "size",
        "within_mean",
        "neighbor_mean",
        "contrast",
        "n_unique_genes",
        "top_gene_fraction",
        "gene_entropy",
        "n_unique_cell_types",
        "top_genes",
        "top_cells",
    ]
    member_cols = ["module_id", "position", "label", "cell_type", "gene"]
    write_tsv(LOWRES_SUMMARY, summaries, summary_cols)
    write_tsv(LOWRES_MEMBERS, member_rows, member_cols)
    module_browser.write_browser()

    print(f"Generated {len(summaries)} unsupervised low-resolution modules")
    print(LOWRES_SUMMARY)
    print(LOWRES_HTML)
    for row in summaries:
        print(
            f"module {row['module_id']:>2} {row['start']:>4}:{row['end']:<4} "
            f"n={row['size']:<3} within={row['within_mean']:.3f} "
            f"theme={row['auto_theme']} genes={row['top_genes']}"
        )


if __name__ == "__main__":
    main()
