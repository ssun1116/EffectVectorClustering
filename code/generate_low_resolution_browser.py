import argparse
import csv
import math
import pathlib

from helpers import clustering_helpers as unbiased
from helpers import module_browser


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
    parser = argparse.ArgumentParser()
    parser.add_argument("--clustered-html", type=pathlib.Path, required=True)
    parser.add_argument("--summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--members-output", type=pathlib.Path, required=True)
    parser.add_argument("--browser-output", type=pathlib.Path, required=True)
    args = parser.parse_args()
    for output in (args.summary_output, args.members_output, args.browser_output):
        output.parent.mkdir(parents=True, exist_ok=True)
    labels, matrix = unbiased.load_plotly_heatmap(args.clustered_html)
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
    write_tsv(args.summary_output, summaries, summary_cols)
    write_tsv(args.members_output, member_rows, member_cols)
    module_browser.write_browser(
        args.clustered_html,
        args.summary_output,
        args.members_output,
        args.browser_output,
    )

    print(f"Generated {len(summaries)} unsupervised low-resolution modules")
    print(args.summary_output.resolve())
    print(args.members_output.resolve())
    print(args.browser_output.resolve())
    for row in summaries:
        print(
            f"module {row['module_id']:>2} {row['start']:>4}:{row['end']:<4} "
            f"n={row['size']:<3} within={row['within_mean']:.3f} "
            f"theme={row['auto_theme']} genes={row['top_genes']}"
        )


if __name__ == "__main__":
    main()
