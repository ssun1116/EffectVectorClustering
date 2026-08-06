import argparse
import pathlib

import pandas as pd

from helpers import clustering_helpers as unbiased
from helpers import module_browser


ROOT = pathlib.Path(__file__).resolve().parents[1]
OUTPUTS = ROOT / "analysis_outputs"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=pathlib.Path, required=True)
    parser.add_argument(
        "--output-tables",
        type=pathlib.Path,
        default=OUTPUTS / "tables" / "low_resolution_modules.xlsx",
    )
    parser.add_argument(
        "--output-html",
        type=pathlib.Path,
        default=OUTPUTS / "html" / "low_resolution_module_browser.html",
    )
    args = parser.parse_args()
    args.output_tables.parent.mkdir(parents=True, exist_ok=True)
    args.output_html.parent.mkdir(parents=True, exist_ok=True)
    labels, matrix = unbiased.load_plotly_heatmap(args.input)
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
    with pd.ExcelWriter(args.output_tables) as writer:
        pd.DataFrame(summaries).to_excel(writer, sheet_name="summary", index=False)
        pd.DataFrame(member_rows).to_excel(writer, sheet_name="members", index=False)
    module_browser.write_browser(
        args.input,
        summaries,
        member_rows,
        args.output_html,
    )

    print(f"Generated {len(summaries)} unsupervised low-resolution modules")
    print(args.output_tables.resolve())
    print(args.output_html.resolve())
    for row in summaries:
        print(
            f"module {row['module_id']:>2} {row['start']:>4}:{row['end']:<4} "
            f"n={row['size']:<3} within={row['within_mean']:.3f} "
            f"theme={row['auto_theme']} genes={row['top_genes']}"
        )


if __name__ == "__main__":
    main()
