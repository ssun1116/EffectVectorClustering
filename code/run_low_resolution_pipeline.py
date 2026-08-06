"""Run the HTML-only final low-resolution module-browser pipeline."""

import argparse
import os
import pathlib
import shutil
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
OUTPUTS = ROOT / "analysis_outputs"
LOCAL_STEPS = [
    CODE / "generate_low_resolution_browser.py",
    CODE / "run_go_enrichment.py",
]
INTERMEDIATES = [
    OUTPUTS / "tables" / "low_resolution_module_summary.tsv",
    OUTPUTS / "tables" / "low_resolution_module_members.tsv",
    OUTPUTS / "html" / "low_resolution_module_browser.html",
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source-html",
        type=pathlib.Path,
        default=ROOT / "pipeline_work" / "clustered_pert_effect_vectors.html",
        help="self-contained clustered Plotly heatmap",
    )
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="keep generated module tables and intermediate browsers",
    )
    args = parser.parse_args()
    child_env = os.environ.copy()
    child_env["PYTHONDONTWRITEBYTECODE"] = "1"
    source = args.source_html.resolve()
    if not source.exists():
        raise FileNotFoundError(source)
    child_env["EFFECT_VECTOR_HTML"] = str(source)
    child_env["PLOTLY_RUNTIME_HTML"] = str(source)
    for script in LOCAL_STEPS:
        print(f"\n==> {script.name}", flush=True)
        subprocess.run(
            [sys.executable, str(script)], cwd=ROOT, check=True, env=child_env
        )
    source_browser = OUTPUTS / "html" / "low_resolution_go_module_browser.html"
    final_browser = OUTPUTS / "html" / "final_low_resolution_module_browser.html"
    shutil.copyfile(source_browser, final_browser)
    shutil.copyfile(
        OUTPUTS / "tables" / "low_resolution_module_go_summary.tsv",
        OUTPUTS / "tables" / "final_low_resolution_module_summary.tsv",
    )
    print("\n==> final low-resolution browser", flush=True)
    if not args.keep_intermediates:
        for path in INTERMEDIATES:
            path.unlink(missing_ok=True)
        print("\nRemoved generated intermediate files", flush=True)
    print(
        "\nOutput:",
        ROOT / "analysis_outputs" / "html" / "final_low_resolution_module_browser.html",
    )


if __name__ == "__main__":
    main()
