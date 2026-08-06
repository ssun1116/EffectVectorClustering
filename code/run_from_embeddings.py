"""Run the complete pipeline from embeddings.h5ad to the final browser."""

import argparse
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--embeddings", type=pathlib.Path, default=pathlib.Path("/data/Input/embeddings.h5ad"))
    parser.add_argument("--work-dir", type=pathlib.Path, default=ROOT / "pipeline_work")
    parser.add_argument("--keep-intermediates", action="store_true")
    args = parser.parse_args()

    subprocess.run(
        [
            sys.executable,
            str(ROOT / "code" / "cluster_perts_latent_vectors.py"),
            "--embeddings", str(args.embeddings),
            "--work-dir", str(args.work_dir),
        ],
        cwd=ROOT,
        check=True,
    )
    command = [
        sys.executable,
        str(ROOT / "code" / "run_low_resolution_pipeline.py"),
        "--source-html",
        str(args.work_dir / "clustered_pert_effect_vectors.html"),
    ]
    if args.keep_intermediates:
        command.append("--keep-intermediates")
    subprocess.run(command, cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
