#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path
from typing import List


def require_exe(name: str) -> str:
    exe = shutil.which(name)
    if not exe:
        raise SystemExit(f"Required executable not found in PATH: {name}")
    return exe


def run_capture(cmd: List[str]) -> None:
    proc = subprocess.run(cmd, text=True, capture_output=True)
    if proc.returncode != 0:
        msg = (proc.stderr or proc.stdout or "").strip()
        raise RuntimeError(
            f"Command failed ({proc.returncode}): {' '.join(cmd)}"
            + (f"\n{msg}" if msg else "")
        )


def run_clustalo(in_fasta: Path, out_fasta: Path, threads: int) -> None:
    clustalo = require_exe("clustalo")
    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    requested = max(1, threads)

    def cmd_for(t: int) -> List[str]:
        return [
            clustalo,
            "-i",
            str(in_fasta),
            "-o",
            str(out_fasta),
            "--outfmt",
            "fasta",
            "--force",
            "--threads",
            str(max(1, t)),
        ]

    try:
        run_capture(cmd_for(requested))
    except RuntimeError as e:
        msg = str(e)
        if requested != 1 and ("without OpenMP" in msg or "Cannot change number of threads" in msg):
            print(
                f"Clustal Omega has no OpenMP support; retrying with --threads 1 "
                f"(requested {requested})."
            )
            run_capture(cmd_for(1))
        else:
            raise SystemExit(msg)


def run_iqtree(aln_fasta: Path, out_prefix: Path, seq_type: str, threads: int) -> None:
    exe = shutil.which("iqtree2") or shutil.which("iqtree")
    if not exe:
        raise SystemExit("Required executable not found in PATH: iqtree2 (or iqtree)")

    cmd = [
        exe,
        "-s",
        str(aln_fasta),
        "-st",
        "AA" if seq_type == "protein" else "DNA",
        "-m",
        "MFP",
        "-bb",
        "1000",
        "-nt",
        str(max(1, threads)),
        "-pre",
        str(out_prefix),
    ]
    try:
        run_capture(cmd)
    except RuntimeError as e:
        raise SystemExit(str(e))


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Align a FASTA with Clustal Omega and build a tree with IQ-TREE."
    )
    ap.add_argument(
        "-i",
        "--input",
        type=Path,
        default=Path("Files/Alignment/For_alignment.fasta"),
        help="Input FASTA. Default: Files/Alignment/For_alignment.fasta",
    )
    ap.add_argument(
        "-o",
        "--output-aln",
        type=Path,
        default=Path("Files/Alignment/ClustalO.aln.fa"),
        help="Output alignment FASTA. Default: Files/Alignment/ClostalO.aln.fa",
    )
    ap.add_argument(
        "--seq-type",
        choices=["protein", "nucleotide"],
        default="protein",
        help="Sequence type for IQ-TREE (-st). Default: protein",
    )
    ap.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Threads. Clustal Omega may auto-fallback to 1 if no OpenMP. Default: 8",
    )
    args = ap.parse_args()

    # IQ-TREE prefix is set to the alignment path so outputs include:
    #   Files/Alignment/ClustalO.aln.fa.contree
    #   Files/Alignment/ClustalO.aln.fa.iqtree
    #   etc.
    run_clustalo(args.input, args.output_aln, args.threads)
    run_iqtree(args.output_aln, args.output_aln, args.seq_type, args.threads)

    print(f"Alignment: {args.output_aln}")
    print(f"Tree (contree): {args.output_aln}.contree")


if __name__ == "__main__":
    main()

