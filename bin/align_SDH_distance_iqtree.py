#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Tuple


def parse_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    """Minimal FASTA parser yielding (id, sequence) with whitespace removed."""
    header: str | None = None
    seq_chunks: list[str] = []
    with path.open("r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks).replace(" ", "")
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if header is not None:
        yield header, "".join(seq_chunks).replace(" ", "")


def wrap(seq: str, width: int = 60) -> Iterable[str]:
    for i in range(0, len(seq), width):
        yield seq[i : i + width]


def write_fasta(path: Path, records: List[Tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for line in wrap(seq, 60):
                f.write(line + "\n")


def require_exe(name: str) -> str:
    exe = shutil.which(name)
    if not exe:
        raise SystemExit(
            f"Required executable not found in PATH: {name}\n"
            f"Install it, then re-run."
        )
    return exe


def run(cmd: List[str], cwd: Path | None = None) -> None:
    proc = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        text=True,
        capture_output=True,
    )
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
        run(cmd_for(requested))
    except RuntimeError as e:
        # Common failure when Clustal Omega is built without OpenMP:
        # "FATAL: Cannot change number of threads to X. Clustal Omega was build without OpenMP support."
        msg = str(e)
        if requested != 1 and ("without OpenMP" in msg or "Cannot change number of threads" in msg):
            print(
                f"Clustal Omega has no OpenMP support; retrying with --threads 1 "
                f"(requested {requested})."
            )
            run(cmd_for(1))
        else:
            raise SystemExit(msg)


def run_iqtree(aln_fasta: Path, out_prefix: Path, seq_type: str, threads: int) -> None:
    # Prefer iqtree2 if available, else iqtree.
    exe = shutil.which("iqtree2") or shutil.which("iqtree")
    if not exe:
        raise SystemExit(
            "Required executable not found in PATH: iqtree2 (or iqtree)\n"
            "Install it, then re-run."
        )

    # Keep it simple: model finder + ultrafast bootstrap.
    # can tune flags later without changing the pipeline structure.
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
        run(cmd)
    except RuntimeError as e:
        raise SystemExit(str(e))


def load_alignment(path: Path) -> Tuple[List[str], List[str]]:
    records = list(parse_fasta(path))
    if not records:
        raise SystemExit(f"No records found in alignment: {path}")
    ids = [r[0].split()[0] for r in records]
    seqs = [r[1].upper() for r in records]
    aln_len = {len(s) for s in seqs}
    if len(aln_len) != 1:
        raise SystemExit(f"Alignment sequences have different lengths in: {path}")
    return ids, seqs


@dataclass(frozen=True)
class PairwiseDistance:
    distance: float
    compared_sites: int
    mismatches: int


def p_distance_no_gaps(a: str, b: str) -> PairwiseDistance:
    """
    Gap-ignored p-distance on an alignment:
    - Skip columns where either sequence has a gap '-'
    - Distance = mismatches / compared_sites
    - If compared_sites == 0, distance is NaN (written as empty)
    """
    mismatches = 0
    compared = 0
    for ca, cb in zip(a, b):
        if ca == "-" or cb == "-":
            continue
        compared += 1
        if ca != cb:
            mismatches += 1
    if compared == 0:
        return PairwiseDistance(distance=float("nan"), compared_sites=0, mismatches=0)
    return PairwiseDistance(distance=mismatches / compared, compared_sites=compared, mismatches=mismatches)


def write_distance_matrix_tsv(
    out_path: Path, ids: List[str], seqs: List[str]
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n = len(ids)
    # Precompute for reproducibility and to also enable optional QC outputs later.
    mat: List[List[PairwiseDistance]] = [[PairwiseDistance(0.0, 0, 0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        mat[i][i] = PairwiseDistance(0.0, len(seqs[i].replace("-", "")), 0)
        for j in range(i + 1, n):
            d = p_distance_no_gaps(seqs[i], seqs[j])
            mat[i][j] = d
            mat[j][i] = d

    with out_path.open("w", encoding="utf-8") as f:
        f.write("id\t" + "\t".join(ids) + "\n")
        for i in range(n):
            row = [ids[i]]
            for j in range(n):
                val = mat[i][j].distance
                row.append("" if val != val else f"{val:.6f}")  # NaN check: val != val
            f.write("\t".join(row) + "\n")


def write_distance_qc_tsv(out_path: Path, ids: List[str], seqs: List[str]) -> None:
    """
    Additional QC table with compared sites and mismatches (gap-ignored).
    This helps validate that distances were computed on substantial overlap.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n = len(ids)
    with out_path.open("w", encoding="utf-8") as f:
        f.write("id1\tid2\tcompared_sites\tmismatches\tdistance\n")
        for i in range(n):
            for j in range(i + 1, n):
                d = p_distance_no_gaps(seqs[i], seqs[j])
                dist = "" if d.distance != d.distance else f"{d.distance:.6f}"
                f.write(f"{ids[i]}\t{ids[j]}\t{d.compared_sites}\t{d.mismatches}\t{dist}\n")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Run Clustal Omega alignment, compute gap-ignored distance matrix, "
            "and build an IQ-TREE tree (protein + nucleotide)."
        )
    )
    ap.add_argument(
        "--protein-fasta",
        type=Path,
        default=Path("Files/DQD_SDH/DQD_SDH_possible_pr.fasta"),
        help="Protein FASTA input. Default: Files/DQD_SDH/DQD_SDH_possible_pr.fasta",
    )
    ap.add_argument(
        "--nt-fasta",
        type=Path,
        default=Path("Files/DQD_SDH/DQD_SDH_possible_nt.fasta"),
        help="Nucleotide FASTA input. Default: Files/DQD_SDH/DQD_SDH_possible_nt.fasta",
    )
    ap.add_argument(
        "--outdir",
        type=Path,
        default=Path("Files/DQD_SDH/align_tree"),
        help="Output directory. Default: Files/DQD_SDH/align_tree",
    )
    ap.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Threads for clustalo/iqtree. Default: 4",
    )
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Protein
    pr_aln = args.outdir / "DQD_SDH_possible_pr.clustalo.fasta"
    run_clustalo(args.protein_fasta, pr_aln, args.threads)
    pr_ids, pr_seqs = load_alignment(pr_aln)
    write_distance_matrix_tsv(args.outdir / "DQD_SDH_possible_pr.dist_nogaps.tsv", pr_ids, pr_seqs)
    write_distance_qc_tsv(args.outdir / "DQD_SDH_possible_pr.dist_nogaps.qc.tsv", pr_ids, pr_seqs)
    run_iqtree(pr_aln, args.outdir / "DQD_SDH_possible_pr.iqtree", "protein", args.threads)

    # Nucleotide
    nt_aln = args.outdir / "DQD_SDH_possible_nt.clustalo.fasta"
    run_clustalo(args.nt_fasta, nt_aln, args.threads)
    nt_ids, nt_seqs = load_alignment(nt_aln)
    write_distance_matrix_tsv(args.outdir / "DQD_SDH_possible_nt.dist_nogaps.tsv", nt_ids, nt_seqs)
    write_distance_qc_tsv(args.outdir / "DQD_SDH_possible_nt.dist_nogaps.qc.tsv", nt_ids, nt_seqs)
    run_iqtree(nt_aln, args.outdir / "DQD_SDH_possible_nt.iqtree", "nucleotide", args.threads)

    print(f"Done. Outputs written under: {args.outdir}")


if __name__ == "__main__":
    main()

