#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, Iterator, Tuple


def parse_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    """
    Minimal FASTA parser.
    Yields (header_without_>, sequence_with_no_whitespace).
    """
    header: str | None = None
    seq_chunks: list[str] = []

    with path.open("r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks).replace(" ", "").upper()
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)

    if header is not None:
        yield header, "".join(seq_chunks).replace(" ", "").upper()


def wrap_fasta_sequence(seq: str, width: int = 60) -> Iterable[str]:
    for i in range(0, len(seq), width):
        yield seq[i : i + width]


def remove_1based_inclusive_interval(seq: str, start_1: int, end_1: int) -> str:
    """Remove seq[start_1..end_1] (1-based inclusive)."""
    if start_1 < 1 or end_1 < 1 or start_1 > end_1:
        raise ValueError(f"Invalid 1-based inclusive interval: {start_1}-{end_1}")
    if end_1 > len(seq):
        raise ValueError(
            f"Interval {start_1}-{end_1} exceeds sequence length {len(seq)}"
        )
    return seq[: start_1 - 1] + seq[end_1:]


def main() -> None:
    p = argparse.ArgumentParser(
        description="Remove a 1-based inclusive interval from FASTA sequences."
    )
    p.add_argument("-i", "--input", type=Path, required=True, help="Input FASTA file")
    p.add_argument("-o", "--output", type=Path, required=True, help="Output FASTA file")
    p.add_argument(
        "--start",
        type=int,
        default=1586,
        help="Interval start (1-based inclusive). Default: 1586",
    )
    p.add_argument(
        "--end",
        type=int,
        default=2031,
        help="Interval end (1-based inclusive). Default: 2031",
    )
    p.add_argument(
        "--only-first",
        action="store_true",
        help="Process only the first FASTA record.",
    )
    args = p.parse_args()

    records = list(parse_fasta(args.input))
    if not records:
        raise SystemExit(f"No FASTA records found in: {args.input}")

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open("w", encoding="utf-8") as out:
        for idx, (header, seq) in enumerate(records):
            spliced = remove_1based_inclusive_interval(seq, args.start, args.end)
            new_header = f"{header} splice:{args.start}-{args.end}"

            out.write(f">{new_header}\n")
            for line in wrap_fasta_sequence(spliced, width=60):
                out.write(line + "\n")

            print(
                f"{header}: original_len={len(seq)} spliced_len={len(spliced)} "
                f"removed={args.end - args.start + 1}"
            )

            if args.only_first:
                break


if __name__ == "__main__":
    main()

