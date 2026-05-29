#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Tuple

CODON_TABLE_STD = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def parse_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    """Minimal FASTA parser yielding (header_without_>, sequence_no_whitespace)."""
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


def wrap(seq: str, width: int = 60) -> Iterable[str]:
    for i in range(0, len(seq), width):
        yield seq[i : i + width]


def slice_1based_inclusive(seq: str, start_1: int, end_1: int) -> str:
    if start_1 < 1 or end_1 < 1 or start_1 > end_1:
        raise ValueError(f"Invalid 1-based inclusive range: {start_1}-{end_1}")
    if end_1 > len(seq):
        raise ValueError(f"Range {start_1}-{end_1} exceeds sequence length {len(seq)}")
    return seq[start_1 - 1 : end_1]


def trim_cds_to_multiple_of_3(cds: str) -> Tuple[str, int]:
    r = len(cds) % 3
    if r == 0:
        return cds, 0
    return cds[: len(cds) - r], r


def translate_cds_standard(cds: str) -> str:
    """
    Translate CDS using the standard genetic code.
    - Any codon containing non-ACGT gets 'X'
    - Caller should ensure len(cds) % 3 == 0 (or trim beforehand)
    """
    protein_chars: list[str] = []
    for i in range(0, len(cds) - (len(cds) % 3), 3):
        codon = cds[i : i + 3]
        aa = CODON_TABLE_STD.get(codon, "X")
        protein_chars.append(aa)
    return "".join(protein_chars)


@dataclass(frozen=True)
class QC:
    transcript_len: int
    cds_start_1: int
    cds_end_1: int
    utr5_len: int
    cds_len_raw: int
    cds_len_trimmed: int
    cds_trimmed_bp: int
    utr3_len: int
    cds_mod3_raw: int
    cds_mod3_trimmed: int
    cds_start_codon: str
    cds_stop_codon: str
    protein_len: int
    internal_stop_count: int
    ends_with_stop: bool


def count_internal_stops(protein: str) -> int:
    # Internal stops are "*" excluding the last position.
    if not protein:
        return 0
    return protein[:-1].count("*")


def evaluate_one(
    header: str, transcript: str, cds_start_1: int, cds_end_1: int, trim_cds: bool
) -> Tuple[str, str, str, str, QC]:
    utr5 = transcript[: cds_start_1 - 1]
    cds_raw = slice_1based_inclusive(transcript, cds_start_1, cds_end_1)
    utr3 = transcript[cds_end_1 :]

    cds = cds_raw
    trimmed_bp = 0
    if trim_cds:
        cds, trimmed_bp = trim_cds_to_multiple_of_3(cds_raw)

    protein = translate_cds_standard(cds) if cds else ""
    internal_stops = count_internal_stops(protein)
    ends_with_stop = protein.endswith("*") if protein else False

    cds_start_codon = cds[:3] if len(cds) >= 3 else ""
    cds_stop_codon = cds[-3:] if len(cds) >= 3 else ""

    qc = QC(
        transcript_len=len(transcript),
        cds_start_1=cds_start_1,
        cds_end_1=cds_end_1,
        utr5_len=len(utr5),
        cds_len_raw=len(cds_raw),
        cds_len_trimmed=len(cds),
        cds_trimmed_bp=trimmed_bp,
        utr3_len=len(utr3),
        cds_mod3_raw=len(cds_raw) % 3,
        cds_mod3_trimmed=len(cds) % 3 if cds else 0,
        cds_start_codon=cds_start_codon,
        cds_stop_codon=cds_stop_codon,
        protein_len=len(protein),
        internal_stop_count=internal_stops,
        ends_with_stop=ends_with_stop,
    )

    # Keep headers short but informative for downstream tools.
    base_id = header.split()[0]
    suffix = f"cds:{cds_start_1}-{cds_end_1}"
    utr5_h = f"{base_id} {suffix} feature:five_prime_UTR"
    cds_h = f"{base_id} {suffix} feature:CDS"
    utr3_h = f"{base_id} {suffix} feature:three_prime_UTR"
    prot_h = f"{base_id} {suffix} feature:protein"

    return utr5_h, utr5, cds_h + (f" trimmed:{trimmed_bp}bp" if trimmed_bp else ""), cds, utr3_h, utr3, prot_h, protein, qc


def write_fasta(path: Path, header: str, seq: str, width: int = 60) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write(f">{header}\n")
        for line in wrap(seq, width=width):
            f.write(line + "\n")


def write_qc(path: Path, header: str, qc: QC) -> None:
    # Simple, grep-friendly key=value lines.
    with path.open("w", encoding="utf-8") as f:
        f.write(f"record={header}\n")
        for k, v in qc.__dict__.items():
            f.write(f"{k}={v}\n")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Evaluate TransDecoder CDS coordinates: extract UTRs/CDS, translate CDS, and QC."
    )
    ap.add_argument("-i", "--input", type=Path, required=True, help="Input FASTA (spliced transcript)")
    ap.add_argument("--cds-start", type=int, required=True, help="CDS start (1-based, inclusive)")
    ap.add_argument("--cds-end", type=int, required=True, help="CDS end (1-based, inclusive)")
    ap.add_argument(
        "-o",
        "--outdir",
        type=Path,
        required=True,
        help="Output directory for UTR/CDS/protein/QC files",
    )
    ap.add_argument(
        "--trim-cds",
        action="store_true",
        help="If CDS length is not divisible by 3, trim trailing bases (recommended for translation).",
    )
    ap.add_argument(
        "--only-first",
        action="store_true",
        help="Process only the first FASTA record.",
    )
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    records = list(parse_fasta(args.input))
    if not records:
        raise SystemExit(f"No FASTA records found in: {args.input}")

    for idx, (header, transcript) in enumerate(records):
        (
            utr5_h,
            utr5_seq,
            cds_h,
            cds_seq,
            utr3_h,
            utr3_seq,
            prot_h,
            prot_seq,
            qc,
        ) = evaluate_one(header, transcript, args.cds_start, args.cds_end, args.trim_cds)

        base = header.split()[0]
        safe_base = base.replace("/", "_")
        write_fasta(
            args.outdir / f"{safe_base}.five_prime_UTR.fasta", utr5_h, utr5_seq
        )
        write_fasta(args.outdir / f"{safe_base}.CDS.fasta", cds_h, cds_seq)
        write_fasta(
            args.outdir / f"{safe_base}.three_prime_UTR.fasta", utr3_h, utr3_seq
        )
        write_fasta(args.outdir / f"{safe_base}.protein.fasta", prot_h, prot_seq)
        write_qc(args.outdir / f"{safe_base}.qc.txt", header, qc)

        print(
            f"{base}: transcript_len={qc.transcript_len} "
            f"utr5={qc.utr5_len} cds_raw={qc.cds_len_raw} cds_trimmed={qc.cds_len_trimmed} "
            f"utr3={qc.utr3_len} protein_len={qc.protein_len} internal_stops={qc.internal_stop_count}"
        )

        if args.only_first:
            break


if __name__ == "__main__":
    main()

