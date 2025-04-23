import argparse
import os
from typing import Optional
from src.aligner.plot import plot_matrix
from src.aligner.core import (
    build_score_matrix,
    traceback as single_traceback,
    trace_all_paths,
)
from src.aligner.io import (
    create_output_dict,
    format_multi_report,
    format_report,
    read_manual,
    read_fasta,
    write_json,
    write_matrix,
    write_report,
)


def parse_args(args=None):
    """
    Parse command-line arguments for the Needleman–Wunsch aligner.
    """
    parser = argparse.ArgumentParser(
        prog="needleman-wunsch",
        description="Global alignment of two sequences using Needleman–Wunsch algorithm",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--input",
        nargs=2,
        metavar=("FASTA1", "FASTA2"),
        help="Paths to two input FASTA files",
    )
    group.add_argument(
        "--manual",
        action="store_true",
        help="Enter sequences manually via prompts",
    )
    parser.add_argument(
        "--match",
        type=int,
        default=1,
        help="Score for a match (default: 1)",
    )
    parser.add_argument(
        "--mismatch",
        type=int,
        default=-1,
        help="Score for a mismatch (default: -1)",
    )
    parser.add_argument(
        "--gap",
        type=int,
        default=-2,
        help="Penalty for a gap (default: -2)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Optional custom filename for text output report",
    )

    parser.add_argument(
        "--all-paths",
        action="store_true",
        help="Enumerate and output all optimal alignments",
    )

    parser.add_argument(
        "--plot",
        type=str,
        default=None,
        help="Optional filename for PNG heatmap output",
    )

    parser.add_argument(
        "--matrix-out",
        type=str,
        default=None,
        help="Filename for raw score matrix output (CSV)",
    )

    parser.add_argument(
        "--json",
        dest="json_out",
        type=str,
        default=None,
        help="Filename for structured JSON output",
    )

    return parser.parse_args(args)


def main():
    args = parse_args()

    # Ensure reports/ directory exists if an output file is specified
    if args.output:
        os.makedirs(os.path.dirname(args.output), exist_ok=True)

    # Ensure plots/ directory exists if a plot file is specified
    if args.plot:
        os.makedirs(os.path.dirname(args.plot), exist_ok=True)

    # 1) Load sequences (manual vs. FASTA)
    if args.manual:
        seq1, seq2 = read_manual()
    else:
        seqs1 = read_fasta(args.input[0])
        seqs2 = read_fasta(args.input[1])
        if len(seqs1) != 1 or len(seqs2) != 1:
            raise ValueError("Each FASTA file must contain exactly one record")
        seq1, seq2 = seqs1[0], seqs2[0]

    # 2) Build the scoring matrix
    matrix = build_score_matrix(
        seq1,
        seq2,
        match=args.match,
        mismatch=args.mismatch,
        gap=args.gap,
    )

    if args.matrix_out:
        os.makedirs(os.path.dirname(args.matrix_out) or ".", exist_ok=True)
        write_matrix(args.matrix_out, matrix)
        print(f"Score matrix CSV written to {args.matrix_out}")

    # 3a) If --all-paths, enumerate and output every optimal alignment
    if args.all_paths:
        all_alignments = trace_all_paths(
            matrix, seq1, seq2, match=args.match, mismatch=args.mismatch, gap=args.gap
        )
        report = format_multi_report(
            seq1, seq2, all_alignments, args.match, args.mismatch, args.gap
        )
        print(report)
        if args.output:
            write_report(args.output, report)
        return

    # 3b) Otherwise, just produce the single optimal path
    aln1, aln2 = single_traceback(
        matrix,
        seq1,
        seq2,
        match=args.match,
        mismatch=args.mismatch,
        gap=args.gap,
    )

    # 4) Format, print, and save the text report
    report = format_report(
        seq1,
        seq2,
        aln1,
        aln2,
        args.match,
        args.mismatch,
        args.gap,
    )

    print(report)
    if args.output:
        write_report(args.output, report)

    # 5) Generate and save the heatmap if requested
    if args.plot:
        plot_matrix(matrix, args.plot)
        print(f"Heatmap saved to {args.plot}")

    # 6) Structured JSON export if requested
    if args.json_out:
        # collect single-path or multi-path alignments into a list
        if args.all_paths:
            alignments = trace_all_paths(
                matrix,
                seq1,
                seq2,
                match=args.match,
                mismatch=args.mismatch,
                gap=args.gap,
            )
        else:
            alignments = [
                single_traceback(
                    matrix,
                    seq1,
                    seq2,
                    match=args.match,
                    mismatch=args.mismatch,
                    gap=args.gap,
                )
            ]

        # make sure output folder exists
        os.makedirs(os.path.dirname(args.json_out) or ".", exist_ok=True)

        data_dict = create_output_dict(
            seq1, seq2, matrix, alignments, args.match, args.mismatch, args.gap
        )
        write_json(args.json_out, data_dict)
        print(f"Structured JSON written to {args.json_out}")


if __name__ == "__main__":
    main()
