import argparse
import os
from typing import Optional

from src.aligner.plot import plot_matrix
from src.aligner.core import build_score_matrix, traceback
from src.aligner.io import format_report, read_manual, read_fasta, write_report


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
        "--plot",
        type=str,
        default=None,
        help="Optional filename for PNG heatmap output",
    )
    return parser.parse_args(args)


def main():
    """
    Entry point for console script. Parses args, loads sequences, and prepares for alignment.
    """
    args = parse_args()
    # Default alphabet is DNA
    alphabet = "dna"

    # Load sequences from manual input or FASTA files
    if args.manual:
        seq1, seq2 = read_manual(alphabet)
    else:
        fasta1, fasta2 = args.input  # type: ignore
        seqs1 = read_fasta(fasta1, alphabet)
        seqs2 = read_fasta(fasta2, alphabet)
        if len(seqs1) != 1 or len(seqs2) != 1:
            raise ValueError("Each FASTA file must contain exactly one record")
        seq1, seq2 = seqs1[0], seqs2[0]

    # For now, just print loaded sequences; next we'll run alignment
    print("Loaded sequences:")
    print(f"  {seq1.id}: {seq1.sequence}")
    print(f"  {seq2.id}: {seq2.sequence}")

    # Run alignment
    matrix = build_score_matrix(seq1, seq2, args.match, args.mismatch, args.gap)
    aln1, aln2 = traceback(matrix, seq1, seq2, args.match, args.mismatch, args.gap)

    # Display alignment
    print("\nOptimal Alignment:")
    print(aln1)
    print(aln2)

    # Format report
    report = format_report(seq1, seq2, aln1, aln2, args.match, args.mismatch, args.gap)

    # Prepare output filename
    # Ensure reports directory exists
    reports_dir = "reports"
    os.makedirs(reports_dir, exist_ok=True)
    if args.output:
        output_file = (
            os.path.join(reports_dir, args.output)
            if not args.output.startswith(reports_dir)
            else args.output
        )
    else:
        # default name: seq1_seq2_m<M>_mm<MM>_g<gap>.txt
        fname = f"{seq1.id}_{seq2.id}_m{args.match}_mm{args.mismatch}_g{args.gap}.txt"
        output_file = os.path.join(reports_dir, fname)

    # Print and write report
    print(report)
    write_report(output_file, report)
    print(f"Report written to {output_file}")

    # Handle plotting
    if args.plot:
        plots_dir = "plots"
        os.makedirs(plots_dir, exist_ok=True)
        plot_file = (
            args.plot
            if os.path.isabs(args.plot) or args.plot.startswith(plots_dir)
            else os.path.join(plots_dir, args.plot)
        )
        plot_matrix(matrix, plot_file)
        print(f"Heatmap saved to {plot_file}")


if __name__ == "__main__":
    main()
