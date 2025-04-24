# Needleman–Wunsch Aligner
bIOINFORMATICS @WUST

# Installation
git clone https://github.com/Karo555/needleman-wunsch.git
cd needleman-wunsch-aligner
pip install -e .[html,pdf]

# Usage
1. Basic DNA alignment (manual input)
needleman-wunsch --manual

You’ll be prompted to enter:
- Sequence 1 ID
- Sequence 2 ID
- Sequence 1
- Sequence 2

2. FASTA file input
needleman-wunsch --input data/seq1.fasta data/seq2.fasta --output reports/alignment.txt

data/seq1.fasta and data/seq2.fasta each contain exactly one record.
The text report is saved to reports/alignment.txt

3. Custom scoring: override match, mismatch, and gap scores
needleman-wunsch --input data/seq1.fasta data/seq2.fasta --match 2 --mismatch -1 --gap -2

4. Enumerate all optimal alignments
needleman-wunsch --input data/seq1.fasta data/seq2.fasta --all-paths --output reports/all_paths.txt

Lists every equally optimal alignment and writes them to reports/all_paths.txt.

5. Export raw DP matrix as CSV
needleman-wunsch --input data/seq1.fasta data/seq2.fasta --matrix-out reports/matrix.csv

The full (N+1)×(M+1) score matrix is saved in comma‐separated format for downstream analysis

6. Structured JSON output
needleman-wunsch --input data/seq1.fasta data/seq2.fasta --json reports/alignment.json

Produces a JSON file containing:

Original sequences and IDs
Scoring parameters
Full DP matrix
All alignments with per-path statistics

7. Heatmap visualization
needleman-wunsch --input data/seq1.fasta data/seq2.fasta --plot plots/heatmap.png

Renders the DP matrix as a heatmap PNG saved to plots/heatmap.png.

8. HTML summary report
needleman-wunsch --input data/seq1.fasta data/seq2.fasta --plot plots/heatmap.png --html reports/alignment.html

Generates a standalone HTML page embedding:

Parameters & sequences
All alignments (or just the single optimal one)
Inline statistics
The heatmap image (correctly linked)

Open it in your browser via:
python3 -m http.server 8000

then visit http://localhost:8000/reports/alignment.html

9. PDF summary report
needleman-wunsch --input data/seq1.fasta data/seq2.fasta --plot plots/heatmap.png --pdf reports/alignment.pdf

10. one shot
needleman-wunsch --input data/seq1.fasta data/seq2.fasta --all-paths --match 2 --mismatch -1 --gap -2 --output reports/alignment.txt --matrix-out reports/matrix.csv --json reports/alignment.json --plot plots/heatmap.png --html reports/alignment.html --pdf reports/alignment.pdf

reports/alignment.txt (text report)
reports/matrix.csv (raw DP matrix)
reports/alignment.json (structured JSON)
plots/heatmap.png (heatmap image)
reports/alignment.html (HTML summary)
reports/alignment.pdf (PDF summary)


## 📄 License
This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
