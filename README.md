# Needleman–Wunsch Aligner
bIOINFORMATICS @WUST

Python tool for global alignment of DNA and protein sequences using the classic Needleman–Wunsch algorithm. 

**Customizable Scoring**: Pick your own match, mismatch, and gap penalties to fit your study.  
**Multiple Paths**: Not just one “best” alignment—see all optimal ways your sequences can line up.  
**Stats on Demand**: Alignment length, percent identity, total gaps, and more.  
**Visual Insight**: Heatmap of the scoring matrix with the optimal traceback path overlaid.

Installation
1. Clone this repo:
   ```bash
   git clone https://github.com/Karo555/needleman-wunsch-aligner.git
   cd needleman-wunsch-aligner
   ```
2. Create & activate a virtual environment:
   ```bash
   python3 -m venv .venv   # or `python -m venv .venv`
   source .venv/bin/activate
   ```
3. Install the package plus dev-tools:
   ```bash
   pip install .
   pip install black pytest matplotlib
   ```

Basic Usage
### Align two FASTA files:
```bash
needleman-wunsch --input seq1.fasta seq2.fasta --match 1 --mismatch -1 --gap -2
```
This will print a text report to the console and save `alignment.txt` by default.

Interactive mode:
```bash
needleman-wunsch --manual
```
You'll be prompted to paste or type in your sequences and scoring parameters.

Save a heatmap:
```bash
needleman-wunsch --input seq1.fasta seq2.fasta --plot heatmap.png
```

Development & Testing
- **Auto-format** with [Black](https://github.com/psf/black):  
  ```bash
  black src tests
  ```
- **Run tests** with [pytest](https://pytest.org/):  
  ```bash
  pytest
  ```
- **CI** is set up via GitHub Actions to lint, type check, and test on every push.

## 📄 License
This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.


