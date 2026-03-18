# iTOL Publication-Quality Annotator

A Python tool to generate publication-ready annotation files for the
[Interactive Tree of Life (iTOL)](https://itol.embl.de/) viewer.
Automatically colors phylogenetic tree branches by phylum, with special
highlighting for key taxa such as *Leptospira*, *Escherichia*, *Shigella*,
and *Bacillus subtilis*.

---

## Features

- Classifies sequences by phylum using a curated hint table covering
  Proteobacteria (α/β/γ/δ/ε), Firmicutes, Actinobacteria, Spirochaetes,
  Chloroflexi, and many more
- Dynamically assigns a distinct color to each phylum present in the data
- Highlights star genera (*Leptospira*, *Escherichia*) with bold colored branches
- Supports hardcoded special sequences with unique colors and tip symbols
- Outputs seven ready-to-upload iTOL annotation files in one run

---

## Requirements

- Python 3.7+
- No external dependencies (standard library only)

---

## Usage

The tool runs in three sequential steps.

**Step 1** — Filter the input CSV against the tree and generate label files:

```bash
python3 itol_publication5.py step1 <csv> <tree> <outdir>
```

**Step 2** — Count sequences per phylum, assign colors, and write the color group file:

```bash
python3 itol_publication5.py step2 <outdir>
```

**Step 3** — Generate all seven iTOL annotation files:

```bash
python3 itol_publication5.py step3 <outdir>
```

---

## Input Files

| File | Description |
|------|-------------|
| `<csv>` | Sequence metadata CSV. Expected columns: ID, ?, organism name, label |
| `<tree>` | Newick tree file (e.g. `.treefile` from IQ-TREE) |
| `<outdir>` | Output directory (created if it does not exist) |

Supported sequence ID formats in the tree:

- `GCF_000005845.2|NP_415951.1` (NCBI RefSeq)
- `sp|O34527|CYMR_BACSU` (UniProtKB/Swiss-Prot with `sp|` prefix)
- `O34527|CYMR_BACSU` (bare UniProt accession)

---

## Output Files

Upload these to iTOL in the order listed for best results.

| # | File | Description |
|---|------|-------------|
| 1 | `itol_labels.txt` | Display labels for all leaves |
| 2 | `itol_branch_colors.txt` | Branch colors by phylum; key taxa bold and thicker |
| 3 | `itol_colorstrip.txt` | Outer taxonomy color ring |
| 4 | `itol_label_background.txt` | Colored label backgrounds for key taxa |
| 5 | `itol_symbols.txt` | Tip symbols (circles, stars, squares, etc.) |
| 6 | `itol_clade_highlight.txt` | Clade highlight strip |
| 7 | `itol_special_seqs.txt` | Dedicated strip for special sequences |
| 8 | `itol_binary_highlight.txt` | Optional binary presence/absence columns |

A plain-text `UPLOAD_ORDER.txt` is also written to the output directory
as a quick reference.

---

## Special Sequences

The following sequences receive individual colors and tip symbols that
override all phylum-level coloring:

| Sequence ID | Label | Color | Symbol |
|-------------|-------|-------|--------|
| `GCF_000005845.2\|NP_415951.1` | E. coli K-12 MG1655 | `#00897B` | square |
| `GCF_000008865.2\|NP_310064.3` | E. coli O157:H7 Sakai | `#FF6F00` | square |
| `GCF_002290485.1\|WP_000429155.1` | Shigella boydii | `#AD1457` | diamond |
| `GCF_002950395.1\|WP_000429142.1` | Shigella sonnei | `#AD1457` | diamond |
| `GCF_022354085.1\|WP_021577428.1` | Shigella dysenteriae | `#AD1457` | diamond |
| `O34527\|CYMR_BACSU_C` | Bacillus subtilis | `#7B1FA2` | triangle |
| `P0A9F3\|CYSB_ECOLI_C` | E. coli CysB | `#0277BD` | circle |

To add or modify special sequences, edit the `SPECIAL` dictionary near
the top of the script.

---

## Extending the Phylum Table

If sequences remain unclassified (`Unknown`), add a pattern to the
`PHYLUM_HINTS` list in the script:

```python
("your_genus_name", "Phylum Label"),
```

More specific patterns should be placed before broader ones. Step 2
reports all unclassified sequences to guide this process.

---

## iTOL Display Settings

- Display mode: Circular
- Labels: 7–9 pt, aligned
- Bootstrap: show values ≥ 70 only
- Branch width: 1 (overridden per-branch by the annotation files)

---

## Citation / Reuse

If you use or adapt this pipeline, please credit the author.

---

## Author
---

**Nadeem Khan, PhD**
Bioinformatician — INRS–Centre Armand-Frappier Santé-Biotechnologie, Laval, QC, Canada
nkhan119@uottawa.ca
[@nkhan119](https://github.com/nkhan119)
