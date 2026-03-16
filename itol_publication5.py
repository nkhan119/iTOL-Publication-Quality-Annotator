#!/usr/bin/env python3
"""
Publication-Quality iTOL Generator

Colored branches by phylum + Leptospira (red) and E. coli (blue) highlighting,
with special highlighting for specific E. coli K12/O157, Shigella, and
Bacillus subtilis sequences.

Usage:
  python3 itol_publication5.py step1 <csv> <tree> <outdir>
  python3 itol_publication5.py step2 <outdir>
  python3 itol_publication5.py step3 <outdir>
"""

import sys
import csv
import re
from pathlib import Path
from collections import defaultdict, OrderedDict


# ─────────────────────────────────────────────────────────────
# STAR GENERA — fixed highlight colors, override phylum coloring
# ─────────────────────────────────────────────────────────────

STARS = {
    "Leptospira":  "#E63946",  # vivid red
    "Escherichia": "#2196F3",  # strong blue
}


# ─────────────────────────────────────────────────────────────
# PHYLUM HINT TABLE
# Maps known phylum/class strings (lower-case, partial match) to a
# canonical phylum label. More-specific patterns must come before
# broader ones. The dynamic palette in step2 assigns actual colors.
# ─────────────────────────────────────────────────────────────

PHYLUM_HINTS = [
    # Spirochaetes — before Proteobacteria (Sphaerochaeta has spirochetes label)
    ("spirochaet",          "Spirochaetes"),
    ("sphaerochaet",        "Spirochaetes"),
    ("leptospir",           "Spirochaetes"),
    ("borreli",             "Spirochaetes"),
    ("treponem",            "Spirochaetes"),
    ("spirochet",           "Spirochaetes"),
    ("sphaerochaeta",       "Spirochaetes"),
    ("parasphaerochaeta",   "Spirochaetes"),
    ("salinispira",         "Spirochaetes"),

    # Chloroflexi — before Proteobacteria
    ("gns bacteria",        "Chloroflexi"),
    ("chloroflex",          "Chloroflexi"),
    ("chloroflexi",         "Chloroflexi"),
    ("roseiflexus",         "Chloroflexi"),
    ("bellilinea",          "Chloroflexi"),
    ("anaerolinea",         "Chloroflexi"),
    ("caldilinea",          "Chloroflexi"),

    # Proteobacteria — specific classes before generic
    ("gammaproteobacteri",  "Proteobacteria (γ)"),
    ("gamma-proteobacteri", "Proteobacteria (γ)"),
    ("g-proteobacteri",     "Proteobacteria (γ)"),
    ("enterobacteri",       "Proteobacteria (γ)"),
    ("e. coli",             "Proteobacteria (γ)"),
    ("cdc enteric group",   "Proteobacteria (γ)"),
    ("enteric group",       "Proteobacteria (γ)"),
    ("nih group",           "Proteobacteria (γ)"),
    ("yokenella",           "Proteobacteria (γ)"),
    ("atlantibacter",       "Proteobacteria (γ)"),
    ("enterobacter ",       "Proteobacteria (γ)"),
    ("kosakonia",           "Proteobacteria (γ)"),
    ("pseudescherichia",    "Proteobacteria (γ)"),
    ("betaproteobacteri",   "Proteobacteria (β)"),
    ("beta-proteobacteri",  "Proteobacteria (β)"),
    ("b-proteobacteri",     "Proteobacteria (β)"),
    ("alphaproteobacteri",  "Proteobacteria (α)"),
    ("alpha-proteobacteri", "Proteobacteria (α)"),
    ("a-proteobacteri",     "Proteobacteria (α)"),
    ("deltaproteobacteri",  "Proteobacteria (δ)"),
    ("delta-proteobacteri", "Proteobacteria (δ)"),
    ("d-proteobacteri",     "Proteobacteria (δ)"),
    ("desulfovibrio",       "Proteobacteria (δ)"),
    ("pseudodesulfovibrio", "Proteobacteria (δ)"),
    ("desulfobacteri",      "Proteobacteria (δ)"),
    ("desulfuromonas",      "Proteobacteria (δ)"),
    ("desulfomicrobium",    "Proteobacteria (δ)"),
    ("desulfonatronum",     "Proteobacteria (δ)"),
    ("oceanidesulfovibrio", "Proteobacteria (δ)"),
    ("desulfospira",        "Proteobacteria (δ)"),
    ("desulfosarcina",      "Proteobacteria (δ)"),
    ("desulfovulcanus",     "Proteobacteria (δ)"),
    ("halodesulfovibrio",   "Proteobacteria (δ)"),
    ("limisalsivibrio",     "Proteobacteria (δ)"),
    ("bilophila",           "Proteobacteria (δ)"),
    ("denitrovibrio",       "Proteobacteria (δ)"),
    ("malonomonas",         "Proteobacteria (δ)"),
    ("geobacter",           "Proteobacteria (δ)"),
    ("myxococc",            "Proteobacteria (δ)"),
    ("sorangium",           "Proteobacteria (δ)"),
    ("cystobacter",         "Proteobacteria (δ)"),
    ("chondromyces",        "Proteobacteria (δ)"),
    ("sandaracinus",        "Proteobacteria (δ)"),
    ("enhygromyxa",         "Proteobacteria (δ)"),
    ("corallococcus",       "Proteobacteria (δ)"),
    ("polyangium",          "Proteobacteria (δ)"),
    ("hyalangium",          "Proteobacteria (δ)"),
    ("pendulispora",        "Proteobacteria (δ)"),
    ("epsilonproteobacteri","Proteobacteria (ε)"),
    ("campylobacter",       "Proteobacteria (ε)"),
    ("helicobacter",        "Proteobacteria (ε)"),
    ("sulfurimonas",        "Proteobacteria (ε)"),
    ("sulfurovum",          "Proteobacteria (ε)"),
    ("nautilia",            "Proteobacteria (ε)"),
    ("proteobacteri",       "Proteobacteria"),

    # Firmicutes
    ("firmicute",           "Firmicutes"),
    ("bacilli",             "Firmicutes"),
    ("clostridia",          "Firmicutes"),
    ("lactococc",           "Firmicutes"),
    ("staphylococc",        "Firmicutes"),
    ("streptococc",         "Firmicutes"),
    ("anthrax bacterium",   "Firmicutes"),
    ("bacillus ",           "Firmicutes"),
    ("clostridium ",        "Firmicutes"),
    ("listeria ",           "Firmicutes"),
    ("enterococcus ",       "Firmicutes"),

    # Actinobacteria
    ("actinobacteri",       "Actinobacteria"),
    ("actinomycet",         "Actinobacteria"),
    ("high g+c gram",       "Actinobacteria"),
    ("high gc gram",        "Actinobacteria"),
    ("mycobacteri",         "Actinobacteria"),
    ("streptomyces ",       "Actinobacteria"),
    ("corynebacteri",       "Actinobacteria"),
    ("acidimicrobium",      "Actinobacteria"),
    ("kribbella",           "Actinobacteria"),

    # Other phyla
    ("cyanobacteri",        "Cyanobacteria"),
    ("fusobacteri",         "Fusobacteria"),
    ("ilyobacter",          "Fusobacteria"),
    ("psychrilyobacter",    "Fusobacteria"),
    ("bacteroidet",         "Bacteroidetes"),
    ("bacteroid",           "Bacteroidetes"),
    ("deinococcus",         "Deinococcus-Thermus"),
    ("thermus",             "Deinococcus-Thermus"),
    ("planctomycet",        "Planctomycetes"),
    ("planctopirus",        "Planctomycetes"),
    ("gimesia",             "Planctomycetes"),
    ("verrucomicrobi",      "Verrucomicrobia"),
    ("verrucomicrobium",    "Verrucomicrobia"),
    ("pedosphaera",         "Verrucomicrobia"),
    ("roseimicrobium",      "Verrucomicrobia"),
    ("chlamydi",            "Chlamydiae"),
    ("tenericute",          "Tenericutes"),
    ("mycoplasm",           "Tenericutes"),
    ("spiroplasm",          "Tenericutes"),
    ("archaea",             "Archaea"),
    ("eukaryo",             "Eukaryota"),
    ("deferribacter",       "Deferribacterota"),
    ("deferrivibrio",       "Deferribacterota"),
    ("seleniivibrio",       "Deferribacterota"),
    ("thermodesulfobium",   "Nitrospirae"),
    ("thermodesulfatator",  "Nitrospirae"),
    ("dethiosulfovibrio",   "Nitrospirae"),
    ("chloracidobacterium", "Acidobacteria"),
    ("pyrinomonas",         "Acidobacteria"),
    ("acetomicrobium",      "Synergistota"),
]


# ─────────────────────────────────────────────────────────────
# PHYLUM COLOR PALETTE — Nature/Cell journal style
# Assigned dynamically in step2 based on what appears in the data.
# ─────────────────────────────────────────────────────────────

PHYLUM_PALETTE = [
    "#FF9F1C",  # amber-orange
    "#4CAF50",  # green
    "#00BCD4",  # cyan
    "#9C27B0",  # purple
    "#FF5722",  # deep orange
    "#673AB7",  # deep purple
    "#CDDC39",  # lime
    "#009688",  # teal-green
    "#E91E63",  # pink
    "#8BC34A",  # light green
    "#0288D1",  # light blue
    "#795548",  # brown
    "#607D8B",  # blue-grey
    "#FFC107",  # amber
    "#CE93D8",  # lavender
    "#78909C",  # steel grey
    "#6D4C41",  # dark brown
    "#F06292",  # light pink
    "#80CBC4",  # pale teal
    "#A5D6A7",  # pale green
    "#FF8F00",  # deep amber
    "#4DB6AC",  # medium teal
    "#AED581",  # yellow-green
    "#F48FB1",  # soft pink
]

PHYLUM_FALLBACK = "#AAAAAA"


# ─────────────────────────────────────────────────────────────
# SPECIAL SEQUENCES — override all other coloring
# Format: tree_id -> (display_label, hex_color, symbol_shape)
# iTOL symbol shapes: 1=circle, 2=square, 3=triangle, 4=star, 5=diamond
# ─────────────────────────────────────────────────────────────

SPECIAL = {
    "GCF_000005845.2|NP_415951.1":    ("E. coli K-12 MG1655",   "#00897B", 2),  # teal   / square
    "GCF_000008865.2|NP_310064.3":    ("E. coli O157:H7 Sakai", "#FF6F00", 2),  # amber  / square
    "GCF_002290485.1|WP_000429155.1": ("Shigella boydii",        "#AD1457", 5),  # rose   / diamond
    "GCF_002950395.1|WP_000429142.1": ("Shigella sonnei",        "#AD1457", 5),  # rose   / diamond
    "GCF_022354085.1|WP_021577428.1": ("Shigella dysenteriae",   "#AD1457", 5),  # rose   / diamond
    "O34527|CYMR_BACSU_C":            ("Bacillus subtilis",      "#7B1FA2", 3),  # violet / triangle
    "P0A9F3|CYSB_ECOLI_C":            ("E. coli CysB",           "#0277BD", 1),  # blue   / circle
}

# Label for each special color (used in legends)
SPECIAL_COLOR_LABEL = {
    "#FF6F00": "E. coli O157 (special)",
    "#AD1457": "Shigella (special)",
    "#00897B": "E. coli K-12 MG1655 (special)",
    "#7B1FA2": "B. subtilis (special)",
    "#0277BD": "E. coli CysB (special)",
}


# ─────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────

def get_genus(name: str) -> str:
    """Return the first word of an organism name string."""
    return name.strip().split()[0]


def extract_phylum(name: str) -> str:
    """
    Best-effort phylum classification from an organism name string.

    Priority:
      1. If the genus is a STARS entry, return the genus name so it keeps
         its own highlight color in step3.
      2. Scan all parenthetical groups (innermost first) against PHYLUM_HINTS.
      3. Scan the full name string against PHYLUM_HINTS.
      4. Return "Unknown".
    """
    genus = get_genus(name)
    if genus in STARS:
        return genus

    parens = re.findall(r'\(([^)]+)\)', name)
    candidates = list(reversed(parens)) + [name]

    for candidate in candidates:
        cl = candidate.lower()
        for pattern, phylum in PHYLUM_HINTS:
            if pattern in cl:
                return phylum

    return "Unknown"


def assign_phylum_colors(phylum_counts: dict) -> dict:
    """
    Assign a distinct color to each phylum present in the data.
    Phyla are sorted by descending count so the largest groups get the
    most prominent palette slots. STARS genera keep their fixed colors.
    """
    color_map = {}
    palette_idx = 0
    for phylum, _ in sorted(phylum_counts.items(), key=lambda x: -x[1]):
        if phylum in STARS:
            color_map[phylum] = STARS[phylum]
        elif phylum == "Unknown":
            color_map[phylum] = PHYLUM_FALLBACK
        else:
            color_map[phylum] = PHYLUM_PALETTE[palette_idx % len(PHYLUM_PALETTE)]
            palette_idx += 1
    return color_map


def load_tree_ids(tree_file):
    """
    Extract all leaf IDs from the tree file.

    Returns:
      ids      — raw set of ID strings as they appear in the tree
      norm_map — mapping from bare CSV-style IDs to full tree IDs

    Supported formats:
      GCF_000005845.2|NP_415951.1       (NCBI RefSeq)
      sp|O34527|CYMR_BACSU_CymR_OS_...  (UniProtKB with sp| prefix)
      O34527|CYMR_BACSU_C               (bare UniProt)
    """
    pat = re.compile(
        r'('
        r'GCF_[0-9]+\.[0-9]+\|[A-Z0-9_]+\.[0-9]+'
        r'|sp\|[A-Z][A-Z0-9]{5}\|[A-Za-z0-9_]+'
        r'|[A-Z][A-Z0-9]{5}\|[A-Z0-9_]+'
        r')'
    )
    ids = set()
    with open(tree_file) as f:
        for line in f:
            ids.update(pat.findall(line))

    norm_map = {}
    sp_re = re.compile(r'^sp\|([A-Z][A-Z0-9]{5})\|([A-Z0-9_]+)')
    for tid in ids:
        m = sp_re.match(tid)
        if m:
            bare = f"{m.group(1)}|{m.group(2).rstrip('_')}"
            norm_map[bare] = tid
            continue
        n = tid.replace('OS_', 'OS=')
        if n != tid:
            norm_map[n] = tid

    return ids, norm_map


def safe_load_colors(groups_file):
    """Read the phylum→color mapping written by step2."""
    colors = {}
    with open(groups_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or '|' not in line:
                continue
            group, rest = line.split('|', 1)
            color = rest.strip().split()[0]
            if re.match(r'^#[0-9A-Fa-f]{6}$', color):
                colors[group.strip()] = color
    return colors


def _build_special_legend(found_special_dict):
    """
    Return (shapes, colors, labels) for legend blocks, deduplicated by color.
    Accepts found_special dict {tree_id: (label, color, shape)}.
    """
    seen = {}
    for _, col, sym in found_special_dict.values():
        if col not in seen:
            seen[col] = (str(sym), SPECIAL_COLOR_LABEL.get(col, col))
    shapes = [v[0] for v in seen.values()]
    colors = list(seen.keys())
    labels = [v[1] for v in seen.values()]
    return shapes, colors, labels


def _bare_id(tid):
    """Normalise a sp|ACC|ENTRY tree ID to bare ACC|ENTRY form."""
    m = re.match(r'^sp\|([A-Z][A-Z0-9]{5})\|([A-Z0-9_]+)', tid)
    return f"{m.group(1)}|{m.group(2).rstrip('_')}" if m else tid


# ─────────────────────────────────────────────────────────────
# STEP 1 — Filter CSV and generate labels
# ─────────────────────────────────────────────────────────────

def step1(csv_file, tree_file, output_dir):
    csv_file   = Path(csv_file)
    tree_file  = Path(tree_file)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 60)
    print("  STEP 1: Filter & Labels")
    print("=" * 60)

    tree_ids, norm_map = load_tree_ids(tree_file)
    print(f"\n  Tree IDs  : {len(tree_ids)}")

    matched, skipped = [], []
    with open(csv_file, newline='') as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            if len(row) < 4:
                continue
            cid = row[0]
            if cid in tree_ids:
                matched.append([cid] + row[1:4])
            elif cid in norm_map:
                matched.append([norm_map[cid]] + row[1:4])
            else:
                skipped.append(row)

    print(f"  Matched   : {len(matched)}")
    print(f"  Skipped   : {len(skipped)}")

    # Write filtered CSV
    fcsv = output_dir / "filtered.csv"
    with open(fcsv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(matched)

    # Write skipped IDs with diagnostic hints
    if skipped:
        gcf_tree = defaultdict(set)
        for tid in tree_ids:
            if tid.startswith('GCF_'):
                gcf_tree[tid.split('|')[0]].add(tid)
        with open(output_dir / "skipped_ids.txt", 'w') as f:
            f.write(f"SKIPPED ({len(skipped)})\n{'─' * 50}\n\n")
            for row in skipped:
                f.write(f"• {row[0]}\n  {row[2]}\n")
                acc = row[0].split('|')[0]
                if acc in gcf_tree:
                    for alt in sorted(gcf_tree[acc]):
                        f.write(f"  ⚠ Same genome → {alt}\n")
                f.write("\n")

    # Write iTOL labels file
    lf = output_dir / "itol_labels.txt"
    with open(lf, 'w') as f:
        f.write("LABELS\nSEPARATOR TAB\nDATA\n")
        for row in matched:
            label = row[0].replace('|', '_') + '_' + row[3]
            f.write(f"{row[0]}\t{label}\n")

    # Report which special IDs were found
    matched_bare = {_bare_id(r[0]) for r in matched} | {r[0] for r in matched}
    print(f"\n  Special sequences found in filtered set:")
    for sid, (lbl, col, _) in SPECIAL.items():
        status = "✓" if sid in matched_bare else "✗ MISSING"
        print(f"    {status}  {sid}  ({lbl})  {col}")

    print(f"\n  ✓ itol_labels.txt ({len(matched)} labels)")
    print(f"\n  Next: python3 {Path(sys.argv[0]).name} step2 {output_dir}")


# ─────────────────────────────────────────────────────────────
# STEP 2 — Census, classify, and assign colors
# ─────────────────────────────────────────────────────────────

def step2(output_dir):
    output_dir = Path(output_dir)
    fcsv = output_dir / "filtered.csv"
    if not fcsv.exists():
        print("❌ Run step1 first")
        sys.exit(1)

    print("\n" + "=" * 60)
    print("  STEP 2: Count → Classify → Assign Colors")
    print("=" * 60)

    counts     = defaultdict(int)
    raw_groups = defaultdict(list)
    unknown_examples = []

    with open(fcsv, newline='') as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            if len(row) >= 3:
                name   = row[2]
                phylum = extract_phylum(name)
                counts[phylum] += 1
                if len(raw_groups[phylum]) < 3:
                    raw_groups[phylum].append(name)
                if phylum == "Unknown":
                    unknown_examples.append(name)

    sorted_phylums = sorted(counts.items(), key=lambda x: -x[1])
    total = sum(counts.values())
    color_map = assign_phylum_colors(counts)

    print(f"\n  Total sequences : {total}")
    print(f"  Distinct phyla  : {len(counts)}\n")
    print(f"  {'#':<4} {'Phylum':<36} {'N':>5}  {'%':>5}   {'Color'}    Note")
    print(f"  {'─' * 4} {'─' * 36} {'─' * 5}  {'─' * 5}   {'─' * 9}   {'─' * 20}")

    for i, (phylum, n) in enumerate(sorted_phylums, 1):
        c    = color_map[phylum]
        pct  = 100 * n / total if total else 0
        note = ""
        if phylum in STARS:
            note = "⭐ highlighted genus"
        elif phylum == "Unknown":
            note = "⚠ unrecognised — check PHYLUM_HINTS"
        print(f"  {i:<4} {phylum:<36} {n:>5}  {pct:>4.1f}%   {c}   {note}")

    print(f"\n  Example organism names per phylum (up to 3 each):")
    for phylum, _ in sorted_phylums:
        print(f"    [{phylum}]")
        for ex in raw_groups[phylum]:
            print(f"      • {ex}")

    if unknown_examples:
        print(f"\n  ⚠  {counts['Unknown']} sequences could not be classified.")
        print(f"     Add matching patterns to PHYLUM_HINTS to fix this.")
        print(f"     First few unrecognised names:")
        for ex in unknown_examples[:5]:
            print(f"       • {ex}")

    gf = output_dir / "color_groups.txt"
    with open(gf, 'w') as f:
        f.write("# Color groups — generated by step2\n")
        f.write("# Edit hex values here if desired; step3 reads this file.\n")
        f.write("# Format: phylum_label|#RRGGBB  (N seqs)\n#\n")
        for phylum, n in sorted_phylums:
            c   = color_map[phylum]
            tag = "HIGHLIGHTED " if phylum in STARS else ""
            f.write(f"{phylum}|{c}  {tag}({n} seqs)\n")

    print(f"\n  ✓ {gf}  ({len(sorted_phylums)} phyla)")
    print(f"\n  Special sequences (hardcoded, always override phylum color):")
    shape_names = {1: "circle", 2: "square", 3: "triangle", 4: "star", 5: "diamond"}
    for sid, (lbl, col, sym) in SPECIAL.items():
        print(f"    {col}  shape={shape_names.get(sym, sym):<8}  {lbl}  ({sid})")
    print(f"\n  Next: python3 {Path(sys.argv[0]).name} step3 {output_dir}")


# ─────────────────────────────────────────────────────────────
# STEP 3 — Generate all iTOL annotation files
# ─────────────────────────────────────────────────────────────

def step3(output_dir):
    output_dir = Path(output_dir)
    fcsv = output_dir / "filtered.csv"
    gf   = output_dir / "color_groups.txt"
    if not fcsv.exists():
        print("❌ Run step1 first")
        sys.exit(1)
    if not gf.exists():
        print("❌ Run step2 first")
        sys.exit(1)

    print("\n" + "=" * 60)
    print("  STEP 3: Generating All iTOL Files")
    print("=" * 60)

    colors = safe_load_colors(gf)
    rows = []
    with open(fcsv, newline='') as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            if len(row) >= 4:
                rows.append(row)

    # Groups present in data (for accurate legends)
    used = OrderedDict()
    for row in rows:
        g = extract_phylum(row[2])
        if g not in used:
            used[g] = colors.get(g, "#AAAAAA")

    lepto = colors.get("Leptospira",  "#C0392B")
    ecoli = colors.get("Escherichia", "#1565C0")
    lepto_ids = [r[0] for r in rows if get_genus(r[2]) == "Leptospira"]
    ecoli_ids = [r[0] for r in rows if get_genus(r[2]) == "Escherichia"]

    # Build special lookup: tree_id -> (label, color, shape)
    found_special = {}
    for r in rows:
        tid  = r[0]
        bare = _bare_id(tid)
        if bare in SPECIAL:
            found_special[tid] = SPECIAL[bare]
        elif tid in SPECIAL:
            found_special[tid] = SPECIAL[tid]

    sp_shapes, sp_colors, sp_labels = _build_special_legend(found_special)

    print(f"\n  Groups       : {len(used)}")
    print(f"  Leptospira   : {len(lepto_ids)} seqs  {lepto}")
    print(f"  Escherichia  : {len(ecoli_ids)} seqs  {ecoli}")
    print(f"  Special seqs : {len(found_special)} / {len(SPECIAL)} found in data")
    for sid, (lbl, col, _) in found_special.items():
        print(f"    ★  {col}  {lbl}  ({sid})")
    print()

    # [1] Branch colors
    # Priority: special > Leptospira > Escherichia > phylum
    with open(output_dir / "itol_branch_colors.txt", 'w') as f:
        f.write("TREE_COLORS\nSEPARATOR TAB\nDATA\n")
        for row in rows:
            nid   = row[0]
            genus = get_genus(row[2])
            c     = colors.get(extract_phylum(row[2]), "#AAAAAA")
            if nid in found_special:
                _, sc, _ = found_special[nid]
                f.write(f"{nid}\tbranch\t{sc}\tnormal\t6\n")
                f.write(f"{nid}\tlabel\t{sc}\tbold-italic\t1.4\n")
            elif genus == "Leptospira":
                f.write(f"{nid}\tbranch\t{lepto}\tnormal\t5\n")
                f.write(f"{nid}\tlabel\t{lepto}\tbold-italic\t1\n")
            elif genus == "Escherichia":
                f.write(f"{nid}\tbranch\t{ecoli}\tnormal\t5\n")
                f.write(f"{nid}\tlabel\t{ecoli}\tbold-italic\t1\n")
            else:
                f.write(f"{nid}\tbranch\t{c}\tnormal\t2\n")
                f.write(f"{nid}\tlabel\t{c}\tnormal\t1\n")
    print(f"  ✓ [1] itol_branch_colors.txt")

    # [2] Color strip
    all_legend_labels = list(used.keys()) + sp_labels
    all_legend_colors = list(used.values()) + sp_colors
    shapes_str = "\t".join(["1"] * len(all_legend_labels))

    with open(output_dir / "itol_colorstrip.txt", 'w') as f:
        f.write("DATASET_COLORSTRIP\nSEPARATOR TAB\n")
        f.write("DATASET_LABEL\tTaxonomy\nCOLOR\t#000000\n")
        f.write("STRIP_WIDTH\t35\nMARGIN\t10\nBORDER_WIDTH\t0\nSHOW_INTERNAL\t0\n\n")
        f.write("LEGEND_TITLE\tPhylum / Genus\n")
        f.write("LEGEND_SHAPES\t" + shapes_str + "\n")
        f.write("LEGEND_COLORS\t" + "\t".join(all_legend_colors) + "\n")
        f.write("LEGEND_LABELS\t" + "\t".join(all_legend_labels) + "\n\n")
        f.write("DATA\n")
        for row in rows:
            nid = row[0]
            if nid in found_special:
                _, sc, _ = found_special[nid]
                f.write(f"{nid}\t{sc}\t{found_special[nid][0]}\n")
            else:
                c   = colors.get(extract_phylum(row[2]), "#AAAAAA")
                grp = extract_phylum(row[2])
                f.write(f"{nid}\t{c}\t{grp}\n")
    print(f"  ✓ [2] itol_colorstrip.txt")

    # [3] Label backgrounds
    with open(output_dir / "itol_label_background.txt", 'w') as f:
        f.write("DATASET_STYLE\nSEPARATOR TAB\n")
        f.write("DATASET_LABEL\tLabel Highlights\nCOLOR\t#000000\nDATA\n")
        for row in rows:
            nid   = row[0]
            genus = get_genus(row[2])
            if nid in found_special:
                _, sc, _ = found_special[nid]
                f.write(f"{nid}\tlabel\tnode\t{sc}\t#FFFFFF\tbold-italic\t1.3\n")
            elif genus == "Leptospira":
                f.write(f"{nid}\tlabel\tnode\t{lepto}\t#FFFFFF\tbold-italic\t1.2\n")
            elif genus == "Escherichia":
                f.write(f"{nid}\tlabel\tnode\t{ecoli}\t#FFFFFF\tbold-italic\t1.2\n")
    print(f"  ✓ [3] itol_label_background.txt")

    # [4] Tip symbols
    with open(output_dir / "itol_symbols.txt", 'w') as f:
        f.write("DATASET_SYMBOL\nSEPARATOR TAB\n")
        f.write("DATASET_LABEL\tKey Taxa\nCOLOR\t#000000\n")
        leg_shapes = ["1", "4"] + sp_shapes
        leg_colors = [lepto, ecoli] + sp_colors
        leg_labels = ["Leptospira", "Escherichia"] + sp_labels
        f.write("LEGEND_TITLE\tKey Taxa\n")
        f.write("LEGEND_SHAPES\t" + "\t".join(leg_shapes) + "\n")
        f.write("LEGEND_COLORS\t" + "\t".join(leg_colors) + "\n")
        f.write("LEGEND_LABELS\t" + "\t".join(leg_labels) + "\n\n")
        f.write("DATA\n")
        for nid in lepto_ids:
            if nid not in found_special:
                f.write(f"{nid}\t1\t12\t{lepto}\t1\t0\n")
        for nid in ecoli_ids:
            if nid not in found_special:
                f.write(f"{nid}\t4\t12\t{ecoli}\t1\t0\n")
        for sid, (lbl, col, sym) in found_special.items():
            f.write(f"{sid}\t{sym}\t16\t{col}\t1\t0\n")
    print(f"  ✓ [4] itol_symbols.txt")

    # [5] Clade highlight strip
    with open(output_dir / "itol_clade_highlight.txt", 'w') as f:
        f.write("DATASET_COLORSTRIP\nSEPARATOR TAB\n")
        f.write("DATASET_LABEL\tKey Clades\nCOLOR\t#000000\n")
        f.write("STRIP_WIDTH\t25\nMARGIN\t2\nBORDER_WIDTH\t0\nSHOW_INTERNAL\t0\n\n")
        leg_shapes = ["1", "1"] + ["1"] * len(sp_colors)
        leg_colors = [lepto, ecoli] + sp_colors
        leg_labels = ["Leptospira", "Escherichia"] + sp_labels
        f.write("LEGEND_TITLE\tHighlighted Clades\n")
        f.write("LEGEND_SHAPES\t" + "\t".join(leg_shapes) + "\n")
        f.write("LEGEND_COLORS\t" + "\t".join(leg_colors) + "\n")
        f.write("LEGEND_LABELS\t" + "\t".join(leg_labels) + "\n\n")
        f.write("DATA\n")
        for nid in lepto_ids:
            f.write(f"{nid}\t{lepto}\tLeptospira\n")
        for nid in ecoli_ids:
            if nid not in found_special:
                f.write(f"{nid}\t{ecoli}\tEscherichia\n")
        for sid, (lbl, col, _) in found_special.items():
            f.write(f"{sid}\t{col}\t{lbl}\n")
    print(f"  ✓ [5] itol_clade_highlight.txt")

    # [6] Binary highlight
    special_color_cols = list(dict.fromkeys(col for _, (_, col, _) in SPECIAL.items()))
    field_labels = ["Leptospira", "Escherichia"] + \
                   [SPECIAL_COLOR_LABEL.get(c, c) for c in special_color_cols]
    field_colors = [lepto, ecoli] + special_color_cols
    shape_map    = {c: str(sym) for _, (_, c, sym) in SPECIAL.items()}
    field_shapes = ["1", "2"] + [shape_map.get(c, "1") for c in special_color_cols]

    with open(output_dir / "itol_binary_highlight.txt", 'w') as f:
        f.write("DATASET_BINARY\nSEPARATOR TAB\n")
        f.write("DATASET_LABEL\tKey Taxa Presence\nCOLOR\t#000000\n")
        f.write("FIELD_LABELS\t"  + "\t".join(field_labels) + "\n")
        f.write("FIELD_COLORS\t"  + "\t".join(field_colors) + "\n")
        f.write("FIELD_SHAPES\t"  + "\t".join(field_shapes) + "\n")
        f.write("MARGIN\t5\nSTRIP_WIDTH\t20\n\n")
        f.write("LEGEND_TITLE\tKey Taxa\n")
        f.write("LEGEND_SHAPES\t" + "\t".join(field_shapes) + "\n")
        f.write("LEGEND_COLORS\t" + "\t".join(field_colors) + "\n")
        f.write("LEGEND_LABELS\t" + "\t".join(field_labels) + "\n\n")
        f.write("DATA\n")
        lepto_set    = set(lepto_ids)
        ecoli_set    = set(ecoli_ids)
        color_to_ids = defaultdict(set)
        for sid, (_, col, _) in SPECIAL.items():
            color_to_ids[col].add(sid)
        for row in rows:
            nid  = row[0]
            vals = [
                "1" if nid in lepto_set else "0",
                "1" if nid in ecoli_set else "0",
            ]
            vals += ["1" if nid in color_to_ids[c] else "0" for c in special_color_cols]
            if "1" in vals:
                f.write(nid + "\t" + "\t".join(vals) + "\n")
    print(f"  ✓ [6] itol_binary_highlight.txt")

    # [7] Special sequences strip
    with open(output_dir / "itol_special_seqs.txt", 'w') as f:
        f.write("DATASET_COLORSTRIP\nSEPARATOR TAB\n")
        f.write("DATASET_LABEL\tSpecial Sequences\nCOLOR\t#000000\n")
        f.write("STRIP_WIDTH\t30\nMARGIN\t5\nBORDER_WIDTH\t1\nSHOW_INTERNAL\t0\n\n")
        f.write("LEGEND_TITLE\tSpecial Sequences\n")
        f.write("LEGEND_SHAPES\t" + "\t".join(["1"] * len(sp_colors)) + "\n")
        f.write("LEGEND_COLORS\t" + "\t".join(sp_colors) + "\n")
        f.write("LEGEND_LABELS\t" + "\t".join(sp_labels) + "\n\n")
        f.write("DATA\n")
        for sid, (lbl, col, _) in found_special.items():
            f.write(f"{sid}\t{col}\t{lbl}\n")
    print(f"  ✓ [7] itol_special_seqs.txt")

    # Upload guide
    shape_names = {1: "circle", 2: "square", 3: "triangle", 4: "star", 5: "diamond"}
    with open(output_dir / "UPLOAD_ORDER.txt", 'w') as f:
        f.write("=" * 65 + "\n")
        f.write("  iTOL Upload Guide\n")
        f.write("=" * 65 + "\n\n")
        f.write("1. Go to https://itol.embl.de/upload.cgi\n")
        f.write("   Upload your .treefile\n\n")
        f.write("2. Drag files onto the tree in this order:\n\n")
        f.write("   [1] itol_labels.txt\n")
        f.write("   [2] itol_branch_colors.txt\n")
        f.write("   [3] itol_colorstrip.txt\n")
        f.write("   [4] itol_label_background.txt\n")
        f.write("   [5] itol_symbols.txt\n")
        f.write("   [6] itol_clade_highlight.txt\n")
        f.write("   [7] itol_special_seqs.txt\n")
        f.write("   [8] itol_binary_highlight.txt  (optional)\n\n")
        f.write("3. Recommended settings:\n")
        f.write("   Display: Circular | Labels: 7-9pt aligned\n")
        f.write("   Bootstrap: show >= 70 only | Branch width: 1\n\n")
        f.write("Special sequences:\n")
        for sid, (lbl, col, sym) in SPECIAL.items():
            status = "✓" if sid in found_special else "✗ not in data"
            f.write(f"  {status}  {col}  {shape_names.get(sym, sym):<8}  {lbl}  ({sid})\n")
        f.write(f"\nColors:\n")
        f.write(f"  Leptospira   {lepto}  ({len(lepto_ids)} seqs)\n")
        f.write(f"  Escherichia  {ecoli}  ({len(ecoli_ids)} seqs)\n")
        for col, lbl in SPECIAL_COLOR_LABEL.items():
            f.write(f"  {lbl:<28}{col}\n")

    print(f"\n  ✓ UPLOAD_ORDER.txt")
    print("\n" + "=" * 60)
    print("  ALL FILES READY")
    print("=" * 60)
    print(f"""
  Upload order:
    1. itol_labels.txt
    2. itol_branch_colors.txt     all branches colored
    3. itol_colorstrip.txt        outer taxonomy ring
    4. itol_label_background.txt  colored label backgrounds
    5. itol_symbols.txt           tip markers
    6. itol_clade_highlight.txt   clade strips
    7. itol_special_seqs.txt      special sequences strip
    8. itol_binary_highlight.txt  (optional)

  Leptospira   → {lepto}  bold x5  ({len(lepto_ids)} seqs)
  Escherichia  → {ecoli}  bold x5  ({len(ecoli_ids)} seqs)
""")


# ─────────────────────────────────────────────────────────────
# ENTRY POINT
# ─────────────────────────────────────────────────────────────

def main():
    cmds = {
        "step1": (step1, 4, "step1 <csv> <tree> <outdir>"),
        "step2": (step2, 2, "step2 <outdir>"),
        "step3": (step3, 2, "step3 <outdir>"),
    }
    if len(sys.argv) < 2 or sys.argv[1] not in cmds:
        print(__doc__)
        sys.exit(1)
    fn, nargs, usage = cmds[sys.argv[1]]
    if len(sys.argv) < nargs + 1:
        print(f"Usage: python3 {Path(sys.argv[0]).name} {usage}")
        sys.exit(1)
    fn(*sys.argv[2:nargs + 1])


if __name__ == "__main__":
    main()
