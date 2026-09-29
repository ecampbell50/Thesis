# ======================================================================
# Thesis Chapter 5, Figures S5.6, S5.7
# QC table of all proviruses; flags ICE/composite over-calls.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/Chapter3_Coevolution/04_genomad/check_phages.sh
# ======================================================================
cat > qc_proviruses.py <<'EOF'
#!/usr/bin/env python3
"""Scan all geNomad outputs in CWD -> one QC table of every virus/provirus.
Flags likely ICE/composite over-calls (like the 212 kbp one) vs real prophages.
Run from the dir containing the *_genomad/ folders:  python qc_proviruses.py"""
import glob, os, re, csv, sys

def read_tsv(p):
    with open(p) as fh:
        rd = csv.reader(fh, delimiter="\t")
        h = next(rd, None)
        return (h or []), [r for r in rd]

def idx(h, name):            # column index by header name, or None
    return h.index(name) if name in h else None

rows = []
sfiles = sorted(glob.glob("*_genomad/*_summary/*_virus_summary.tsv"))
sys.stderr.write(f"found {len(sfiles)} virus_summary.tsv files\n")

for sf in sfiles:
    genome = os.path.basename(sf).replace("_virus_summary.tsv", "")
    gf = sf.replace("_virus_summary.tsv", "_virus_genes.tsv")

    # gene-level counts per provirus (from *_virus_genes.tsv)
    gs = {}
    if os.path.exists(gf):
        gh, grows = read_tsv(gf)
        gi = {c: idx(gh, c) for c in
              ["gene","marker","virus_hallmark","annotation_conjscan","annotation_amr"]}
        for r in grows:
            if gi["gene"] is None or gi["gene"] >= len(r): continue
            seq = re.sub(r"_\d+$", "", r[gi["gene"]])
            d = gs.setdefault(seq, dict(n=0, vv=0, pv=0, hm=0, conj=0, amr=0))
            d["n"] += 1
            mk = r[gi["marker"]] if gi["marker"] is not None and gi["marker"] < len(r) else "NA"
            if mk not in ("NA", ""):
                suf = mk.rsplit(".", 1)[-1]
                if suf == "VV": d["vv"] += 1
                elif suf == "PV": d["pv"] += 1
            def g(col):
                j = gi[col]; return r[j] if j is not None and j < len(r) else "NA"
            if g("virus_hallmark") == "1": d["hm"] += 1
            if g("annotation_conjscan") not in ("NA", ""): d["conj"] += 1
            if g("annotation_amr") not in ("NA", ""): d["amr"] += 1

    # summary rows
    sh, srows = read_tsv(sf)
    si = {c: idx(sh, c) for c in
          ["seq_name","length","topology","n_genes","n_hallmarks","virus_score","taxonomy"]}
    for r in srows:
        def s(col):
            j = si[col]; return r[j] if j is not None and j < len(r) else "NA"
        seq = s("seq_name")
        L = int(s("length")) if s("length").isdigit() else 0
        c = gs.get(seq, dict(n=0, vv=0, pv=0, hm=0, conj=0, amr=0))
        nh = int(s("n_hallmarks")) if s("n_hallmarks").isdigit() else c["hm"]
        rows.append(dict(genome=genome, seq_name=seq, length_bp=L, topology=s("topology"),
                         n_genes=s("n_genes"), n_hallmarks=nh, virus_score=s("virus_score"),
                         taxonomy=s("taxonomy"), n_conjscan=c["conj"], n_marker_VV=c["vv"],
                         n_marker_PV=c["pv"], n_amr=c["amr"],
                         genes_per_kbp=round(c["n"]/(L/1000), 2) if L else 0,
                         hall_per_100kb=round(nh/(L/1e5), 1) if L else 0))

for d in rows:
    f = []
    if d["length_bp"] >= 80000: f.append("LARGE")
    if d["n_conjscan"] >= 2: f.append("ICE_conj")
    if d["n_marker_PV"] > d["n_marker_VV"] and d["n_marker_VV"] >= 0: f.append("plasmid>=virus")
    if d["length_bp"] >= 80000 and d["hall_per_100kb"] < 12: f.append("sparse_hallmarks")
    d["flag"] = ";".join(f) or "ok"

cols = ["genome","seq_name","length_bp","topology","n_genes","n_hallmarks","virus_score",
        "taxonomy","n_conjscan","n_marker_VV","n_marker_PV","n_amr","genes_per_kbp",
        "hall_per_100kb","flag"]
with open("all_proviruses_QC.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader()
    for d in sorted(rows, key=lambda x: -x["length_bp"]): w.writerow(d)

print(f"\nwrote all_proviruses_QC.tsv  ({len(rows)} viruses / {len(sfiles)} genomes)\n")
print("=== 30 longest (inspect these) ===")
print(f"{'genome':15}{'len_kb':>7}{'genes':>6}{'hall':>5}{'VV':>4}{'PV':>4}{'conj':>5}  flag")
for d in sorted(rows, key=lambda x: -x["length_bp"])[:30]:
    print(f"{d['genome']:15}{d['length_bp']/1000:>7.1f}{str(d['n_genes']):>6}"
          f"{str(d['n_hallmarks']):>5}{d['n_marker_VV']:>4}{d['n_marker_PV']:>4}"
          f"{d['n_conjscan']:>5}  {d['flag']}")
EOF
python qc_proviruses.py
