"""Functional characterisation of the Core Drought Signature (CDS).

1. GO enrichment with annotations propagated to all ancestors (QuickGO), BH-FDR.
2. Arabidopsis orthologs (DIAMOND reciprocal best hits): enrichment of genes with experimental evidence
   (IMP/IGI/IDA/IPI) for 'response to water deprivation' (GO:0009414).
3. Transcription-factor families among CDS genes (Arabidopsis ortholog symbol + AHRD description).
4. Promoter cis-elements (1.5 kb upstream of the start codon, Ensembl Cla97 genome): ABRE, DRE/CRT,
   MYB, MYC2, W-box, G-box; enrichment in CDS-up / CDS-down promoters vs all tested genes.
Outputs in results/function/.
"""
import os
import gzip
import json
import re
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy.stats import fisher_exact, hypergeom

M = Path(os.environ.get("PROJECT_DIR", ".") + "")
R = M / "results"
F = R / "function"
F.mkdir(exist_ok=True)
REF = M / "ref"
CUGEN = Path(os.environ.get("PROJECT_DIR", ".") + "/reference")


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


# ---------------- GO with propagation ----------------
def go_ancestors(ids):
    cache = REF / "go_ancestors.json"
    known = json.loads(cache.read_text()) if cache.exists() else {}
    todo = [i for i in ids if i not in known]
    for k in range(0, len(todo), 50):
        chunk = todo[k:k + 50]
        for attempt in range(3):
            try:
                r = requests.get("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"
                                 + ",".join(chunk) + "/ancestors", params={"relations": "is_a,part_of"},
                                 headers={"Accept": "application/json"}, timeout=120)
                r.raise_for_status()
                got = {t["id"]: t.get("ancestors", []) for t in r.json().get("results", [])}
                for i in chunk:
                    known[i] = got.get(i, [i])
                break
            except requests.RequestException:
                time.sleep(5)
        time.sleep(0.2)
    cache.write_text(json.dumps(known))
    return known


def go_info(ids):
    cache = REF / "go_info.json"
    known = json.loads(cache.read_text()) if cache.exists() else {}
    todo = [i for i in ids if i not in known]
    for k in range(0, len(todo), 100):
        chunk = todo[k:k + 100]
        r = requests.get("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/" + ",".join(chunk),
                         headers={"Accept": "application/json"}, timeout=120)
        if r.ok:
            for t in r.json().get("results", []):
                known[t["id"]] = dict(name=t.get("name", ""), aspect=t.get("aspect", ""),
                                      obsolete=t.get("isObsolete", False))
        time.sleep(0.2)
    cache.write_text(json.dumps(known))
    return known


def go_enrichment(sets, background):
    g = pd.read_csv(CUGEN / "97103_GO_v2.5.gaf.txt.gz", sep="\t", comment="!", header=None, dtype=str)
    g = g[[1, 4]].dropna()
    g.columns = ["gene", "go"]
    g["gene"] = g.gene.str.split(".").str[0]
    g = g[g.gene.isin(background)].drop_duplicates()
    anc = go_ancestors(sorted(g.go.unique()))
    gene2go = defaultdict(set)
    for gene, go in zip(g.gene, g.go):
        gene2go[gene].update(anc.get(go, [go]))
        gene2go[gene].add(go)
    term2gene = defaultdict(set)
    for gene, gos in gene2go.items():
        for t in gos:
            term2gene[t].add(gene)
    universe = set(gene2go)
    Nu = len(universe)
    info = go_info(sorted(term2gene))
    rows = []
    for name, genes in sets.items():
        hits = set(genes) & universe
        n = len(hits)
        for t, tg in term2gene.items():
            K = len(tg)
            if K < 5 or K > 2000:
                continue
            k = len(hits & tg)
            if k < 3:
                continue
            i = info.get(t, {})
            if i.get("obsolete"):
                continue
            rows.append(dict(set=name, GO=t, term=i.get("name", ""), aspect=i.get("aspect", ""), k=k, n=n, K=K,
                             N=Nu, fold=(k / n) / (K / Nu), p=hypergeom.sf(k - 1, Nu, K, n),
                             genes=",".join(sorted(hits & tg))))
    df = pd.DataFrame(rows)
    df["FDR"] = df.groupby("set").p.transform(lambda p: bh(p.values))
    df.sort_values(["set", "p"]).to_csv(F / "go_enrichment.tsv", sep="\t", index=False)
    for s in sets:
        top = df[(df.set == s) & (df.FDR < 0.05) & (df.aspect == "biological_process")].nsmallest(8, "p")
        print(f"  GO {s}: {((df.set == s) & (df.FDR < 0.05)).sum()} terms FDR<0.05; top BP: "
              + "; ".join(top.term), flush=True)


# ---------------- Arabidopsis ----------------
def arabidopsis(cds, background):
    c = ["q", "s", "pid", "len", "mm", "gap", "qs", "qe", "ss", "se", "ev", "bits"]
    a = pd.read_csv(REF / "wm_at.tsv", sep="\t", names=c).sort_values("bits", ascending=False).drop_duplicates("q")
    b = pd.read_csv(REF / "at_wm.tsv", sep="\t", names=c).sort_values("bits", ascending=False).drop_duplicates("q")
    rbh = a.merge(b[["q", "s"]], left_on=["q", "s"], right_on=["s", "q"], suffixes=("", "_b"))
    best = a[a.ev < 1e-10][["q", "s", "pid"]].rename(columns={"q": "watermelon", "s": "AGI"})
    best["RBH"] = best.watermelon.isin(rbh.q)
    syms = {}
    with gzip.open(REF / "Arabidopsis_thaliana.TAIR10.pep.all.fa.gz", "rt") as f:
        for l in f:
            if l.startswith(">"):
                gid = re.search(r"gene:(\S+)", l).group(1)
                m = re.search(r"gene_symbol:(\S+)", l)
                d = re.search(r"description:(.*?) \[", l)
                syms[gid] = (m.group(1) if m else "", d.group(1) if d else "")
    best["at_symbol"] = best.AGI.map(lambda x: syms.get(x, ("", ""))[0])
    best["at_description"] = best.AGI.map(lambda x: syms.get(x, ("", ""))[1])
    dro = pd.read_csv(REF / "at_water_deprivation.tsv", sep="\t")
    func = set(dro[dro.functional].AGI)
    anyev = set(dro.AGI)
    best["at_drought_functional"] = best.AGI.isin(func)
    best["at_drought_any"] = best.AGI.isin(anyev)
    best.to_csv(F / "arabidopsis_orthologs.tsv", sep="\t", index=False)
    bmap = best.set_index("watermelon")
    bg = [g for g in background if g in bmap.index]
    rows = []
    for name, genes in [("CDS_all", list(cds.index)), ("CDS_up", list(cds.index[cds.direction == "up"])),
                        ("CDS_down", list(cds.index[cds.direction == "down"]))]:
        gs = [g for g in genes if g in bmap.index]
        for col in ["at_drought_functional", "at_drought_any"]:
            a_ = bmap.loc[gs, col].sum()
            b_ = len(gs) - a_
            c_ = bmap.loc[bg, col].sum() - a_
            d_ = len(bg) - len(gs) - c_
            OR, p = fisher_exact([[a_, b_], [c_, d_]], alternative="greater")
            rows.append(dict(set=name, evidence=col, hits=int(a_), n=len(gs), bg_rate=bmap.loc[bg, col].mean(),
                             set_rate=a_ / max(1, len(gs)), odds_ratio=OR, p=p))
    res = pd.DataFrame(rows)
    res.to_csv(F / "arabidopsis_enrichment.tsv", sep="\t", index=False)
    print(res[["set", "evidence", "hits", "n", "set_rate", "bg_rate", "odds_ratio", "p"]].round(4).to_string(index=False))
    return best


# ---------------- TF families ----------------
TF = [("NAC", r"\bNAC|NO APICAL MERISTEM|ANAC|ATAF"), ("WRKY", r"WRKY"), ("MYB", r"\bMYB"),
      ("bZIP/ABF", r"bZIP|\bABF\d|ABA[- ]responsive element[- ]binding|ABI5|AREB"),
      ("AP2/ERF/DREB", r"AP2|ERF\d|ETHYLENE[- ]RESPONSIVE|DREB|RAP2"), ("HSF", r"\bHSF|heat stress transcription"),
      ("HD-ZIP", r"HD-?ZIP|homeobox[- ]leucine zipper|ATHB"), ("bHLH", r"bHLH|basic helix"),
      ("C2H2/ZAT", r"\bZAT\d|C2H2|zinc finger protein ZAT"), ("NF-Y", r"NF-?Y|nuclear transcription factor Y"),
      ("GRAS", r"GRAS|SCARECROW"), ("TCP", r"\bTCP\d")]


def tf_families(cds, ortho):
    ahrd = pd.read_csv(CUGEN / "97103_AHRD_v2.5.txt.gz", sep="\t", header=None, usecols=[0, 1],
                       names=["gene", "desc"], dtype=str).drop_duplicates("gene").set_index("gene").desc
    o = ortho.set_index("watermelon")
    rows = []
    for g, r in cds.iterrows():
        txt = " ".join([str(ahrd.get(g, "")), str(o.at_symbol.get(g, "")) if g in o.index else "",
                        str(o.at_description.get(g, "")) if g in o.index else ""])
        fam = next((f for f, pat in TF if re.search(pat, txt, re.I)), None)
        if fam:
            rows.append(dict(gene=g, family=fam, direction=r.direction, meta_log2FC=r.meta_log2FC,
                             description=ahrd.get(g, ""), at_symbol=o.at_symbol.get(g, "") if g in o.index else ""))
    t = pd.DataFrame(rows)
    t.to_csv(F / "cds_transcription_factors.tsv", sep="\t", index=False)
    if len(t):
        print("  TFs in CDS:", t.groupby(["family", "direction"]).size().unstack(fill_value=0).to_dict(), flush=True)


# ---------------- promoters ----------------
MOTIFS = {"ABRE": "ACGTG[GT]C", "DRE/CRT": "[AG]CCGAC", "MYB (MYB2AT)": "TAACTG", "MYC2 (E-box)": "CACATG",
          "W-box": "TTGAC[CT]", "G-box": "CACGTG"}


def rc_regex(pat):
    """Reverse complement of a simple DNA regex made of bases and [..] classes."""
    comp = str.maketrans("ACGT", "TGCA")
    toks = re.findall(r"\[[ACGT]+\]|[ACGT]", pat)
    out = []
    for t in reversed(toks):
        out.append("[" + "".join(sorted(t[1:-1].translate(comp))) + "]" if t.startswith("[") else t.translate(comp))
    return "".join(out)


def revcomp(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def promoters(cds, background, length=1500):
    out = F / "promoters_1500.tsv"
    if not out.exists():
        starts = {}
        with gzip.open(REF / "Citrullus_lanatus.Cla97_v1.63.chr.gff3.gz", "rt") as f:
            for l in f:
                p = l.split("\t")
                if len(p) < 9 or p[2] != "CDS":
                    continue
                gid = re.search(r"Parent=transcript:([^;\s]+)", p[8])
                if not gid:
                    continue
                g = re.sub(r"\.\d+$", "", gid.group(1))
                s, e, strand = int(p[3]), int(p[4]), p[6]
                pos = s if strand == "+" else e
                if g not in starts:
                    starts[g] = (p[0], pos, strand)
                else:
                    c, q, st = starts[g]
                    starts[g] = (c, min(q, pos) if strand == "+" else max(q, pos), st)
        seqs, name, buf = {}, None, []
        with gzip.open(REF / "Citrullus_lanatus.Cla97_v1.dna.toplevel.fa.gz", "rt") as f:
            for l in f:
                if l.startswith(">"):
                    if name:
                        seqs[name] = "".join(buf)
                    name, buf = l[1:].split()[0], []
                else:
                    buf.append(l.strip().upper())
            seqs[name] = "".join(buf)
        rows = []
        for g, (c, pos, strand) in starts.items():
            if c not in seqs:
                continue
            if strand == "+":
                s = seqs[c][max(0, pos - 1 - length):pos - 1]
            else:
                s = revcomp(seqs[c][pos:pos + length])
            rows.append((g, s))
        pd.DataFrame(rows, columns=["gene", "seq"]).to_csv(out, sep="\t", index=False)
    pr = pd.read_csv(out, sep="\t").set_index("gene").seq.dropna()
    counts = {}
    for name, pat in MOTIFS.items():
        fwd = re.compile(f"(?=({pat}))")
        rcp = rc_regex(pat)
        rev = re.compile(f"(?=({rcp}))")
        palin = rcp == pat
        counts[name] = pr.map(lambda s, f=fwd, r=rev, pal=palin: len(f.findall(s)) + (0 if pal else len(r.findall(s))))
    cnt = pd.DataFrame(counts)
    cnt.to_csv(F / "promoter_motif_counts.tsv", sep="\t")
    bg = [g for g in background if g in cnt.index]
    rows = []
    for name, genes in [("CDS_up", cds.index[cds.direction == "up"]), ("CDS_down", cds.index[cds.direction == "down"])]:
        gs = [g for g in genes if g in cnt.index]
        for mname in MOTIFS:
            thr = 2 if mname == "ABRE" else 1        # ABA response needs >=2 ABREs (coupling)
            a = (cnt.loc[gs, mname] >= thr).sum()
            c = (cnt.loc[bg, mname] >= thr).sum() - a
            OR, p = fisher_exact([[a, len(gs) - a], [c, len(bg) - len(gs) - c]])
            rows.append(dict(set=name, motif=mname, threshold=f">={thr}", set_rate=a / len(gs),
                             bg_rate=(cnt.loc[bg, mname] >= thr).mean(), odds_ratio=OR, p=p, n=len(gs)))
    res = pd.DataFrame(rows)
    res["FDR"] = bh(res.p.values)
    res.to_csv(F / "promoter_enrichment.tsv", sep="\t", index=False)
    print(res.round(4).to_string(index=False), flush=True)


def main():
    cds = pd.read_csv(R / "cds.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(R / "meta.tsv.gz", sep="\t", index_col=0)
    background = list(meta.index)
    go_enrichment({"CDS_up": list(cds.index[cds.direction == "up"]),
                   "CDS_down": list(cds.index[cds.direction == "down"])}, background)
    ortho = arabidopsis(cds, background)
    tf_families(cds, ortho)
    promoters(cds, background)


if __name__ == "__main__":
    main()
