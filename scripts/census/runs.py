import requests, pandas as pd, io, time
c=pd.read_csv("candidates_raw.tsv",sep="\t",dtype=str).fillna("")
c=c.groupby("acc").agg(crop=("crop","first"),stress=("stress",lambda x: ",".join(sorted(set(x)))),title=("title",lambda x: max(x,key=len))).reset_index()
F="run_accession,study_accession,study_title,sample_title,experiment_title,sample_alias,library_layout,read_count,base_count,instrument_platform,first_public,fastq_bytes"
out=[]
for r in c.itertuples():
    for t in range(3):
        try:
            x=requests.get("https://www.ebi.ac.uk/ena/portal/api/filereport",params={"accession":r.acc,"result":"read_run","fields":F,"format":"tsv"},timeout=60)
            break
        except Exception as e: time.sleep(5)
    if x.status_code!=200 or not x.text.strip():
        out.append(dict(acc=r.acc,crop=r.crop,stress=r.stress,n=0,title=r.title)); continue
    d=pd.read_csv(io.StringIO(x.text),sep="\t",dtype=str).fillna("")
    if d.empty:
        out.append(dict(acc=r.acc,crop=r.crop,stress=r.stress,n=0,title=r.title)); continue
    d["acc"]=r.acc; d["crop"]=r.crop; d["stress_q"]=r.stress
    d.to_csv(f"runs_{r.acc}.tsv",sep="\t",index=False)
    ex=" | ".join(sorted(set((d.experiment_title.str.replace(r"^.*?sequencing: ?","",regex=True)+" / "+d.sample_title).head(40)))[:6])
    out.append(dict(acc=r.acc,crop=r.crop,stress=r.stress,n=len(d),layout=",".join(sorted(set(d.library_layout))),
        Mreads=round(pd.to_numeric(d.read_count,errors="coerce").median()/1e6,1),year=d.first_public.str[:4].min(),
        title=d.study_title.iloc[0][:110],examples=ex[:300]))
    time.sleep(0.3)
s=pd.DataFrame(out).sort_values(["crop","stress","acc"]); s.to_csv("census_summary.tsv",sep="\t",index=False)
pd.set_option("display.width",400); pd.set_option("display.max_colwidth",120)
for r in s.itertuples():
    print(f"\n[{r.crop}/{r.stress}] {r.acc} n={r.n} {getattr(r,'layout','')} med {getattr(r,'Mreads','')}M yr {getattr(r,'year','')}\n  {r.title}\n  e.g. {getattr(r,'examples','')}")
