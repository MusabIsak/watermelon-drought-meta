import requests, pandas as pd, io, time, re
E="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
orgs={"watermelon":'"Citrullus"[Organism]',"melon":'"Cucumis melo"[Organism]',"cucumber":'"Cucumis sativus"[Organism]'}
kw=re.compile(r"drought|water[ -]?(deficit|stress|withh|limit)|dehydrat|osmotic|PEG|polyethylene glycol|mannitol|desiccat|salt|salin|NaCl|ABA|abscisic",re.I)
allrows=[]
for crop,o in orgs.items():
    r=requests.get(E+"esearch.fcgi",params={"db":"bioproject","term":o+' AND ("transcriptome or gene expression"[Project Data Type] OR transcriptome OR RNA-seq)',"retmax":5000,"retmode":"json"},timeout=60).json()
    ids=r["esearchresult"]["idlist"]; print(crop,len(ids))
    for k in range(0,len(ids),200):
        s=requests.post(E+"esummary.fcgi",data={"db":"bioproject","id":",".join(ids[k:k+200]),"retmode":"json"},timeout=120).json()["result"]
        for i in ids[k:k+200]:
            d=s.get(i,{}); txt=(d.get("project_title","")+" "+d.get("project_description",""))
            allrows.append(dict(crop=crop,acc=d.get("project_acc"),title=d.get("project_title",""),desc=d.get("project_description","")[:600],hit=bool(kw.search(txt))))
        time.sleep(0.4)
a=pd.DataFrame(allrows); a.to_csv("all_projects.tsv",sep="\t",index=False)
h=a[a.hit]; print(h.groupby("crop").size())
pd.set_option("display.width",300); pd.set_option("display.max_colwidth",230)
for r in h.itertuples(): print(f"\n{r.crop} {r.acc} | {r.title}\n   {r.desc[:300]}")
