import requests, time, json, pandas as pd, io
E="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
queries = {
 ("watermelon","drought"): '"Citrullus"[Organism] AND (drought OR "water deficit" OR "water stress" OR PEG OR osmotic OR dehydration OR mannitol)',
 ("watermelon","salt"): '"Citrullus"[Organism] AND (salt OR salinity OR NaCl)',
 ("watermelon","temperature"): '"Citrullus"[Organism] AND (heat OR cold OR chilling OR "low temperature" OR "high temperature")',
 ("melon","drought"): '"Cucumis melo"[Organism] AND (drought OR "water deficit" OR "water stress" OR PEG OR osmotic OR dehydration)',
 ("cucumber","drought"): '"Cucumis sativus"[Organism] AND (drought OR "water deficit" OR "water stress" OR PEG OR osmotic OR dehydration)',
}
rows=[]
for (crop,stress),q in queries.items():
    for db in ["bioproject","sra"]:
        r=requests.get(E+"esearch.fcgi",params={"db":db,"term":q+(' AND "transcriptomic"[Filter]' if False else ''),"retmax":500,"retmode":"json"},timeout=60).json()
        ids=r["esearchresult"]["idlist"]; time.sleep(0.4)
        if db=="bioproject" and ids:
            s=requests.get(E+"esummary.fcgi",params={"db":"bioproject","id":",".join(ids),"retmode":"json"},timeout=60).json()["result"]
            for i in ids:
                d=s[i]; rows.append(dict(crop=crop,stress=stress,src="bioproject",acc=d.get("project_acc"),title=d.get("project_title"),type=d.get("project_data_type")))
        if db=="sra" and ids:
            # runinfo gives BioProject per run
            for k in range(0,len(ids),200):
                t=requests.post(E+"efetch.fcgi",data={"db":"sra","id":",".join(ids[k:k+200]),"rettype":"runinfo","retmode":"text"},timeout=120).text
                df=pd.read_csv(io.StringIO(t),dtype=str)
                df=df[df.LibraryStrategy=="RNA-Seq"] if "LibraryStrategy" in df else df
                for bp,g in df.groupby("BioProject"):
                    rows.append(dict(crop=crop,stress=stress,src="sra",acc=bp,title="",type=f"{len(g)} RNA-Seq runs"))
                time.sleep(0.4)
        time.sleep(0.4)
c=pd.DataFrame(rows)
c.to_csv("candidates_raw.tsv",sep="\t",index=False)
print(c.groupby(["crop","stress"]).acc.nunique())
