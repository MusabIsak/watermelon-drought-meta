import requests, pandas as pd, io, time, xml.etree.ElementTree as ET, sys
accs=sys.argv[1:]
F="run_accession,sample_accession,sample_title,sample_alias,experiment_title,library_layout,read_count,fastq_bytes,fastq_ftp,fastq_md5"
for acc in accs:
    x=requests.get("https://www.ebi.ac.uk/ena/portal/api/filereport",params={"accession":acc,"result":"read_run","fields":F,"format":"tsv"},timeout=60)
    d=pd.read_csv(io.StringIO(x.text),sep="\t",dtype=str).fillna("")
    att=[]
    for s in d.sample_accession:
        try:
            r=requests.get(f"https://www.ebi.ac.uk/ena/browser/api/xml/{s}",timeout=60); root=ET.fromstring(r.text)
            att.append("; ".join(f"{a.findtext('TAG')}={a.findtext('VALUE')}" for a in root.iter("SAMPLE_ATTRIBUTE") if not (a.findtext('TAG') or '').startswith(('ENA-','INSDC','NCBI','BioSample','geo_loc','organism'))))
        except Exception as e: att.append("ERR")
        time.sleep(0.15)
    d["attributes"]=att; d.to_csv(f"runsattr_{acc}.tsv",sep="\t",index=False)
    print(f"\n===== {acc} ({len(d)} runs)")
    for r in d.itertuples():
        print(f"  {r.run_accession} {r.library_layout[:2]} {int(r.read_count)/1e6:5.1f}M | {r.sample_alias[:28]:28} | {r.experiment_title[-60:]:60} | {r.attributes[:150]}")
