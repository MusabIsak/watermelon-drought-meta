import os
import pandas as pd, re
def L(acc): return pd.read_csv(f"runsattr_{acc}.tsv" if __import__('os').path.exists(f"runsattr_{acc}.tsv") else f"runs_{acc}.tsv",sep="\t",dtype=str).fillna("")
rows=[]
def add(d,study,crop,stress,role,f):
    for r in d.itertuples():
        x=f(r)
        if x is None: continue
        cond,tissue,geno,tp=x
        rows.append(dict(run=r.run_accession,study=study,crop=crop,stress=stress,role=role,condition=cond,tissue=tissue,genotype=geno,timepoint=tp,
                         layout=r.library_layout,read_count=r.read_count,fastq_ftp=r.fastq_ftp))
# --- watermelon drought (discovery)
old=pd.read_csv(os.environ.get("PROJECT_DIR", ".") + "/metadata/metadata.tsv",sep="\t",dtype=str).fillna("")
lab={"PRJNA454040":"WM_Leaf_M20Y34_2018","PRJNA604984":"WM_Leaf_CS_2020","PRJNA326331":"WM_Root_PEG_2016"}
for acc,st in lab.items():
    d=L("PRJNA454040") if False else None
for r in old.itertuples():
    st={"Leaf_M20vsY34_2018":"WM_Leaf_M20Y34_2018","Leaf_8d_2020":"WM_Leaf_CS_2020","Root_PEG6h_2016":"WM_Root_PEG_2016"}[r.study]
    rows.append(dict(run=r.run,study=st,crop="watermelon",stress="drought",role="discovery",condition="stress" if r.condition=="drought" else "control",
                     tissue=r.tissue,genotype=r.genotype,timepoint=r.timepoint,layout=r.library_layout,read_count=r.read_count,fastq_ftp=r.fastq_ftp))
d=L("PRJNA1055702")
add(d,"WM_Leaf_PEGtime_2024","watermelon","drought","discovery",
    lambda r:(("control" if r.sample_alias.startswith("RNA0h") else "stress"),"leaf","",re.search(r"RNA(\d+h)",r.sample_alias).group(1)))
# --- cucumber drought (external validation)
d=L("PRJDB37115")
add(d,"CU_Leaf_PEG_2025","cucumber","drought","validation",
    lambda r: None if "mRNA" not in r.experiment_title else (("control" if "CTRL" in r.experiment_title else "stress"),"leaf","",""))
d=L("PRJNA1150737")
add(d,"CU_Leaf_PEG48h_2024","cucumber","drought","validation",
    lambda r:(("stress" if "drought" in r.attributes else "control"),"leaf",("tolerant" if r.sample_alias.endswith("T") else "sensitive"),"48h"))
d=L("PRJNA720243")
add(d,"CU_Leaf_graft_2022","cucumber","drought","validation",
    lambda r:(("stress" if "WS" in r.sample_alias else "control"),"leaf",("pumpkin_rootstock" if r.sample_alias.startswith("P") else "self_grafted"),""))
# --- watermelon salt / cold (specificity)
d=L("PRJNA609260")
add(d,"WM_Leaf_salt7h_2020","watermelon","salt","specificity",
    lambda r:(("control" if "Control" in r.attributes else "stress"),"leaf","CrimsonSweet","7h"))
d=L("PRJNA844416")
add(d,"WM_salt_LR_2025","watermelon","salt","specificity",
    lambda r:(("control" if r.sample_alias[1]=="C" else "stress"),("leaf" if r.sample_alias[2]=="L" else "root"),("geno_"+r.sample_alias[0]),""))
d=L("PRJNA328189")
add(d,"WM_Leaf_cold_2017","watermelon","cold","specificity",
    lambda r: None if (r.sample_alias.endswith("-mi") or "MT" in r.sample_alias) else (("control" if r.sample_alias.startswith("CK") else "stress"),"leaf","Y134","4C"))
D=pd.DataFrame(rows)
D.to_csv("design.tsv",sep="\t",index=False)
print(D.groupby(["role","study","condition"]).size().unstack(fill_value=0)); print("total",len(D))
