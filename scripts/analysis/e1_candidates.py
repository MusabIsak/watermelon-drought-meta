"""prioritised candidate genes for functional validation (Table S12).
Output: results/extended/candidate_priority.tsv"""
import os
os.chdir(os.environ.get("PROJECT_DIR", ".") + "/results")
import pandas as pd, numpy as np
a=pd.read_csv('cds_annotated.tsv',sep='\t',index_col=0)
o=pd.read_csv('function/arabidopsis_orthologs.tsv',sep='\t',index_col=0)
gm=pd.read_csv('network/gene_modules.tsv',sep='\t').set_index('gene')
mt=pd.read_csv('network/modules.tsv',sep='\t'); IND=mt.loc[mt.r_drought.idxmax(),'module']; REP=mt.loc[mt.r_drought.idxmin(),'module']
cnt=pd.read_csv('function/promoter_motif_counts.tsv',sep='\t',index_col=0)
rbh=pd.read_csv('../ref/rbh_wm_cu.tsv',sep='\t'); oc=dict(zip(rbh.watermelon,rbh.cucumber))
cu=pd.read_csv('de/CU_Leaf_PEG_2025.tsv',sep='\t',index_col=0)
sp=pd.read_csv('validation/specificity.tsv',sep='\t',index_col=0)
tc=pd.read_csv('validation/time_course.tsv',sep='\t',index_col=0)
tf=pd.read_csv('function/cds_transcription_factors.tsv',sep='\t',index_col=0) if True else None
print(tf.columns.tolist()[:8])
d=a[['meta_log2FC','meta_fdr','I2','n_same_dir','direction','desc','At','AGI']].copy()
d['module']=gm.module.reindex(d.index); d['kME']=gm.kME.reindex(d.index)
d['ABRE_n']=cnt.ABRE.reindex(d.index); d['Gbox_n']=cnt['G-box'].reindex(d.index)
d['At_drought_any']=o.at_drought_any.reindex(d.index); d['At_drought_functional']=o.at_drought_functional.reindex(d.index)
d['cucumber_ortholog']=[oc.get(g) for g in d.index]
d['cucumber_PEG_log2FC']=[cu.log2FoldChange.get(oc.get(g),np.nan) for g in d.index]
d['cucumber_PEG_padj']=[cu.padj.get(oc.get(g),np.nan) for g in d.index]
d['salt_cold_class']=sp['class'].reindex(d.index); d['first_response_PEG']=tc.first_response.reindex(d.index)
d['TF_family']=tf.iloc[:,0].reindex(d.index) if len(tf.columns) else None
conc=(np.sign(d.cucumber_PEG_log2FC)==np.sign(d.meta_log2FC))&(d.cucumber_PEG_padj<0.05)
inmod=((d.direction=='up')&(d.module==IND))|((d.direction=='down')&(d.module==REP))
score=(d.meta_log2FC.abs().rank(pct=True)+conc.astype(float)+inmod.astype(float)*d.kME.fillna(0)+
       d.At_drought_any.fillna(False).astype(float)*0.5+d.At_drought_functional.fillna(False).astype(float)*0.5+
       ((d.ABRE_n>=2)&(d.direction=='up')).astype(float)*0.5+(d.n_same_dir==4).astype(float)*0.5)
d['priority_score']=score
d['reproducible_in_cucumber']=conc; d['in_drought_module']=inmod
d=d.sort_values('priority_score',ascending=False)
d.to_csv('extended/candidate_priority.tsv',sep='\t')
pd.set_option('display.width',250)
print(d.head(30)[['At','direction','meta_log2FC','module','kME','reproducible_in_cucumber','ABRE_n','At_drought_any','At_drought_functional','first_response_PEG','priority_score']].round(2).to_string())
