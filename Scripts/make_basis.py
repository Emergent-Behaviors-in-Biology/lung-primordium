import pandas as pd
import numpy as np
import pickle
import lap_v2_py3 as lap_v2

reprocess_new_basis = True

#Folder with data:
source_folder = '../Source_Data/'
dest_folder = '../Processed_Data/'


if reprocess_new_basis:
    #Load in the conversion table
    conv = pd.read_csv(source_folder+'ann.csv',usecols=[1,2],index_col=0,squeeze=True)
    #Throw out ambiguous reads
    conv = conv[~conv.index.duplicated(keep=False)]
    #Load in the new basis data
    LU = pd.read_excel(source_folder+'InVivo.xlsx',sheet_name='LU_POS vs FG',index_col=0,usecols=[0,5,6,7]).reindex(index=conv.keys()).mean(axis=1)
    FG = pd.read_excel(source_folder+'InVivo.xlsx',sheet_name='LU_POS vs FG',index_col=0,usecols=[0,2,3,4]).reindex(index=conv.keys()).mean(axis=1)
    TH = pd.read_excel(source_folder+'InVivo.xlsx',sheet_name='TH_POS vs FG',index_col=0,usecols=[0,5,6,7]).reindex(index=conv.keys()).mean(axis=1)
    BR = pd.read_excel(source_folder+'InVivo.xlsx',sheet_name='BR_POS vs ECT',index_col=0,usecols=[0,2,3,4]).reindex(index=conv.keys()).mean(axis=1)
    ECT = pd.read_excel(source_folder+'InVivo.xlsx',sheet_name='BR_POS vs ECT',index_col=0,usecols=[0,5,6,7]).reindex(index=conv.keys()).mean(axis=1)
    newdict = {'E.9 Lung Prim Nkx+':LU,'E.13 Thyroid':TH,'E.8.25 Foregut Endoderm':FG,'E.9 Forebrain':BR,'E.8.25 Ectoderm':ECT}

    #Reindex using Entrez ID's, add name, and average duplicates
    for name in newdict.keys():
        newdict[name].index=conv
        newdict[name].dropna(inplace=True)
        newdict[name].name = name
        temp = newdict[name].copy()
        newdict[name] = newdict[name][~newdict[name].index.duplicated()]
        for item in newdict[name].index:
            newdict[name].loc[item] = temp.loc[item].mean()
        del temp
    
    f = open(dest_folder+'NewBasis.dat','wb')
    pickle.dump(newdict,f)
    f.close()
    
else:
    f = open(dest_folder+'NewBasis.dat','rb')
    newdict = pickle.load(f)
    f.close()

#%% Load in the basis data
basis = pd.read_csv(source_folder+'Mouse_Basis_Data.txt',sep='\t',index_col=0,usecols=[0]+list(range(3,64)))

#####################
# Append new basis data
#thresh = 1
basis_new = basis.copy()
newdict_log = {}
#for name in newdict.keys():
#    newdict_log[name] = np.log2(newdict[name]+1)
#    basis_new = basis_new.join(newdict_log[name][newdict_log[name]>thresh],how='inner')
for name in newdict.keys():
    basis_new = basis_new.join(newdict[name],how='inner')
basis_new.dropna(inplace=True)

####################
#Load Keri's data
keri_index = pd.read_csv(source_folder+'entrez_id.txt',index_col=None,header=None).squeeze().values
keri_label = pd.read_csv(source_folder+'keri_cell_lbl.txt',index_col=None,header=None).squeeze()
keri = pd.read_csv(source_folder+'keri_ranknorm_data_corr.txt',index_col=None,header=None,sep='\t')
keri.index=keri_index
keri = keri.rename(columns=keri_label)

####################
# Load original project data 
data = pd.read_excel(source_folder+'InVitroNew.xlsx',index_col=0,usecols=[0]+list(range(38,50)))

####################
#Load new data, reindex using Entrez ID's, and average duplicates
#Load in the conversion table
conv = pd.read_csv(source_folder+'AnnotationTable_mogene10sttranscriptcluster.txt',index_col=1,usecols=[2,3],squeeze=True,sep='\t')
#Throw out ambiguous reads
conv = conv[~conv.index.duplicated(keep=False)]
data_new = pd.read_excel(source_folder+'fpkm_7-2018.xlsx',index_col=0).reindex(index=conv.keys())
data_new.index=conv
data_new.dropna(inplace=True)
#data_new = data_new[(data_new.T != 0).any()]
temp = data_new.copy()
data_new = data_new[~data_new.index.duplicated()]
#Right now, not averaging duplicates, because it is slow and doesn't matter
#for label in data_new.keys():
#    for item in data_new.index:
#        data_new.loc[item,label] = temp[label].loc[item].mean()
del temp

####################
#Load 13.5 hepatoblast
hepatoblast = pd.read_csv(source_folder+'E13.5.hepatoblast.csv',index_col=0).mean(axis=1).squeeze().reindex(conv.index)
hepatoblast.index = conv.values
hepatoblast.dropna(inplace=True)
hepatoblast = hepatoblast[~hepatoblast.index.duplicated()]
hepatoblast.name = 'E.13.5 Hepatoblast'

####################
#Load 16.5 lung
lung = pd.read_excel(source_folder+'GSE57391_count.xlsx',index_col=0).mean(axis=1).squeeze()
lung.name = 'E.16.5 Lung'

####################
#Load old test data, reindex using Entrez ID's, and average duplicates
conv_probe = pd.read_csv(source_folder+'AnnotationTable_mogene10sttranscriptcluster.txt',index_col=0,usecols=[0,2],squeeze=True,sep='\t')
data_old = pd.read_csv(source_folder+'InVitroOld.txt',index_col=0,skipfooter=1,engine='python',sep='\t').reindex(index=conv_probe.keys())
data_old.index = conv_probe
data_old_dedup = data_old[~data_old.index.duplicated()].copy()
data_old_dedup.dropna(inplace=True)
#for item in data_old_dedup.index:
#    data_old_dedup.loc[item] = data_old.loc[item].mean()

###################
#Load single-cell data and reindex using Entrez ID's, dropping duplicates
TimePoints = ['E9','E13','E15.5','E17.5']
df0=pd.read_csv(source_folder+TimePoints[0]+'/res.0.5.clusterAverages.tsv',sep='\t')
SingleCell = pd.DataFrame(index=df0.index)
conv = pd.read_csv(source_folder+'AnnotationTable_mogene10sttranscriptcluster.txt',index_col=0,usecols=[1,2],squeeze=True,sep='\t')
conv = conv[~conv.duplicated()]
conv = conv[~conv.index.duplicated()]
for TimePoint in TimePoints:
    df=pd.read_csv(source_folder+TimePoint+'/res.0.5.clusterAverages.tsv',sep='\t')
    for item in df:
        SingleCell[TimePoint+' cluster '+str(item)] = df[item]
SingleCell = SingleCell.reindex(conv.index)
SingleCell.index = conv.values
#SingleCell = SingleCell[(SingleCell.T != 0).any()].dropna()
SingleCell.dropna(inplace=True)

####################
# Combine into master table
master_table = basis_new.join(hepatoblast,how='inner')
master_table = master_table.join(lung,how='inner')
master_table = master_table.join(data,how='inner')
master_table = master_table.join(data_new,how='inner')
master_table = master_table.join(data_old_dedup,how='inner')
master_table = master_table.join(keri,how='inner')
master_table = master_table.join(SingleCell,how='inner')
master_table_old = master_table.copy()

# Rank-norm data
N = np.shape(master_table)[0]
for item in master_table:
    master_table[item] = lap_v2.rank_norm(np.asarray(master_table[item]), dist='normal',norm=N)

data_keys = data.keys()
data_keys = data_keys.append(data_old.keys())
data_keys = data_keys.append(data_new.keys())
data_keys = data_keys.append(keri.keys())

f = open(dest_folder+'master_table.dat','wb')
pickle.dump([master_table,data_keys],f)
f.close()