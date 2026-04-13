import pandas as pd
import numpy as np
import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm
import scanpy as sc
from anndata import AnnData
from scipy.io import mmread
from scipy.sparse import csr_matrix
import pyarrow.parquet as pq
import matplotlib
import matplotlib.pyplot as plt
import importlib

import warnings
warnings.simplefilter(action='ignore')
import pycisTopic
pycisTopic.__version__

import scenicplus

bs_mtx = mmread(os.path.join("/data/users/tianxiaojun2/tianxiaojun2_1fd87811099740d497851fec2551d35f/online/RG_scplus/rg.peak.counts.mtx")).astype("float32")
bs_mtx = csr_matrix(bs_mtx)
cell_ids = pd.read_csv(os.path.join("/data/users/tianxiaojun2/tianxiaojun2_1fd87811099740d497851fec2551d35f/online/RG_scplus/rg.atac.barcode.txt"), header = None).values[:,0]
region_ids = pd.read_csv(os.path.join("/data/users/tianxiaojun2/tianxiaojun2_1fd87811099740d497851fec2551d35f/online/RG_scplus/rg.atac.region.txt"), header = None).values[:, 0]
meta_info = pd.read_csv(os.path.join("/data/users/tianxiaojun2/tianxiaojun2_1fd87811099740d497851fec2551d35f/online/RG_scplus/rg.atac.meta.csv"),index_col=0)

import anndata as ad
adata  = ad.AnnData(bs_mtx)
# 设置行名（细胞名）
adata.obs_names = region_ids
adata.var_names = cell_ids

peak_counts = bs_mtx.toarray()
df = pd.DataFrame(peak_counts, index=region_ids, columns=cell_ids)
df

# Output directory
outDir = '/data/work/0414_rg/scplus/'
import os
if not os.path.exists(outDir):
    os.makedirs(outDir)
    
# Create cisTopic object
from pycisTopic.cistopic_class import *
#count_matrix=df
path_to_blacklist='/data/users/sunyunong/sunyunong_c1aeb974968047c5bd6c186742dbed88/online/diabetes/2501/0115/black_list.bed'
blacklist = pr.read_bed(path_to_blacklist)
blacklist

cistopic_obj = create_cistopic_object(fragment_matrix=df,path_to_blacklist=path_to_blacklist)

# Adding cell information
#cell_data =  pd.read_csv('/data/users/sunyunong/sunyunong_f0d1bd5b33894fcbb9c5c1fd7059b0cb/online/hypothalamus_atac/2501/0115/bs_atac_5w_meta.csv', sep=',')
meta_info.index = meta_info.index + '___cisTopic'
meta_info

import os
#path_to_mallet_binary='mallet'
# Run models
models=run_cgs_models(cistopic_obj,
                    n_topics=[20,30,40],
                    n_cpu=64,
                    n_iter=500, 
                    random_state=555,
                    alpha=50,
                    alpha_by_topic=True,
                    eta=0.1,
                    eta_by_topic=False, #Use SCRATCH if many models or big data set
                    save_path=None,
  #                  tmp_path='/data/work/brain_mic/tmp/',
                    ignore_reinit_error=True)

# Save
import pickle 
path=os.path.join(outDir,'Mallet_models_500.pkl')

with open(path, 'wb') as f:
  pickle.dump(models, f)

#os.mkdir(outDir+'models')
path=os.path.join(outDir,'model_selection.pdf')
model=evaluate_models(models,
                     select_model=None, 
                     return_model=True, 
                     metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                     plot_metrics=False,
                     save= path)
cistopic_obj.add_LDA_model(model)

# Save
import pickle 
path=os.path.join(outDir,'bm_atac_cistopic_obj2.pkl')

with open(path, 'wb') as f:
  pickle.dump(cistopic_obj, f)

import pandas as pd

# 创建一个示例DataFrame
data = cistopic_obj.cell_data
df = pd.DataFrame(data)

# 将DataFrame输出为CSV文件
df.to_csv('/data/work/0414_rg/scplus/cistopic_meta.csv', index=True)

meta_info_sorted = meta_info.reindex(data.index)
meta_info_sorted

data_reduced = data.iloc[:, :6]

# 将df1的全部列添加到df2_reduced的右侧
cisobj_meta = pd.concat([data_reduced, meta_info_sorted], axis=1)
cisobj_meta

cistopic_obj.add_cell_data(cisobj_meta)

from pycisTopic.clust_vis import *
find_clusters(cistopic_obj,
                 target  = 'cell',
                 k = 10,
                 res = [0.6], 
                 prefix = 'pycisTopic_', 
                 scale = True,
                 split_pattern = '-')

run_umap(cistopic_obj,
                 target  = 'cell', scale=True)

visual_path = outDir + '/visualization'
if not os.path.exists(visual_path):
    os.mkdir(visual_path)

plot_metadata(cistopic_obj, 
             reduction_name = 'UMAP',
             variables = ['celltype2','pycisTopic_leiden_10_0.6'], 
             target = 'cell',  num_columns = 1,
             text_size = 10, dot_size = 5,
             figsize = (5,5),
             save = outDir + '/visualization/dimensionality_reduction_label.pdf')

cell_topic_heatmap(cistopic_obj, 
                  variables = ['celltype2'], 
                  scale = True, 
                  legend_loc_x = 1.05, 
                  legend_loc_y = -1.2, 
                  legend_dist_y = -1,
                  figsize = (10, 10), 
                  save = outDir + 'visualization/heatmap_topic_contr.pdf')
                  
with open(outDir + 'macaca_cluster_cistopicObject.pkl', 'wb') as f:
    pickle.dump(cistopic_obj, f)

run_tsne(cistopic_obj,
                 target  = 'cell', scale=True)

#os.mkdir(outDir+'/visualization')
plot_metadata(cistopic_obj,
                 reduction_name='tSNE',
                 variables=[ 'celltype2', 'pycisTopic_leiden_10_0.6'], # Labels from RNA and new clusters
                 target='cell', num_columns=3,
                 text_size=10, 
                 dot_size=5,
                 figsize=(15,5),
                 save= outDir + 'visualization/dimensionality_reduction_label.pdf')

plot_topic(cistopic_obj,
            reduction_name = 'tSNE',
            target = 'cell',
            num_columns=5,
            save= outDir + 'visualization/dimensionality_reduction_topic_contr.pdf')

cell_topic_heatmap(cistopic_obj,
                     variables = ['celltype2'],
                     scale = False,
                     legend_loc_x = 1.05,
                     legend_loc_y = -1.2,
                     legend_dist_y = -1,
                     figsize=(10,10),
                     save = outDir + 'visualization/heatmap_topic_contr.pdf')

# Load cisTopic object
import pickle
infile = open(outDir + 'bm_atac_cisTopicObject4.pkl', 'rb')
cistopic_obj = pickle.load(infile)
infile.close()

os.mkdir(outDir+'topic_binarization')
from pycisTopic.topic_binarization import *
region_bin_topics = binarize_topics(cistopic_obj, method='otsu', 
                                    ntop=3000, plot=True, num_columns=5, save= outDir + 'topic_binarization/otsu.pdf')

binarized_cell_topic = binarize_topics(cistopic_obj, target='cell', method='li', plot=True, num_columns=5, nbins=60)
# Save

with open(outDir + 'topic_binarization/binarized_cell_topic.pkl', 'wb') as f:
    pickle.dump(binarized_cell_topic, f)
with open(outDir + 'topic_binarization/binarized_topic_region.pkl', 'wb') as f:
    pickle.dump(region_bin_topics, f)

from pycisTopic.diff_features import *
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
