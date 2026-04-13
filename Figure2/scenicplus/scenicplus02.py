import warnings
warnings.simplefilter(action = 'ignore')
import pycisTopic
pycisTopic.__version__
import numpy as np
import pandas as pd

projDir = '/data/work/0414_rg/scplus0704/'
outDir = projDir + 'output/'
import os
if not os.path.exists(outDir):
    os.makedirs(outDir)
tmpDir = '/data/work/0414_rg/scplus0704/output/tmp'
if not os.path.exists(tmpDir):
    
    os.makedirs(tmpDir)
    
saveDir = outDir + 'scenicplus/'
#
from scenicplus.scenicplus_class import SCENICPLUS, create_SCENICPLUS_object
from scenicplus.preprocessing.filtering import *

# load data
# load data
import pickle
infile = open('/data/input/Files/tianxiaojun2/embryo/data/scenic/RG_scplus_0630/output/bm_atac_cisTopicObject4.pkl', 'rb')
cistopic_obj = pickle.load(infile)
infile.close()


import pickle
infile = open('/data/input/Files/tianxiaojun2/embryo/data/scenic/RG_scplus_0630/output/DARs/Imputed_accessibility.pkl', 'rb')
imputed_acc_obj = pickle.load(infile)
infile.close()

infile = open('/data/input/Files/tianxiaojun2/embryo/data/scenic/RG_scplus_0630/output/menr.pkl', 'rb') # change here
menr = pickle.load(infile)
infile.close()

## RNA 
from loomxpy.loomxpy import SCopeLoom
from pycisTopic.loom import *
import itertools
import anndata
import scanpy as sc
import numpy as np
rna_anndata = sc.read_loom('/data/users/tianxiaojun2/tianxiaojun2_1fd87811099740d497851fec2551d35f/online/0624_RG_scplus/rg_rna2.loom')
rna_anndata

## Create Scenicplus object.
scplus_obj = create_SCENICPLUS_object(GEX_anndata = rna_anndata, 
                                     cisTopic_obj = cistopic_obj,
                                     imputed_acc_obj = imputed_acc_obj,
                                     menr = menr,
                                     multi_ome_mode = False,
                                     key_to_group_by = 'celltype2',
                                     nr_cells_per_metacells = 1)

with open('/data/work/0414_rg/scplus0704/output/scplus_obj1.pkl', 'wb') as f:
  pickle.dump(scplus_obj, f)

saveDir = outDir + 'scenicplus/'
if not os.path.exists(saveDir):
    os.mkdir(saveDir)
import scenicplus
from scenicplus.wrappers.run_scenicplus import run_scenicplus

annot=pd.read_csv('/data/users/sunyunong/sunyunong_f0d1bd5b33894fcbb9c5c1fd7059b0cb/online/hypothalamus_atac/2502/0227/macaque_gene_anno.csv')
import pyranges as pr
annot = pr.PyRanges(annot.dropna(axis=0))
chromsizes = pd.read_csv('/data/users/sunyunong/sunyunong_2f97cdfaee8345b8b1a7c5ce3b5fb531/online/code_test/scenic_plus/database/macFas5.fa.sizes', sep='\t', header=None)
chromsizes.columns = ['Chromosome', 'End']
chromsizes['Start'] = [0]*chromsizes.shape[0]
chromsizes = chromsizes.loc[:, ['Chromosome', 'Start', 'End']]

chromsizes = pr.PyRanges(chromsizes)

from scenicplus.enhancer_to_gene import get_search_space
get_search_space(
        scplus_obj,
        species = None,
        assembly = None,
        pr_annot = annot,
        pr_chromsizes = chromsizes,
        upstream = [1000, 150000],
        downstream = [1000, 150000])
with open('/data/work/0414_rg/scplus0704/output/scplus_obj2.pkl', 'wb') as f:
  pickle.dump(scplus_obj, f)

import pickle
infile = open('/data/work/0414_rg/scplus0704/output/scplus_obj2.pkl', 'rb')
scplus_obj = pickle.load(infile)
infile.close()

from scenicplus.enhancer_to_gene import get_search_space, calculate_regions_to_genes_relationships, GBM_KWARGS
calculate_regions_to_genes_relationships(scplus_obj,
                    ray_n_cpu = 20,
                    _temp_dir = '/data/work/0414_rg/scplus0704/output/tmp',
                    importance_scoring_method = 'GBM',
                    importance_scoring_kwargs = GBM_KWARGS)

from scenicplus.TF_to_gene import *
#tf_file = '/staging/leuven/stg_00002/lcb/cflerin/resources/allTFs_hg38.txt'
tf_file='/data/users/sunyunong/sunyunong_2f97cdfaee8345b8b1a7c5ce3b5fb531/online/code_test/scenic_plus/database/allTFs_macaca.txt'
calculate_TFs_to_genes_relationships(scplus_obj,
                    tf_file = tf_file,
                    ray_n_cpu = 60,
                    method = 'GBM',
                    _temp_dir = '/data/work/0414_rg/scplus0704/output/tmp',
                    key= 'TF2G_adj')

# Merge cistromes (all)
from scenicplus.cistromes import *
import time
start_time = time.time()
merge_cistromes(scplus_obj)
time = time.time()-start_time
print(time/60)

# Load functions
from scenicplus.grn_builder.gsea_approach import build_grn

build_grn(scplus_obj,
         min_target_genes = 1,
         adj_pval_thr = 1,
         min_regions_per_gene = 0,
         quantiles = (0.85, 0.90, 0.95),
         top_n_regionTogenes_per_gene = (5, 10, 15),
         top_n_regionTogenes_per_region = (),
         binarize_using_basc = True,
         rho_dichotomize_tf2g = True,
         rho_dichotomize_r2g = True,
         rho_dichotomize_eregulon = True,
         rho_threshold = 0.01,
         keep_extended_motif_annot = True,
         merge_eRegulons = True,
         order_regions_to_genes_by = 'importance',
         order_TFs_to_genes_by = 'importance',
         key_added = 'eRegulons_importance',
         cistromes_key = 'Unfiltered',
         disable_tqdm = False, #If running in notebook, set to True
         ray_n_cpu = 80,
         _temp_dir = '/data/work/0414_rg/scplus0704/output/tmp')

from scenicplus.utils import format_egrns
format_egrns(scplus_obj, eregulons_key = 'eRegulons_importance', TF2G_key = 'TF2G_adj', key_added = 'eRegulon_metadata')
scplus_obj.uns.keys()

# Format eRegulons
from scenicplus.eregulon_enrichment import *
get_eRegulons_as_signatures(scplus_obj, eRegulon_metadata_key='eRegulon_metadata', key_added='eRegulon_signatures')

## Score chromatin layer
# Region based raking
from scenicplus.cistromes import *
import time
start_time = time.time()
region_ranking = make_rankings(scplus_obj, target='region')
# Score region regulons
score_eRegulons(scplus_obj,
                ranking = region_ranking,
                eRegulon_signatures_key = 'eRegulon_signatures',
                key_added = 'eRegulon_AUC',
                enrichment_type= 'region',
                auc_threshold = 0.05,
                normalize = False,
                n_cpu = 1)
time = time.time()-start_time
print(time/60)

## Score transcriptome layer
# Gene based raking
from scenicplus.cistromes import *
import time
start_time = time.time()
gene_ranking = make_rankings(scplus_obj, target='gene')
# Score gene regulons
score_eRegulons(scplus_obj,
                gene_ranking,
                eRegulon_signatures_key = 'eRegulon_signatures',
                key_added = 'eRegulon_AUC',
                enrichment_type = 'gene',
                auc_threshold = 0.05,
                normalize= False,
                n_cpu = 1)
time = time.time()-start_time
print(time/60)

# Generate pseudobulks
import time
start_time = time.time()
generate_pseudobulks(scplus_obj,
                         variable = 'celltype2',
                         auc_key = 'eRegulon_AUC',
                         signature_key = 'Gene_based',
                         nr_cells = 5,
                         nr_pseudobulks = 100,
                         seed=555)
generate_pseudobulks(scplus_obj,
                         variable = 'celltype2',
                         auc_key = 'eRegulon_AUC',
                         signature_key = 'Region_based',
                         nr_cells = 5,
                         nr_pseudobulks = 100,
                         seed=555)
time = time.time()-start_time
print(time/60)

# Correlation between TF and eRegulons
import time
start_time = time.time()
TF_cistrome_correlation(scplus_obj,
                        variable = 'celltype2',
                        auc_key = 'eRegulon_AUC',
                        signature_key = 'Gene_based',
                        out_key = 'ACC_Seurat_cell_type_eGRN_gene_based')
TF_cistrome_correlation(scplus_obj,
                        variable = 'celltype2',
                        auc_key = 'eRegulon_AUC',
                        signature_key = 'Region_based',
                        out_key = 'ACC_Seurat_cell_type_eGRN_region_based')
time = time.time()-start_time
print(time/60)

import pandas as pd

# 复制数据
df1 = scplus_obj.uns['eRegulon_AUC']['Gene_based'].copy()
df2 = scplus_obj.uns['eRegulon_AUC']['Region_based'].copy()

# 清理列名
df1.columns = [x.split('_(')[0] for x in df1.columns]
df2.columns = [x.split('_(')[0] for x in df2.columns]

# 计算 Gene-based 和 Region-based 的相关性
correlations = df1.corrwith(df2, axis = 0)
correlations = correlations[abs(correlations) > 0.3] 

# 只保留 R2G 正相关的
keep = [x for x in correlations.index if '+_+' in x] + \
       [x for x in correlations.index if '-_+' in x]

# 只保留 extended 如果没有对应的 direct
extended = [x for x in keep if 'extended' in x]
direct = [x for x in keep if not 'extended' in x]
keep_extended = [x for x in extended if not x.replace('extended_', '') in direct]
keep = direct + keep_extended

# 只保留靶基因数>10 的
keep_gene = [x for x in scplus_obj.uns['eRegulon_AUC']['Gene_based'].columns 
             if x.split('_(')[0] in keep]
keep_gene = [x for x in keep_gene 
             if (int(x.split('_(')[1].replace('g)', '')) > 10)]

keep_all = [x.split('_(')[0] for x in keep_gene]

keep_region = [x for x in scplus_obj.uns['eRegulon_AUC']['Region_based'].columns 
               if x.split('_(')[0] in keep]

# 存储选中的 eRegulon
scplus_obj.uns['selected_eRegulons'] = {}
scplus_obj.uns['selected_eRegulons']['Gene_based'] = keep_gene
scplus_obj.uns['selected_eRegulons']['Region_based'] = keep_region

print(len(keep_gene))  

from scenicplus.plotting.dotplot import heatmap_dotplot
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt

heatmap_dotplot(
        scplus_obj = scplus_obj,
        size_matrix = scplus_obj.uns['eRegulon_AUC_filtered']['Gene_based'], #specify what to plot as dot sizes, target region enrichment in this case
        color_matrix = scplus_obj.to_df('EXP'), #specify  what to plot as colors, TF expression in this case
        scale_size_matrix = True,
        scale_color_matrix = True,
        group_variable = 'celltype2',
        subset_eRegulons = scplus_obj.uns['selected_eRegulon']['Gene_based'],
        figsize = (20, 5),
        orientation = 'horizontal',
        split_repressor_activator=True,
        save = "/data/work/0414_rg/scplus0630/output/scenicplus/heatmap_dotplot6.pdf")
