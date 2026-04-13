#0214_neuron npc coembed
library(Seurat)
library(Signac)
library(dplyr)

#加载Peak、genescore、meta
npc.peak<-readRDS('/data/users/tianxiaojun2/tianxiaojun2_41f54c6f5998467badf60418865dad72/online/1105_f.npc.peak.Rds')
npc.genescore<-readRDS('/data/users/tianxiaojun2/tianxiaojun2_41f54c6f5998467badf60418865dad72/online/1105_f.npc.genescore.Rds')

atac_meta <- readRDS("/hwfssz5/ST_SUPERCELLS/P22Z10200N0639/lizihao/projects/embryo/ATAC/out/meta_final.RDS")

neuron.peak<-readRDS("/data/users/tianxiaojun2/tianxiaojun2_41f54c6f5998467badf60418865dad72/online/1104_f.neuron.peak.Rds")
neuron.genescore<-readRDS("/data/users/tianxiaojun2/tianxiaojun2_41f54c6f5998467badf60418865dad72/online/1104_f.neuron.genescore.Rds")

peak<-cbind(neuron.peak,npc.peak)
dim(peak)
gs<-cbind(neuron.genescore,npc.genescore)
dim(gs)

cells<-c(colnames(neuron.peak),colnames(npc.peak))
length(cells)
atac_meta<-atac_meta[cells,]
dim(atac_meta)

pbmc.atac <- CreateSeuratObject(peak, assay = "ATAC", project = "ATAC")
pbmc.atac <- RunTFIDF(pbmc.atac)
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = "q0")
pbmc.atac <- RunSVD(pbmc.atac)
pbmc.atac <- RunUMAP(pbmc.atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
atac_meta <- atac_meta[colnames(pbmc.atac), ]
pbmc.atac <- AddMetaData(pbmc.atac, metadata = atac_meta)
pbmc.atac$tech <- "atac"
activity.matrix <- gs
pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
DefaultAssay(pbmc.atac) <- "ACTIVITY"
pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- ScaleData(pbmc.atac, features = rownames(pbmc.atac))

#加载rna
pbmc.rna<-readRDS("/data/users/tianxiaojun2/tianxiaojun2_fc3f1efb5edc4c3296279ee2326c8949/online/0214_neuron.npc.rna.combined.seu.Rds")
DimPlot(pbmc.rna)
# Perform standard analysis of each modality independently RNA analysis

pbmc.rna$tech <- "rna"
pbmc.rna$rna_celltype <- pbmc.rna$new.cluster.ids
integration_var_gene <- VariableFeatures(pbmc.rna)


# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, query = pbmc.atac, 
features = integration_var_gene, k.anchor = 30,
reference.assay = "integrated", query.assay = "ACTIVITY", reduction = "cca")
#saveRDS(transfer.anchors,'/data/work/0214_neuron.npc.coembed.k30.transferanchors.Rds')

#transfer
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc.rna$rna_celltype,
weight.reduction = pbmc.atac[["lsi"]],dims = 2:30)
pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)
table(celltype.predictions$prediction.score.max > 0.5)

genes.use <- integration_var_gene
refdata <- GetAssayData(pbmc.rna, assay = "integrated", slot = "data")[genes.use, ]

imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = pbmc.atac[["lsi"]],dims = 2:30)
pbmc.atac[["integrated"]] <- imputation
coembed <- merge(x = pbmc.rna, y = pbmc.atac)
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed <- FindNeighbors(coembed, reduction = "pca", dims = 1:30)
coembed <- FindClusters(coembed, resolution = 1)
coembed$celltype_f <- ifelse(is.na(coembed$rna_celltype),coembed$predicted.id,coembed$rna_celltype)
saveRDS(coembed,'/data/work/0214_neuron.npc.coembed.k30.Rds')
