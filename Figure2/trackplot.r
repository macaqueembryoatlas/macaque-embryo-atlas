library(Seurat)
library(ArchR)

meta<-readRDS(path)
colnames(meta)

rg1_seu<-readRDS(path)
DefaultAssay(rg1_seu)<-'RNA'

rg1_atac<-rownames(meta)[meta$celltype_f=='RG_1'&meta$tech=='atac']
length(rg1_atac)
proj<-loadArchRProject(path)
proj.sub<-subsetCells(ArchRProj = proj, cellNames = rg1_atac)

#proj降维
proj.sub <- addIterativeLSI(
    ArchRProj = proj.sub,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI_rg1", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)

#proj聚类
proj.sub <- addClusters(
    input = proj.sub,
    reducedDims = "IterativeLSI_rg1",
    method = "Seurat",
    name = "Clusters_rg1",
    resolution = 0.8
)

#加入rna数据
proj.sub <- addGeneIntegrationMatrix(
    ArchRProj = proj.sub, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI_rg1",
    seRNA = rg1_seu,
    addToArrow = FALSE,
    groupRNA = "new.cluster.ids",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

#建立peak2genelinks
proj.sub <- addPeak2GeneLinks(
    ArchRProj = proj.sub,
    reducedDims = "IterativeLSI_rg1"
)

markerGenes<-c('AKAP12','APOE','BCAN','C1QL2','C1QTNF1','DST','GDPD2','GLUL','HOPX','PAX6','SCRG1','SPARCL1','SPTAN1','TFAP2C','TNC')
p <- plotBrowserTrack(
    ArchRProj = proj.sub, 
    groupBy = "Day", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getPeak2GeneLinks(proj.sub)
)

pdf("./rg_up_trackplot.pdf", width = 6, height = 4)

for (gene in names(p)) {
  tryCatch({
    grid.newpage()
    grid.draw(p[[gene]])  # 绘制图形
    
    # 添加基因标题
    grid.text(
      label = gene,
      x = unit(0.05, "npc"), y = unit(0.95, "npc"),
      just = c("left", "top"),
      gp = gpar(fontsize = 5, col = "blue")
    )
  }, error = function(e) message("Skipped ", gene, ": ", e$message))
}

dev.off()
