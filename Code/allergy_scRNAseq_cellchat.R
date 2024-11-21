
#Cellchat Code 
  
setwd("C:/Users/lmrad/OneDrive - Michigan Medicine/Documents/Allergy/cellchat 21cluster SILP scRNAseq/")

int.SILP = readRDS("C:/Users/lmrad/OneDrive - Michigan Medicine/Documents/Allergy/cellchat 21cluster SILP scRNAseq/updated_21cluster_int.SILP.rds")
#BiocManager::install("BiocNeighbors")
#devtools::install_github("jokergoo/ComplexHeatmap")
#devtools::install_github("jinworks/CellChat")
library(CellChat)
library(patchwork)
library(Seurat)

int.SILP[[]]
int.SILP$samples = factor(int.SILP$sample)
Idents(int.SILP) = "Cell.types3"
DimPlot(int.SILP)

data = UpdateSeuratObject(int.SILP)
data[[]]

table(int.SILP$Cell.types3, int.SILP$sample)

Idents(data) = "Cell.types3"

DimPlot(data)

data@meta.data$cells = data@active.ident

data@meta.data$short = data@meta.data$cells
levels(data@meta.data$short)
levels(Idents(data))
Idents(data) = data@meta.data$short


data = RenameIdents(data, 
                    "Neutrophil" = "Neut",
                    
                    "Stromal Cells" = "Stromal",
                    "Plasma Cells" = "Plasma",
                    "IgG1+ Plasma Cells" = "IgG1+ Plasma"
)
data@meta.data$short = data@active.ident
DimPlot(data, reduction = 'umap', label = T)

Idents(data) = data@meta.data$sample
PBS_Ch7 = subset(data, idents = 'PBS Ch7')

PBS_Ch7@active.ident

NP_Ch7 = subset(data, idents = 'NP Ch7')


BregMarkers = read.csv("C:/Users/lmrad/OneDrive - Michigan Medicine/Documents/Allergy/Bregmarkers.csv", 
                         row.names = NULL, stringsAsFactors = FALSE)


# Make plots.
plot_list = list()
k = 1
for (i in seq.int(1, length(BregMarkers$Genes),5) ){
  j = i + 4
  p = VlnPlot(int.SILP, features = BregMarkers$Genes[i:j], 
          split.by = "sample", assay = "RNA") +
    theme(legend.position = 'right')
  plot_list[[k]] = p
  k = k + 1
}


pdf("BregMarkerVlnPLots.pdf", width = 25, height = 15)

for (i in 1:length(plot_list)) {
  print(plot_list[[i]])
}
dev.off()


VlnPlot(int.SILP, features = BregMarkers$Genes[6:10], 
        split.by = "sample", assay = "RNA") +
  theme(legend.position = 'right')
VlnPlot(int.SILP, features = BregMarkers$Genes[11:15], 
        split.by = "sample", assay = "RNA") +
  theme(legend.position = 'right')
VlnPlot(int.SILP, features = BregMarkers$Genes[16:20], 
        split.by = "sample", assay = "RNA") +
  theme(legend.position = 'right')
Idents(int.SILP) = "Cell.types3"
VlnPlot(int.SILP, features = "Ms4a1", 
        split.by = "sample", assay = "RNA") +
  theme(legend.position = 'right')


#cellchat PBS TP CH7 analysis ----

Idents(PBS_Ch7) = PBS_Ch7@meta.data$short
PBS_Ch7 = SetIdent(PBS_Ch7, value = PBS_Ch7@meta.data$short)

cellchat_PBS7 <- createCellChat(object = PBS_Ch7, group.by = 'short', assay = "RNA")

cell.labels = levels(cellchat_PBS7@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_PBS7@idents)) # number of cells in each cell group
groupSize



CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data

showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
#> Rows: 1,939
#> Columns: 11
#> $ interaction_name   <chr> "TGFB1_TGFBR1_TGFBR2", "TGFB2_TGFBR1_TGFBR2", "TGF…
#> $ pathway_name       <chr> "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "T…
#> $ ligand             <chr> "TGFB1", "TGFB2", "TGFB3", "TGFB1", "TGFB1", "TGFB…
#> $ receptor           <chr> "TGFbR1_R2", "TGFbR1_R2", "TGFbR1_R2", "ACVR1B_TGF…
#> $ agonist            <chr> "TGFb agonist", "TGFb agonist", "TGFb agonist", "T…
#> $ antagonist         <chr> "TGFb antagonist", "TGFb antagonist", "TGFb antago…
#> $ co_A_receptor      <chr> "", "", "", "", "", "", "", "", "", "", "", "", ""…
#> $ co_I_receptor      <chr> "TGFb inhibition receptor", "TGFb inhibition recep…
#> $ evidence           <chr> "KEGG: hsa04350", "KEGG: hsa04350", "KEGG: hsa0435…
#> $ annotation         <chr> "Secreted Signaling", "Secreted Signaling", "Secre…
#> $ interaction_name_2 <chr> "TGFB1 - (TGFBR1+TGFBR2)", "TGFB2 - (TGFBR1+TGFBR2…

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis


# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB)


# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 

# set the used database in the object
cellchat_PBS7@DB <- CellChatDB.use


cellchat_PBS7 <- subsetData(cellchat_PBS7) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
cellchat_PBS7 <- identifyOverExpressedGenes(cellchat_PBS7)
cellchat_PBS7 <- identifyOverExpressedInteractions(cellchat_PBS7)

groupSize
# Compute Communicaiton probability and infer network ---------------------
cellchat_PBS7 <- computeCommunProb(cellchat_PBS7, population.size = T) #remove pop.size = F
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2022-11-26 08:34:38]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2022-11-26 08:35:36]"
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_PBS7 <- filterCommunication(cellchat_PBS7, min.cells = 10)

# Extract inferred network as data frame ----------------------------------
# We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
# 
df.net.PBSch7 <- subsetCommunication(cellchat_PBS7) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
# 
df.net.path.PBSch7 <- subsetCommunication(cellchat_PBS7, slot.name = "netP") # look at level of signaling pathways

#df.net.ccl.cxcl <- subsetCommunication(cellchat, signaling = c("CCL", "CXCL"))

# Infer cell-cell communication at signaling pathway level ----------------
cellchat_PBS7 <- computeCommunProbPathway(cellchat_PBS7)


# Calculate aggregated communication network ------------------------------

cellchat_PBS7 <- aggregateNet(cellchat_PBS7)

groupSize <- as.numeric(table(cellchat_PBS7@idents))
groupSize
#setwd("C:/Users/lmrad/OneDrive - Michigan Medicine/Documents/Allergy/")

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_PBS7@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Number of interactions"
                )

colSide =getPalette(14)

names(colSide) = cell.labels

cairo_pdf(file = "CellChat_Aggregate_Weights_PBS7.pdf",width = 3, height = 3)

png(file = "CellChat_Aggregate_Weights_PBS7.png",
    units = "in",width = 3, height = 3, res = 400)
netVisual_circle(cellchat_PBS7@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength",
                  vertex.label.cex = 0.5
                 ) 

dev.off()


class(getPalette(14))

mat <- cellchat_PBS7@net$weight
getwd()
#pdf(file = paste0("C:/Users/lmrad/OneDrive - Michigan Medicine/Documents/Allergy/signaling_clusters_PBSch7", '.pdf'), onefile = TRUE)
par(mfrow = c(5,5), mar = c(1, 1, 1, 1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

dev.off()


cellchat_PBS7@netP$pathways

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat_PBS7@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat_PBS7@idents)

netVisual_aggregate(cellchat_PBS7, signaling = c("CCL"), layout = "chord")
netVisual_aggregate(cellchat_PBS7, signaling = c("CXCL"), layout = "chord")
netVisual_aggregate(cellchat_PBS7, signaling = c("IL4"), layout = "chord")
netVisual_bubble(cellchat_PBS7,  signaling = c("CCL","CXCL"), remove.isolate = T)
netVisual_bubble(cellchat_PBS7,  signaling = c("IL4"), remove.isolate = T)
netVisual_chord_gene(cellchat_PBS7, targets.use = c(1,2) ,legend.pos.x = 8, thresh = 0.01)


cairo_pdf(file = "CellChat_interactionstoTcells_PBS7.pdf",width = 7, height = 7)

netVisual_chord_gene(cellchat_PBS7,  targets.use = c(3,4,10,19), 
                     #signaling = c("IL16","IL4", "APP"),
                     #slot.name = "netP", 
                     legend.pos.x = 5,
                     show.legend = T,
                     lab.cex = 0.8,
                     scale =T)

dev.off()

pathways.show = c("IL4")
# Compute the network centrality scores
cellchat_PBS7 <- netAnalysis_computeCentrality(cellchat_PBS7, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat_PBS7, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat_PBS7)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat_PBS7, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat_PBS7, pattern = "outgoing",height = 16)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_PBS7, pattern = "incoming")
ht1 + ht2

ht <- netAnalysis_signalingRole_heatmap(cellchat_PBS7, signaling = c("CXCL", "CCL"))

library(NMF)
library(ggalluvial)

ko = selectK(cellchat_PBS7, pattern = "outgoing")
dev.off()
ko #3 or 5

nPatterns = 4

cairo_pdf(file = "CellChat_outgoingPatterns_PBS7.pdf",width = 16, height = 14)
cellchat_PBS7 <- identifyCommunicationPatterns(cellchat_PBS7, 
                                               pattern = "outgoing",
                                               k = nPatterns,
                                               width = 10,
                                               height=20,
                                               font.size = 12)
dev.off()

netAnalysis_river(cellchat_PBS7, pattern = "outgoing")
netAnalysis_dot(cellchat_PBS7, pattern = "outgoing")

dev.off()

ki = selectK(cellchat_PBS7, pattern = "incoming")
ki #4 or 6

nPatterns = 5
cellchat_PBS7 <- identifyCommunicationPatterns(cellchat_PBS7, pattern = "incoming", k = nPatterns)
netAnalysis_river(cellchat_PBS7, pattern = "incoming")
netAnalysis_dot(cellchat_PBS7, pattern = "incoming")


cairo_pdf(file = "CellChat_incomingPatterns_PBS7.pdf",width = 16, height = 14)
cellchat_PBS7 <- identifyCommunicationPatterns(cellchat_PBS7, 
                                               pattern = "incoming",
                                               k = nPatterns,
                                               width = 10,
                                               height=18,
                                               font.size = 12)
dev.off()

saveRDS(cellchat_PBS7, file = "cellchat_PBS7.rds")




#cellchat NP TP CH7 analysis ----

Idents(NP_Ch7) = NP_Ch7@meta.data$short
NP_Ch7 = SetIdent(NP_Ch7, value = NP_Ch7@meta.data$short)

cellchat_NP7 <- createCellChat(object = NP_Ch7, group.by = 'short', assay = "RNA")

cell.labels = levels(cellchat_NP7@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_NP7@idents)) # number of cells in each cell group
groupSize



CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data

showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
#> Rows: 1,939
#> Columns: 11
#> $ interaction_name   <chr> "TGFB1_TGFBR1_TGFBR2", "TGFB2_TGFBR1_TGFBR2", "TGF…
#> $ pathway_name       <chr> "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "T…
#> $ ligand             <chr> "TGFB1", "TGFB2", "TGFB3", "TGFB1", "TGFB1", "TGFB…
#> $ receptor           <chr> "TGFbR1_R2", "TGFbR1_R2", "TGFbR1_R2", "ACVR1B_TGF…
#> $ agonist            <chr> "TGFb agonist", "TGFb agonist", "TGFb agonist", "T…
#> $ antagonist         <chr> "TGFb antagonist", "TGFb antagonist", "TGFb antago…
#> $ co_A_receptor      <chr> "", "", "", "", "", "", "", "", "", "", "", "", ""…
#> $ co_I_receptor      <chr> "TGFb inhibition receptor", "TGFb inhibition recep…
#> $ evidence           <chr> "KEGG: hsa04350", "KEGG: hsa04350", "KEGG: hsa0435…
#> $ annotation         <chr> "Secreted Signaling", "Secreted Signaling", "Secre…
#> $ interaction_name_2 <chr> "TGFB1 - (TGFBR1+TGFBR2)", "TGFB2 - (TGFBR1+TGFBR2…

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis


# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB)


# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 

# set the used database in the object
cellchat_NP7@DB <- CellChatDB.use


cellchat_NP7 <- subsetData(cellchat_NP7) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
cellchat_NP7 <- identifyOverExpressedGenes(cellchat_NP7)
cellchat_NP7 <- identifyOverExpressedInteractions(cellchat_NP7)

groupSize
# Compute Communicaiton probability and infer network ---------------------
cellchat_NP7 <- computeCommunProb(cellchat_NP7, population.size = T) #remove pop.size = F
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2022-11-26 08:34:38]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2022-11-26 08:35:36]"
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_NP7 <- filterCommunication(cellchat_NP7, min.cells = 10)

# Extract inferred network as data frame ----------------------------------
# We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
# 
df.net.NPch7 <- subsetCommunication(cellchat_NP7) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
# 
df.net.path.NPch7 <- subsetCommunication(cellchat_NP7, slot.name = "netP") # look at level of signaling pathways

#df.net.ccl.cxcl <- subsetCommunication(cellchat, signaling = c("CCL", "CXCL"))

# Infer cell-cell communication at signaling pathway level ----------------
cellchat_NP7 <- computeCommunProbPathway(cellchat_NP7)


# Calculate aggregated communication network ------------------------------

cellchat_NP7 <- aggregateNet(cellchat_NP7)

groupSize <- as.numeric(table(cellchat_NP7@idents))
groupSize

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_NP7@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Number of interactions"
)

colSide =getPalette(14)

names(colSide) = cell.labels

cairo_pdf(file = "CellChat_Aggregate_Weights_NP7.pdf",width = 3, height = 3)

png(file = "CellChat_Aggregate_Weights_NP7.png",
    units = "in",width = 3, height = 3, res = 400)
netVisual_circle(cellchat_NP7@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength",
                 vertex.label.cex = 0.5
) 

dev.off()


class(getPalette(14))

mat <- cellchat_NP7@net$weight
getwd()
#pdf(file = paste0("C:/Users/lmrad/OneDrive - Michigan Medicine/Documents/Allergy/signaling_clusters_PBSch7", '.pdf'), onefile = TRUE)
par(mfrow = c(5,5), mar = c(1, 1, 1, 1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

dev.off()


cellchat_NP7@netP$pathways


# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat_NP7@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat_NP7@idents)

netVisual_aggregate(cellchat_NP7, signaling = c("CCL"), layout = "chord")
netVisual_aggregate(cellchat_NP7, signaling = c("CXCL"), layout = "chord")
netVisual_aggregate(cellchat_NP7, signaling = c("IL4"), layout = "chord")
netVisual_aggregate(cellchat_NP7, signaling = c("IL10"), layout = "chord")
netVisual_bubble(cellchat_NP7,  signaling = c("CCL","CXCL"), remove.isolate = T)
netVisual_bubble(cellchat_NP7,  signaling = c("IL10"), remove.isolate = T)
netVisual_chord_gene(cellchat_NP7, targets.use = c(1,2) ,legend.pos.x = 8, thresh = 0.01)

par(mfrow = c(1,2))
netVisual_aggregate(cellchat_PBS7, signaling = c("IL4"), layout = "chord")
netVisual_aggregate(cellchat_NP7, signaling = c("IL4"), layout = "chord")


cairo_pdf(file = "CellChat_interactionstoTcells_NP7.pdf",width = 7, height = 7)


netVisual_chord_gene(cellchat_NP7,  targets.use =  c(3,4,10,19), 
                     #signaling = c("IL16","IL4", "APP"),
                     #slot.name = "netP", 
                     legend.pos.x = 4,
                     show.legend = T,
                     lab.cex = 0.8,
                     scale =T)

dev.off()

pathways.show = c("IL4")
# Compute the network centrality scores
cellchat_NP7 <- netAnalysis_computeCentrality(cellchat_NP7, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat_NP7, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat_NP7)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat_NP7, signaling = c("IL10"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat_NP7, pattern = "outgoing", height = 16)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_NP7, pattern = "incoming")
ht1 + ht2

ht <- netAnalysis_signalingRole_heatmap(cellchat_NP7, signaling = c("CXCL", "CCL"))
ht

library(NMF)
library(ggalluvial)

ko_np = selectK(cellchat_NP7, pattern = "outgoing")
ko_np #4,6

nPatterns = 5

cairo_pdf(file = "CellChat_outgoingPatterns_NP7.pdf",width = 16, height = 14)
cellchat_NP7 <- identifyCommunicationPatterns(cellchat_NP7, 
                                               pattern = "outgoing",
                                               k = nPatterns,
                                               width = 10,
                                               height=20,
                                               font.size = 12)
dev.off()

netAnalysis_river(cellchat_NP7, pattern = "outgoing")
netAnalysis_dot(cellchat_NP7, pattern = "outgoing")

dev.off()

ki_np = selectK(cellchat_NP7, pattern = "incoming")
ki_np #3,6

nPatterns = 6
cellchat_NP7 <- identifyCommunicationPatterns(cellchat_NP7, pattern = "incoming", k = nPatterns)
netAnalysis_river(cellchat_NP7, pattern = "incoming")
netAnalysis_dot(cellchat_NP7, pattern = "incoming")


cairo_pdf(file = "CellChat_incomingPatterns_NP7.pdf",width = 16, height = 14)
cellchat_NP7 <- identifyCommunicationPatterns(cellchat_NP7, 
                                               pattern = "incoming",
                                               k = nPatterns,
                                               width = 10,
                                               height=20,
                                               font.size = 12)
dev.off()
getwd()
saveRDS(cellchat_NP7, file = "cellchat_NP7.rds")

#comparison analysis NP vs PBS ch7 ----
unique(cellchat_PBS7@idents)
unique(cellchat_NP7@idents)
object.list <- list(PBS = cellchat_PBS7, NP = cellchat_NP7)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


# two datasets can be visualized using circle plot, where red (or blue) colored edges represent increased
# (or decreased) signaling in the second dataset compared to the first one.
#red more in NP
# blue less in NP

object.list[2]

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()


netVisual_circle(cellchat_PBS7@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength",
                 vertex.label.cex = 0.5
) 


cairo_pdf(file = "CellChat_Aggregate_DifferentialInterac_NP7.pdf",width = 3, height = 3)
png(file = "CellChat_Aggregate_DifferentialInterac_NP7.png",
    units = "in",width = 3, height = 3, res = 400)

netVisual_diffInteraction(cellchat, weight.scale = T,
                          vertex.label.cex = 0.5, margin = 0.3,
                          title.name = "")
dev.off()


#(B) Heatmap showing differential number of interactions or 
# interaction strength among different cell populations across two datasets

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
#> 
#> 
cairo_pdf(file = "CellChat_Aggregate_DifferentialInterac_NP7_heatmap.pdf",width = 7, height = 4)
png(file = "CellChat_Aggregate_DifferentialInterac_NP7_heatmap.png",
    units = "in",width = 7, height = 4, res = 400)
gg1 + gg2
dev.off()

# Compare the major sources and targets in a 2D space

#(A) Identify cell populations with significant changes in sending or receiving signals
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], 
                                               weight.MinMax = weight.MinMax) + ylim(0, 0.022) +xlim(0, 0.02)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> 
cairo_pdf(file = "CellChat_Aggregate_signaling_scatter.pdf",width = 12, height = 7)
png(file = "CellChat_Aggregate_signaling_scatter.png",
    units = "in",width = 12, height = 7, res = 400)
patchwork::wrap_plots(plots = gg)
dev.off()

cell.labels
groupSize

#Identify altered signaling with distinct interaction strength
# (A) Compare the overall information flow of each signaling pathway or ligand-receptor pair

gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "count", 
               #sources.use = c(10,19), 
               #targets.use = c(10,19), 
               slot.name = "netP",
               #pairLR = c(pairLR.use.up.Bcellsource$interaction_name, pairLR.use.down.Bcellsource$interaction_name),
               stacked = F, do.stat = TRUE,)

gg1 + gg2

cairo_pdf(file = "CellChat_ranknet.pdf",width = 5, height = 8)
png(file = "CellChat_ranknet.png",
    units = "in",width = 5, height = 8, res = 400)
rankNet(cellchat, mode = "comparison", measure = "weight", 
        sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE,
        font.size = 8)
dev.off()

# (B) Compare outgoing (or incoming) signaling patterns associated with each cell population

library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", 
                                        signaling = pathway.union, title = names(object.list)[i], width = 5, height = 20)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", 
                                        signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 20)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", 
                                        signaling = pathway.union, title = names(object.list)[i], width = 5, height = 15, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", 
                                        signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 15, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all",
                                        signaling = pathway.union, title = names(object.list)[i], width = 5, height = 15, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all",
                                        signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 15, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


#Part III: Identify the up-gulated and down-regulated signaling ligand-receptor pairs
#Identify dysfunctional signaling by comparing the communication probabities
netVisual_bubble(cellchat, sources.use = cell.labels[2], targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)

#> Comparing communications on a merged object
cell.labels[4]

cell.labels[c(8, 9,10, 13,14,19)]

gg1 <- netVisual_bubble(cellchat, sources.use = cell.labels[c(8, 9,10, 13,14,19)],
                        targets.use = cell.labels[c(8, 9,10, 13,14,19)],  
                        #comparison = c(1, 2), 
                        max.dataset = 2,
                        title.name = paste("Increased signaling in", object.list[2]$NP@meta$Tx), angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = cell.labels[c(8, 9,10, 13,14,19)], 
                        targets.use = cell.labels[c(8, 9,10, 13,14,19)],
                        #comparison = c(1, 2), 
                        max.dataset = 1, 
                        title.name = paste("Decreased signaling in", object.list[2]$NP@meta$Tx ), angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2

gg1 <- netVisual_bubble(cellchat, #sources.use = 4, 
                        targets.use = 13,  comparison = c(1, 2), max.dataset = 2,
                        title.name = paste("Increased signaling in", object.list[2]$NP@meta$Tx), angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, #sources.use = 4, 
                        targets.use = 13,  comparison = c(1, 2), max.dataset = 1, 
                        title.name = paste("Decreased signaling in", object.list[2]$NP@meta$Tx ), angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2


gg1 <- netVisual_bubble(cellchat, sources.use = 13, 
                        #targets.use = 13,  
                        comparison = c(1, 2), max.dataset = 2,
                        title.name = paste("Increased signaling in", object.list[2]$NP@meta$Tx), angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 13, 
                        #targets.use = 13,  
                        comparison = c(1, 2), max.dataset = 1, 
                        title.name = paste("Decreased signaling in", object.list[2]$NP@meta$Tx ), angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2


#gg1$data

#Identify dysfunctional signaling by using differential expression analysis
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "NP"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".Ch7.merged")

# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, 
# which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1).
# Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 


cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, 
                                       thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
#> Use the joint cell labels from the merged CellChat object

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = F)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.upL <- subsetCommunication(cellchat, net = net, datasets = "NP",ligand.logFC = 0.05)
net.upR <- subsetCommunication(cellchat, net = net, datasets = "NP",receptor.logFC = 0.05)
net.up = rbind(net.upL, net.upR)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.downL <- subsetCommunication(cellchat, net = net, datasets = "PBS",ligand.logFC = -0.05)
net.downR <- subsetCommunication(cellchat, net = net, datasets = "PBS",receptor.logFC = -0.05)
net.down <- rbind(net.downL, net.downR)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

df <- findEnrichedSignaling(object.list[[2]], features = c("Ctla4"), 
                            idents = c(cell.labels), pattern ="outgoing")

# Visualize the identified up-regulated and down-regulated signaling ligand-receptor pairs

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up,
                        sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  
                        angle.x = 90, remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4,
                        targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, 
                        remove.isolate = T,
                        title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

cell.labels
par(mfrow = c(1,2), xpd=TRUE)

#png(file = "CellChat_chord_Comparison_toTcells_up.png",units = "in",width = 6, height = 6, res = 400)

netVisual_chord_gene(object.list[[2]],  
                     #sources.use = C(13),
                     #targets.use = c(2,3), 
                     signaling = "IL10",
                     slot.name = 'net', net = net.up,
                     #lab.cex = 0.7, small.gap = 3, 
                     legend.pos.x = 1,
                     #scale = T,
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

dev.off()

png(file = "CellChat_chord_Comparison_toTcells_down.png",
    units = "in",width = 6, height = 6, res = 400)
par(mfrow = c(1,1), xpd=TRUE, mar =  c(1,1,1,5))
netVisual_chord_gene(object.list[[1]],  
                     #sources.use = 14,
                     #targets.use = c(2,3), 
                     signaling = "IL10",
                     slot.name = 'net', net = net.down, lab.cex = 0.7, small.gap = 3, 
                     legend.pos.x = 10,
                     scale = T,
                     title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

dev.off()
netVisual_chord_gene(object.list[[1]],  targets.use = c(2,3),
                     slot.name = 'net', net = net.down, lab.cex = 0.3, small.gap = 3.5, 
                     title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()

#> You may try the function `netVisual_chord_cell` for visualizing individual signaling pathway







# Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram

pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

dev.off()

cairo_pdf(file = "CirclePlot_IL4_Ch7.pdf",width = 10, height = 6)
pathways.show <- c("IL4") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 5, 
                      #vertex.label.cex = 0.5,
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()


net$interaction_ST = paste0(net$interaction_name, "_",
                               net$source, "_",
                               net$target)

png(file = "Il4_geneChord_PBS_Ch7.png",units = "in",width = 6, height = 6, res = 400)
netVisual_chord_gene(object.list[[1]],  
                     #sources.use = C(13),
                     #targets.use = c(2,3), 
                     signaling = "IL4",
                     slot.name = 'net', net = net.down,
                     #lab.cex = 0.7, small.gap = 3, 
                     legend.pos.x = 10,
                     #scale = T,
                     
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[1])
                     )
dev.off()

png(file = "Il4_geneChord_NP_Ch7.png",units = "in",width = 7, height = 4, res = 400)
netVisual_chord_gene(object.list[[2]],  
                     #sources.use = C(13),
                     #targets.use = c(2,3), 
                     signaling = "IL4",
                     slot.name = 'net', net = net.up,
                     #lab.cex = 0.7, small.gap = 3, 
                     legend.pos.x = 1,
                     #scale = T,
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))



dev.off()
pathways.show <- c("MHC-II") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      #vertex.label.cex = 0.5,
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}


pairLR.MHCII <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = F)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
netVisual_individual(object.list[[1]], signaling = pathways.show, pairLR.use = pairLR.MHCII[1,], layout = "circle")
netVisual_individual(object.list[[2]], signaling = pathways.show, pairLR.use = pairLR.MHCII[1,], layout = "circle")



cairo_pdf(file = "ChordGene_MHCII_Ch7.pdf",width = 15, height = 10)
png(file = "ChordGene_MHCII_Ch7.png",res = 400,units = "in",width = 15, height = 10)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  if(i == 1){
  netVisual_chord_gene(object.list[[i]],  
                       #sources.use = C(13),
                       #targets.use = c(2,3), 
                       signaling = pathways.show,
                       slot.name = 'net', net = net.down,
                       #lab.cex = 0.7, small.gap = 3, 
                       legend.pos.x = 200,
                       scale = T,
                       
                       title.name = paste0("Up-regulated signaling in ", names(object.list)[i])
  )
  } else{
    netVisual_chord_gene(object.list[[i]],  
                         #sources.use = C(13),
                         #targets.use = c(2,3), 
                         signaling = pathways.show,
                         slot.name = 'net', net = net.up,
                         #lab.cex = 0.7, small.gap = 3, 
                         legend.pos.x = 10,
                         scale = T,
                         title.name = paste0("Up-regulated signaling in ", names(object.list)[i]))
  }
}
dev.off()

netVisual_chord_gene(object.list[[1]],  
                     #sources.use = C(13),
                     #targets.use = c(2,3), 
                     signaling = pathways.show,
                     slot.name = 'net', net = net.down,
                     #lab.cex = 0.7, small.gap = 3, 
                     legend.pos.x = 80,
                     #scale = T,
                     
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[1])
)


netVisual_chord_gene(object.list[[2]],  
                     #sources.use = C(13),
                     #targets.use = c(2,3), 
                     signaling = pathways.show,
                     slot.name = 'net', net = net.up,
                     #lab.cex = 0.7, small.gap = 3, 
                     legend.pos.x = 1,
                     scale = T,
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

dev.off()

pathways.show <- c("CD23") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[2]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 5, 
                      #vertex.label.cex = 0.5,
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

cairo_pdf(file = "ChordGene_CD23_Ch7.pdf",width = 15, height = 10)
png(file = "ChordGene_CD23_Ch7.png",res = 400,units = "in",width = 5, height = 5)

    netVisual_chord_gene(object.list[[2]],  
                         #sources.use = C(13),
                         #targets.use = c(2,3), 
                         signaling = pathways.show,
                         slot.name = 'net', net = net.up,
                         #lab.cex = 0.7, small.gap = 3, 
                         legend.pos.x = 10,
                         scale = F,
                         title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))


dev.off()



netVisual_chord_gene(object.list[[1]],  
                     #sources.use = C(13),
                     #targets.use = c(2,3), 
                     signaling = pathways.show,
                     slot.name = 'net', #net = net.down,
                     #lab.cex = 0.7, small.gap = 3, 
                     legend.pos.x = 10,
                     #scale = T,
                     
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[1])
)


netVisual_chord_gene(object.list[[2]],  
                     #sources.use = C(13),
                     #targets.use = c(2,3), 
                     signaling = pathways.show,
                     slot.name = 'net',# net = net.up,
                     #lab.cex = 0.7, small.gap = 3, 
                     legend.pos.x = 1,
                     scale = T,
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
dev.off()


cairo_pdf(file = "CirclePlot_IL10_Ch7.pdf",width = 10, height = 6)
pathways.show <- c( "TGFb") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 15, 
                      #vertex.label.cex = 0.5,
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()


pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))


# Chord diagram
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()
cell.labels

par(mfrow = c(1, 2), xpd=TRUE, cex = 0.5, cex.main = 0.3, cex.axis = 0.3, cex.sub = 0.3, cex.lab =0.3,
    mar = c(5,4,4,10))
# compare all the interactions sending from all cells to T cells
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], 
                       #sources.use = c(14),targets.use = c(2,3), 
                       lab.cex = 0.5, 
                       title.name = paste0("Signaling to T cells in ", names(object.list)[i]),
                       legend.pos.x = 0,
                       show.legend = T,
                      scale = T,
                       annotationTrackHeight = c(0.01)
                       )
}
dev.off()
cell.labels
# show all the significant signaling pathways from fibroblast to immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]],targets.use = c(2,3),slot.name = "netP", 
                       title.name = paste0("Signaling pathways sending to T cells - ", names(object.list)[i]),
                       legend.pos.x = 40,
                       scale = T)
}


# show all the significant signaling pathways from fibroblast to immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]],targets.use = 13,slot.name = "netP", 
                       title.name = paste0("Signaling pathways sending to T cells - ", names(object.list)[i]),
                       legend.pos.x = 40,
                       scale = T)
}
dev.off()
cell.labels

cairo_pdf(file = "chordpath_Ch7_B_Treg_cd103DC.pdf",width = 10, height = 6)
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  if ( i == 1){
  netVisual_chord_gene(object.list[[i]],
                       sources.use = c(10,13,14, 19),
                       targets.use = c(10,13,14, 19),slot.name = "netP", net = net.down, 
                       title.name = paste0("Treg- CD103 DC - Bcell interactions ", names(object.list)[i]),
                       legend.pos.x = 100,
                       scale = T)
  }
  else{
    netVisual_chord_gene(object.list[[i]],
                         sources.use = c(10,13,14, 19),
                         targets.use = c(10,13,14, 19),slot.name = "netP", net = net.up, 
                         title.name = paste0("Treg- CD103 DC - Bcell interactions ", names(object.list)[i]),
                         #legend.pos.x = 10,
                         scale = T)
  }
}

dev.off()

cairo_pdf(file = "chordpath_Ch7_B_allT_cd103DC.pdf",width = 10, height = 6)
png(file = "chordpath_Ch7_B_allT_cd103DC.png",width = 10, height = 6, unit = "in", res = 400)
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  if ( i == 1){
    netVisual_chord_gene(object.list[[i]],
                         sources.use = c(10:15,17:19),
                         targets.use = c(10:15,17:19),slot.name = "netP", net = net.down, 
                         title.name = paste0("T Cell - CD103 DC - B cell interactions ", names(object.list)[i]),
                         legend.pos.x = 100,
                         scale = T)
  }
  else{
    netVisual_chord_gene(object.list[[i]],
                         sources.use = c(10:15,17:19),
                         targets.use = c(10:15,17:19),slot.name = "netP", net = net.up, 
                         title.name = paste0("T Cell - CD103 DC - B cell interactions ", names(object.list)[i]),
                         legend.pos.x = -50,
                         scale = T)
  }
}
dev.off()

cell.labels
i = 1

pairLR.use.up.Bcellsource = net.up[which(net.up$source == "B Cells" | net.up$target == "B Cells" | 
                                           net.up$source == "IgG1+ B Cells" | net.up$target == "IgG1+ B Cells"), ]
pairLR.use.down.Bcellsource = net.down[which(net.down$source == "B Cells" | net.down$target == "B Cells"| 
                                               net.down$source == "IgG1+ B Cells" | net.down$target == "IgG1+ B Cells"), ]


pairLR.use.down.Tcell  = net.down[which(net.down$source %in% cell.labels[c(11:13,15,17:18)] | 
                 net.down$target %in% cell.labels[c(11:13,15,17:18)]) ,]

pairLR.use.up.Tcell  = net.up[which(net.up$source %in% cell.labels[c(11:13,15,17:18)] | 
                                          net.up$target %in% cell.labels[c(11:13,15,17:18)]) ,]

png(file = "chordpath_Ch7_allT.png",width = 14, height = 10, unit = "in", res = 400)
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  if ( i == 1){

netVisual_chord_gene(object.list[[i]],
                     #sources.use = c(10),
                     #targets.use = 10,
                     slot.name = "netP", net = pairLR.use.down.Tcell, 
                     #pairLR.use = pairLR.use.up.Bcellsource[, "interaction_name", drop = F],
                     title.name = paste0("Upregulated T cell interactions in ", names(object.list)[i]),
                     legend.pos.x = 150,
                     scale = T)

  }
  else{
netVisual_chord_gene(object.list[[i]],
                     #sources.use = c(10),
                     #targets.use = 10,
                     slot.name = "netP", net = pairLR.use.up.Tcell, 
                     #pairLR.use = pairLR.use.up.Bcellsource[, "interaction_name", drop = F],
                     title.name = paste0("Upregulated T cell interactions in ", names(object.list)[i]),
                     legend.pos.x = -50,
                     scale = T)
  }}
  dev.off()

pairLR.use.Bcell.NP= net[which((net$source == "B Cells" | net$target == "B Cells" | 
                                             net$source == "IgG1+ B Cells" | net$target == "IgG1+ B Cells")
                            &net$datasets == "NP"), ]

  
pairLR.use.Bcell.PBS = net[which((net$source == "B Cells" | net$target == "B Cells" | 
                                  net$source == "IgG1+ B Cells" | net$target == "IgG1+ B Cells")
                               &net$datasets == "PBS"), ]

net.Bcells = net[which((net$source == "B Cells" | net$target == "B Cells" | 
                          net$source == "IgG1+ B Cells" | net$target == "IgG1+ B Cells")
                       ), ]

net.Bcells$interaction_ST = paste0(net.Bcells$interaction_name, "_",
                                   net.Bcells$source, "_",
                                   net.Bcells$target)

pairLR.use.Bcell.NP$interaction_ST = paste0(pairLR.use.Bcell.NP$interaction_name, "_",
                                            pairLR.use.Bcell.NP$source, "_",
                                            pairLR.use.Bcell.NP$target)

pairLR.use.Bcell.PBS$interaction_ST = paste0(pairLR.use.Bcell.PBS$interaction_name, "_",
                                             pairLR.use.Bcell.PBS$source, "_",
                                             pairLR.use.Bcell.PBS$target)

common_B_int = intersect(pairLR.use.Bcell.PBS$interaction_ST, pairLR.use.Bcell.NP$interaction_ST)
common_B_int_net = net.Bcells[which(net.Bcells$interaction_ST %in% common_B_int ),]

dif_NP = pairLR.use.Bcell.NP[-which(pairLR.use.Bcell.NP$interaction_ST %in% common_B_int ),]
dif_PBS = pairLR.use.Bcell.PBS[-which(pairLR.use.Bcell.PBS$interaction_ST %in% common_B_int ),]

png(file = "chordpath_Ch7_B_allcommon.png",width =8, height =8, unit = "in", res = 400)
netVisual_chord_gene(object.list[[2]],
                     #sources.use = c(10),
                     #targets.use = 10,
                     slot.name = "netP", net = common_B_int_net, 
                     #pairLR.use = pairLR.use.up.Bcellsource[, "interaction_name", drop = F],
                     title.name = paste0("All B cell interactions in both NP and PBS"),
                     legend.pos.x = 1,
                     scale = T)
dev.off()

png(file = "chordpath_Ch7_B_NPonly.png",width =8.5, height =8, unit = "in", res = 400)
netVisual_chord_gene(object.list[[2]],
                     #sources.use = c(10),
                     #targets.use = 10,
                     slot.name = "netP", net = dif_NP, 
                     #pairLR.use = pairLR.use.up.Bcellsource[, "interaction_name", drop = F],
                     title.name = paste0("All B cell interactions in ", names(object.list)[2], " only"),
                     legend.pos.x = 1,
                     scale = T)
dev.off()

png(file = "chordpath_Ch7_B_PBSonly.png",width =8, height =8, unit = "in", res = 400)
netVisual_chord_gene(object.list[[1]],
                     #sources.use = c(10),
                     #targets.use = 10,
                     slot.name = "netP", net = dif_PBS, 
                     #pairLR.use = pairLR.use.up.Bcellsource[, "interaction_name", drop = F],
                     title.name = paste0("All B cell interactions in ", names(object.list)[1], " only"),
                     legend.pos.x = 1,
                     scale = T)
dev.off()

  
png(file = "chordpath_Ch7_B_NP_all.png",width =8, height =8, unit = "in", res = 400)
netVisual_chord_gene(object.list[[2]],
                     #sources.use = c(10),
                     #targets.use = 10,
                     slot.name = "netP", net = pairLR.use.Bcell.NP, 
                     #pairLR.use = pairLR.use.up.Bcellsource[, "interaction_name", drop = F],
                     title.name = paste0("All B cell interactions in ", names(object.list)[2]),
                     legend.pos.x = 1,
                     scale = T)
dev.off()

png(file = "chordpath_Ch7_B_PBS_all.png",width =8, height =8, unit = "in", res = 400)
netVisual_chord_gene(object.list[[1]],
                     #sources.use = c(10),
                     #targets.use = 10,
                     slot.name = "netP", net = pairLR.use.Bcell.PBS, 
                     #pairLR.use = pairLR.use.up.Bcellsource[, "interaction_name", drop = F],
                     title.name = paste0("All B cell interactions in ", names(object.list)[1]),
                     legend.pos.x = 1,
                     scale = T)
dev.off()

png(file = "chordpath_Ch7_B_NP_all.png",width =6, height =5, unit = "in", res = 400)
netVisual_chord_gene(object.list[[2]],
                     #sources.use = c(10),
                     #targets.use = 10,
                     slot.name = "netP", net = pairLR.use.up.Bcellsource, 
                     #pairLR.use = pairLR.use.up.Bcellsource[, "interaction_name", drop = F],
                     title.name = paste0("All B cell interactions in ", names(object.list)[2]),
                     legend.pos.x = 1,
                     scale = T)
dev.off()

png(file = "chordpath_Ch7_B_PBS.png_all",width =6, height =10, unit = "in", res = 400)
netVisual_chord_gene(object.list[[1]],
                     #sources.use = c(10),
                     #targets.use = 10,
                     slot.name = "netP", net = pairLR.use.down.Bcellsource, 
                     #pairLR.use = pairLR.use.up.Bcellsource[, "interaction_name", drop = F],
                     title.name = paste0("All B cell interactions in ", names(object.list)[1]),
                     legend.pos.x = 1,
                     scale = T)
dev.off()
cairo_pdf(file = "chordpath_Ch7_T_APC.pdf",width = 14, height = 10)
png(file = "chordpath_Ch7_T_APC.png",width = 14, height = 6, unit = "in", res = 400)
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  if ( i == 1){

    netVisual_chord_gene(object.list[[i]],
                         #sources.use = c(10),
                         #targets.use = 10,
                         slot.name = "netP", net = pairLR.use.down.Tcell, 
                         signaling = c("IL10", "TGFb", "CD80", "CD86", "MHC-II"),
                         title.name = paste0("Upregulated T cell interactions in ", names(object.list)[i]),
                         legend.pos.x = 170,
                         scale = T)
  }
  else{

    netVisual_chord_gene(object.list[[i]],
                         #sources.use = c(10),
                         #targets.use = c(10),
                         slot.name = "netP", net = pairLR.use.up.Tcell, 
                         signaling = c("IL10", "TGFb", "CD80", "CD86", "MHC-II"),
                         title.name = paste0("Upregulated T cell interactions in ", names(object.list)[i]),
                         legend.pos.x = 1,
                         scale = T)
  }
}
dev.off()




png(file = "chordgene_Ch7_IL10TGFb.png",width = 13, height = 6, unit = "in", res = 400)
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  if (i == 1){
    netVisual_chord_gene(object.list[[i]],
                         signaling = c("IL10", "TGFb"),
                         slot.name = "net", net = net.down,
                         title.name = paste0("IL10 and TGFb Signaling in ", names(object.list)[i]),
                         legend.pos.x = 170,
                         scale = T)
  }else{
  netVisual_chord_gene(object.list[[i]],
                       signaling = c("IL10", "TGFb"),
                       slot.name = "net", net = net.up,
                       title.name = paste0("IL10 and TGFb Signaling in ", names(object.list)[i]),
                       legend.pos.x = 10,
                       scale = T)
  }
}
dev.off()


netAnalysis_signalingRole_network(cellchat_PBS7, signaling = c( "IL10"), width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat_NP7, signaling = c( "TGFb"), width = 8, height = 2.5, font.size = 10)



## ----Top Gene changing over time per Cell Type in PBS mice--------------------------------------------------------------------------------------------
unique(int.SILP$sample)
Idents(int.SILP) = int.SILP$Cell.Type.TP.Tx

Idents(int.SILP) = "sample"

int.SILP$Cell.Type3.TP.Tx = paste0(int.SILP$Cell.types3, "_", int.SILP$TP, "_", int.SILP$Tx)
unique(int.SILP$Cell.Type3.TP.Tx)
Idents(int.SILP) = "Cell.Type3.TP.Tx"

Treg_7_PBS_NP = FindMarkers(
  int.SILP,
  ident.1 = "Treg_Ch7_PBS",
  ident.2 = "Treg_Ch7_NP", 
  only.pos = F
)
Treg_7_PBS_NP = Treg_7_PBS_NP[which(Treg_7_PBS_NP$p_val_adj <0.1),]

saveRDS(cellchat, file = "cellchat_merged_PBS7_NP7.rds")
