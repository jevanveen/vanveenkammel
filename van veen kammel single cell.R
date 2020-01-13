
library(tidyverse)
library(Seurat)
library(calibrate)
library(cowplot)

packageVersion("Seurat")

#get rid of annoying ggplot gray grid
theme_set(theme_classic())


### read raw UMI data ####

male1.data <- as.matrix(read.csv("~/Dropbox/van Veen dropbox/Single Cell/d701PC.csv", row.names = 1))
male2.data <- as.matrix(read.csv("~/Dropbox/van Veen dropbox/Single Cell/d710PC.csv", row.names = 1))
male3.data <- as.matrix(read.csv("~/Dropbox/van Veen dropbox/Single Cell/MRmale1PC.csv", row.names = 1))

female1.data <- as.matrix(read.csv("~/Dropbox/van Veen dropbox/Single Cell/d709PC.csv", row.names = 1))
female2.data <- as.matrix(read.csv("~/Dropbox/van Veen dropbox/Single Cell/MRfem1PC.csv", row.names = 1))
female3.data <- as.matrix(read.csv("~/Dropbox/van Veen dropbox/Single Cell/MRfem2PC.csv", row.names = 1))

male1 <- CreateSeuratObject(raw.data = male1.data, min.cells = 2, min.genes = 250, project = "m1_singles")
male2 <- CreateSeuratObject(raw.data = male2.data, min.cells = 2, min.genes = 250, project = "m2_singles")
male3 <- CreateSeuratObject(raw.data = male3.data, min.cells = 2, min.genes = 250, project = "m3_singles")

female1 <- CreateSeuratObject(raw.data = female1.data, min.cells = 2, min.genes = 250, project = "f2_singles")
female2 <- CreateSeuratObject(raw.data = female2.data, min.cells = 2, min.genes = 250, project = "f3_singles")
female3 <- CreateSeuratObject(raw.data = female3.data, min.cells = 2, min.genes = 250, project = "f4_singles")

male <- MergeSeurat(object1 = male1, object2 = male2, add.cell.id1 = "m1", add.cell.id2 = "m2", project = "m_singles")
table(male@meta.data$orig.ident)

male <- MergeSeurat(object1 = male, object2 = male3, add.cell.id2 = "m3", project = "m_singles")
female <- MergeSeurat(object1 = female1, object2 = female2, add.cell.id1 = "f1", add.cell.id2 = "f2", project = "f_singles")
female <- MergeSeurat(object1 = female, object2 = female3, add.cell.id2 = "f3", project = "f_singles")

male@meta.data$sex <- "Male"
female@meta.data$sex <- "Female"

Singles <- MergeSeurat(object1 = male, object2 = female, add.cell.id1 = "male", add.cell.id2 = "female", project = "sf1_singles")


#basic metrics ####
table(Singles@meta.data$orig.ident)
table(Singles@ident)
table(Singles@meta.data$sex)
head(x = Singles@cell.names)
median(Singles@meta.data$nGene)
mean(Singles@meta.data$nGene)
mean(Singles@meta.data$nUMI)
median(Singles@meta.data$nUMI)

# normalize data ####
mito.genes <- grep(pattern = "^mt-", x = rownames(x = Singles@data), value = TRUE)
percent.mito <- Matrix::colSums(Singles@raw.data[mito.genes, ])/Matrix::colSums(Singles@raw.data)
Singles <- AddMetaData(object = Singles, metadata = percent.mito, col.name = "percent.mito")
par(mfrow = c(1, 2))


VlnPlot(object = Singles, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, group.by = "sex")

GenePlot(object = Singles, gene1 = "nGene", gene2 = "percent.mito", do.hover = TRUE)

GenePlot(object = Singles, gene1 = "nUMI", gene2 = "nGene", do.hover = TRUE)

Singles <- FilterCells(object = Singles, subset.names = c("nGene", "percent.mito"), low.thresholds = c(500, 0), high.thresholds = c(Inf, .17))

Singles <- NormalizeData(object = Singles, normalization.method = "LogNormalize", scale.factor = 10000)

Singles <- ScaleData(object = Singles, vars.to.regress = c("nUMI", "percent.mito"))


#inspect gene expression between preps ####
#how many cells pass filtering? 261 females 269 males. 530 total.
table(Singles@meta.data$orig.ident)
table(Singles@meta.data$sex)


#examine IEG between preps
VlnPlot(object = Singles, features.plot = c("Fos", "Myc", "Jun"), nCol = 3, group.by = "orig.ident")
VlnPlot(object = Singles, features.plot = c("Egr1", "Arc", "Homer1"), nCol = 3, group.by = "orig.ident")

#neuro vs glia


#glia
VlnPlot(object = Singles, features.plot = c("Gfap", "Olig1"), nCol = 2, group.by = "orig.ident")
#neurons
VlnPlot(object = Singles, features.plot = c("Tubb3", "Nefl"), nCol = 2, group.by = "orig.ident")
#glutamatergic
VlnPlot(object = Singles, features.plot = c("Slc17a6", "Grin1", "Grin2b", "Gls", "Glul"), nCol = 3, group.by = "orig.ident")
#Gabaergic
VlnPlot(object = Singles, features.plot = c("Slc32a1", "Gabra1", "Gabbr1", "Gabbr2", "Gad1", "Gad2"), nCol = 3, group.by = "orig.ident")

#cluster cells ####

Singles <- FindVariableGenes(object = Singles, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.3, x.high.cutoff = 10, y.cutoff = 2.75, num.bin = 5)

length(x = Singles@var.genes)
print(x= Singles@var.genes)

Singles <- RunPCA(object = Singles, pc.genes = Singles@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PrintPCA(object = Singles, pcs.print = 1, genes.print = 50, use.full = FALSE)
VizPCA(object = Singles, pcs.use = 1:2)
PCAPlot(object = Singles, dim.1 = 1, dim.2 = 2, do.hover = TRUE)
PCHeatmap(object = Singles, pc.use = 1:6, cells.use = 700, num.genes = 50, do.balanced = TRUE, label.columns = FALSE)

PrintPCA(object = Singles, pcs.print = 1, genes.print = 200, use.full = FALSE)


Singles <- JackStraw(object = Singles, num.replicate = 100, do.print = FALSE)
JackStrawPlot(object = Singles, PCs = 1:20)
PCElbowPlot(object = Singles)

Singles <- FindClusters(object = Singles, reduction.type = "pca", dims.use = 1:8, resolution = 1.4, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
PrintFindClustersParams(object = Singles)
VlnPlot(object = Singles, features.plot = c("Apoe", "Pomc"))

#tsne with all cells, in fig. supp 2
Singles <- RunTSNE(object = Singles, dims.use = 1:8, do.fast = FALSE)
TSNEPlot(object = Singles, pt.size = .75)
save_plot(filename = "~/box/scHypo paper/2020tsne.pdf", plot = p1, base_asp = 1.2, base_height = 3)


#remove apoe and pomc clusters, remake tsne
SinglesVMH <- SubsetData(object = Singles, ident.remove = c("6", "7"), do.clean = F)
SinglesVMH <- RunTSNE(object = SinglesVMH, dims.use = 1:8, do.fast = FALSE)
TSNEPlot(object = SinglesVMH, pt.size = .75)
save_plot(filename = "~/box/scHypo paper/2020vmhtsne.pdf", plot = p1, base_asp = 1.2, base_height = 3)


singles.markers <- FindAllMarkers(Singles, only.pos = FALSE, min.pct = 0.1, logfc.threshold = .1)
singles.markers.sig <- filter(singles.markers, p_val_adj < .05)
singles.markers.logfc <- filter(singles.markers, avg_logFC > 1)

#top markers from each cluster
DoHeatmap(SinglesVMH, genes.use = c("Hpcal1", "Fezf1", "Pcp4", "Tac1", "Ifi27l2a", "Serpine2", "Pdyn", "Tmem35a", "Mt1", "Gal", "Penk", "Synpr", "Sst", "Cartpt", "Nts", "Rprm", "Crym", "Syp"))


#cluster tree
BuildClusterTree(Singles)
BuildClusterTree(SinglesVMH)
