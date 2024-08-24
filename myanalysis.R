library(patchwok)
source("functs.R")
load("All.Rdata")
library(ggpubr)

sce_cca <- sce_cca %>% 
  RunUMAP(dims = pc.num, reduction = "cca")



DimPlot(sce_cca, arg=1, pt.size=0.5, label.size = 5) + theme(aspect.ratio = 1) + scale_col_manual(values=cols25())



featureplot_celltype(sce_cca, "Features_sce_cca","tsne")

featureplot_celltype(sce_cca, "Features_sce_cca","umap")


celltype_VlnPlot(sce_cca, "vln.pdf", overal_vln)


celltype_redblue_dotplot(sce_cca,
                         filename = "dotplot.pdf",
                         features=overal_vln, 
                         w = 13, h=7.5)

 
bone_cca <- bone_cca %>% 
  RunUMAP(dims = pc.num, reduction = "cca")

reduct = "umap"


DimPlot(bone_cca, pt.size=0.5, label.size = 5,  reduction=reduct) + theme(aspect.ratio = 1) + scale_col_manual(values=unname(cols25()))


colCls <- colorRampPalette(brewer.pal(num, "Paired"))(nClust)


pbmc.markers <- read.csv("markers.csv", row.names = 1)

markers <- pbmc.markers %>%
  slice_max(n = 5, order_by = vg)

Htmap(bone_cca, f = markers$gene, grp = pal, hsize = 4, raster = F) 


DimPlot(bone_cca, reduction="umap", pt.size=0.5, label.size = 5, split.by = "orig", group.by = 'typ') + scale_col_manual(values=unname(cols25()))


DimPlot(bone_cca, reduction="umap", pt.size=0.5, label.size = 5, split.by = "orig") + scale_col_manual(values=unname(cols25()))

bone_cca@active.ident


levels(bone_cca)


DimPlot(bone_cca,   pt.size=0.5, label.size = 3, repel = T, split.by = "orig", shape.by = "orig", label.box = TRUE) + scale_col_manual(values=pal) 

DimPlot(bone_cca,  pt.size=0.5, label.size = 3, repel = T, shape.by = "orig") + scale_col_manual(values=pal)

# 
p1 <- DimPlot(sce_ctrl, pt.size=0.5, label.size = 3, repel = T, label.box = TRUE) + scale_col_manual(values=pal) + theme(aspect.ratio = 1, legend.position="none") 


p2 <- DimPlot(sce_oiv, pt.size=0.5, label.size = 3, repel = T, label.box = TRUE) + scale_col_manual(values=pal) + theme(aspect.ratio = 1, legend.position="none") 

ggsave(plot = p1, filename = "see.pdf", w=4, h=5)


ggsave(plot = p1/p2, filename = "split.pdf", w=4, h=9)


bone_cca$cell_type <- ident(bone_cca)

fortime_cell$cell_type <- ident(fortime_cell)

df_sta <- df %>% 
 group_by(Group) %>%
 mutate(Percent = Count/sum(Count)*100)
df_sta


d_wide <- dcast(df, cell - Group, value.var = "counts")



p1 <- pp +
 geom_col(width = 0.6)


p1 + theme(axis.title.x=element_blank(), axis.text.x=element_text(size = 15)) 


p1 <- FeaturePlot(object = object, features = "score1", split.by = "orig") & theme(aspect.ratio = 1)


ggsave(plot = p1, "Score.pdf",w=14,h=4)

ident(bone_cca)

p <- list()

for (i in aoc) { 
p[[i]] <- 
  oc %>%
  Vlnplot(features=i,group.by = cel) +
 theme(axis.title.x = element_blank(), legend.position = "none",
          plot.title = element_text(face="bold", size = 20),
          strip.text.x = element_text(size=5),
          strip.background = element_rect(color = "black", fill = "white", linetype = "solid")) 
}
  


ggarrange(plotlist = p)
dev.off()

osteo_cell



DimPlot(my_osteo, reduction="umap", pt.size=1, label.size = 5, shape.by = "orig", split.by = "orig", repel = T) + scale_col_manual(values=unname(cols25()))


dev.off()





