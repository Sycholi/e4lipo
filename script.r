rm(list = ls());gc();setwd('~/r_repo/USP14')
source('~/r_repo/start_up_toolbox.r');use_condaenv('bioinfo')

seu = readRDS('seu.rds')
seu.nkt = readRDS('seu.nkt.rds')

# dir.create('figures')
# UMAP_Type ----
seu$type %<>% factor(levels = c('B_Plasma','NK_T','Endothelial','Epithelial','Fibroblast','Mast','Mono_Mac','Neutrophil','DC','mDC','pDC'))
Idents(seu) = 'type'
umap1 = DimPlot(seu, cols = 'Spectral');umap1
ggsave(umap1,filename = 'figures/UMAP_Type.pdf',height = 2000,width = 2400,units = 'px',dpi = 300)
# 72949 cells

# Annotate_Dotplot_Type ----
seu$type %<>% factor(levels = c('NK_T','Neutrophil','Mono_Mac','Mast','Fibroblast','Epithelial','Endothelial','pDC','mDC','DC','B_Plasma'))
immune_cellmarker = list(B_Plasma = c('Cd79a','Cd79b','Jchain'),
                         DC = c('Clec10a','Cd86'),
                         mDC = c('Ccr7','Cd40'),
                         pDC = c('Siglech','Ccr9'),
                         Endo = c('Pecam1','Flt1'),
                         Epithelial = c('Epcam','Krt19','Afp','Alb','Aldh1a1'),
                         Fibroblast = c('Col1a2','Cxcl12','Col3a1','Acta2','Col1a1'),
                         Mast = c('Gata2','Ms4a2','Cpa3'),
                         Mono_Mac = c('Csf1r','C1qa','C1qb'),
                         Neutrphil = c('S100a8','S100a9','Csf3r'),
                         NK_T = c('Cd3d','Cd3e','Cd3g','Ncr1','Nkg7')
)
immune_dotplot = DotPlot(seu, features = immune_cellmarker, group.by = 'type') +
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(hjust = 1,vjust = 0.5,angle = 90)) +
  labs(x = NULL, y = NULL) + guides(size = guide_legend(order =  3)) +
  scale_color_gradientn(values = seq(0, 1, 0.2), 
                        colours = c('#2B83BA','#ABDDA4','#FDAE61','#D7191C'))
immune_dotplot
ggsave(immune_dotplot,filename = 'figures/Annotate_Dotplot_Type.pdf',
       height = 1200,width = 3300,units = 'px',dpi = 300)

# FeaturePlot_Type ----
fp = FeaturePlot(seu, c('Cd79a','Clec10a','Ccr7','Siglech','Flt1','Epcam',
                        'Col1a2','Gata2','Csf1r','S100a8','Cd3d','Nkg7'),
                 cols = c('grey90','#A50026'));fp
ggsave(fp,filename = 'figures/FeaturePlot_Type.pdf',
       height = 2500,width = 4000,units = 'px',dpi = 300)


# Boxplot_ROGUE_Type ----
expr = GetAssayData(seu, assay = 'SCT')
rogue.res = rogue(expr, 
                  labels = seu$type,
                  samples = seu$orig.ident,
                  platform = "UMI",
                  span = 0.5)
colnames(rogue.res) %>% sort() %>% paste0(collapse = "','")
df = pivot_longer(rogue.res,
                  cols = c('B_Plasma','DC','mDC','pDC','Endothelial','Epithelial','Fibroblast','Mast','Mono_Mac','Neutrophil','NK_T'),
                  names_to = 'Clusters',
                  values_to = 'ROGUE.Purity')
paste0(rownames(rogue.res), collapse = "','")
df$treat = rep(c('NC1','NC2','NC3','TR1','TR2','TR3'), each = 11)
df$Clusters %<>% as.factor()
df$treat %<>% as.factor()
boxplot.rogue = ggplot(data = df, aes(x = Clusters, y = ROGUE.Purity)) +
  geom_boxplot(aes(group = Clusters, colour = Clusters)) + 
  ylim(c(0,1)) +
  scale_color_manual(values = c("#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B",
                                '#F6CB1D',"#FFFFBF","#E6F598","#ABDDA4","#8AD088",
                                '#66C2A5','#017D6F','#19663E',"#3288BD",'#5E4FA2',
                                "#313695",'#A5678E',
                                "#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B",
                                '#F6CB1D',"#FFFFBF","#E6F598","#ABDDA4","#8AD088",
                                '#66C2A5','#017D6F','#19663E',"#3288BD",'#5E4FA2',
                                "#313695",'#A5678E')) +
  geom_jitter(aes(colour = treat), shape = 16, position = position_jitter(0.2)) + 
  theme_classic() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(hjust = 1,angle = 45))
boxplot.rogue

ggsave(boxplot.rogue,
       filename = 'figures/Boxplot_ROGUE_Type.pdf',
       height = 1500,
       width = 2000,
       units = 'px',
       dpi = 300)
# Barplot_Sample ----
df = table(tibble(clusters = seu$Sample, group = seu$type) %>% group_by(group)) %>% as_tibble()
df$group %<>% factor(levels = c('B_Plasma','NK_T','Endothelial','Epithelial','Fibroblast','Mast','Neutrophil','Mono_Mac','DC','mDC','pDC'))
barplot1 = ggplot(data = df, aes(x = clusters, y = n, fill = group)) +
  geom_bar(stat = "identity", position = 'fill') +
  scale_fill_manual(values = brewer.pal(11,'Spectral'), name = 'Cluster') + 
  labs(y = 'Cell Type Ratio', x = 'Datasets') + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank())
barplot1
ggsave(barplot1, filename = 'figures/Barplot_Sample.pdf',
       height = 1500, width = 1500, units = 'px', dpi = 300)

# Barplot_Cell_Type ----
md = seu@meta.data
df = table(tibble(clusters = seu$type, sample = seu$orig.ident) %>%
           group_by(sample)) %>% as_tibble() %>% arrange(clusters)
df$group = rep(c('NC','NC','NC','TR','TR','TR'),11)
df_sum = tibble(sample = seu$orig.ident) %>% group_by(sample) %>% table() %>% as_tibble() 
df$sum = rep(df_sum$n, 11)
df$percent = df$n / df$sum
df$color = rep(1:22, each = 3) %>% as.factor()
df$clusters %<>% factor(levels = c('Neutrophil','NK_T','Mono_Mac','Epithelial','B_Plasma','DC','Endothelial','Fibroblast','mDC','pDC','Mast'))

bp = ggplot(df, aes(x = clusters, y = percent, group = group, fill = color, colour = color)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(0.8), width = 0.7) +
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = 'grey80',
               width = 0.1, position = position_dodge(0.8)) + 
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                     method = "t.test", 
                     method.args = list(alternative = "two.sided")) + 
  scale_fill_manual(values = c(c("#9E0142",'white',"#D53E4F",'white',"#F46D43",'white',"#FDAE61",'white',"#FEE08B",'white',
                                 "#FFFFBF",'white',"#E6F598",'white',"#ABDDA4",'white',"#66C2A5",'white',"#3288BD",'white',
                                 '#313695','white',"#5E4FA2",'white','#A5678E','white'))) + 
  scale_color_manual(values = c('#9E0142','#9E0142','#D53E4F','#D53E4F','#F46D43','#F46D43','#FDAE61','#FDAE61','#FEE08B','#FEE08B',
                                '#FFFFBF','#FFFFBF','#E6F598','#E6F598','#ABDDA4','#ABDDA4','#66C2A5','#66C2A5','#3288BD','#3288BD',
                                '#313695','#313695','#5E4FA2','#5E4FA2','#A5678E','#A5678E')) + 
  scale_x_discrete("Cell Type")+ 
  ylab('count') + NoLegend() +
  theme(panel.background = element_blank(), axis.line.y = element_line()) + 
  geom_hline(yintercept = 0)
bp

ggsave(bp, filename = 'figures/Barplot_Cell_Type.pdf',
       width = 3000, height = 1800, units = 'px',dpi = 300)

# UMAP_NK_T ----
umap2 = DimPlot(seu.nkt, group.by = 'subtype', pt.size = 0.8,
                cols = c('#9E0142','#D53E4F',"#F46D43","#FDAE61","#FEE08B",
                         '#F6CB1D',"#FFFFBF"));umap2
ggsave(umap2,filename = 'figures/UMAP_NK_T.pdf',height = 2000,width = 2400,units = 'px',dpi = 300)

# Barplot_NK_T ----
df = table(tibble(clusters = seu.nkt$Sample, group = seu.nkt$subtype) %>% group_by(group)) %>% as_tibble()
barplot2 = ggplot(data = df, aes(x = clusters, y = n, fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c('#9E0142','#D53E4F',"#F46D43","#FDAE61","#FEE08B",
                               '#F6CB1D',"#FFFFBF"), name = 'Cluster') + 
  labs(y = 'Cell Type Ratio', x = 'Datasets') + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank())
barplot2
ggsave(barplot2, filename = 'figures/Barplot_NK_T.pdf',
       height = 1500, width = 1500, units = 'px', dpi = 300)

# Annotate_Dotplot_NK_T ----
seu.nkt$subtype %>% unique() %>% paste0(collapse = "','")
seu.nkt$subtype %<>% factor(levels = c('Tn','CD4T','CD8T','Treg','NKT','NK'))
Idents(seu.nkt) = 'subtype'
immune_cellmarker = list(T_lineage = c('Cd3g','Cd4', 'Cd8a','Il2ra'), 
                         Cytotoxic = c('Gzmb','Gzmk','Prf1','Xcl1','Nkg7'),
                         Regulatory = c('Tnfrsf4','Batf','Ltb','Foxp3'),
                         Inhibitory = c('Pdcd1','Lag3','Tigit','Havcr2', 'Tox','Ctla4'),
                         Naive = c('Il7r','Tcf7','Lef1','Cd44'),
                         NK = c('Ncr1','Klrb1c'))
immune_dotplot = DotPlot(seu.nkt, features = immune_cellmarker) +
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(hjust = 1,vjust = 0.5,angle = 90)) +
  labs(x = NULL, y = NULL) + guides(size = guide_legend(order =  3)) +
  scale_color_gradientn(values = seq(0, 1, 0.2), 
                        colours = c('#2C7BB6', '#ABD9E9', '#FDAE61', '#D7191C'))
immune_dotplot
ggsave(immune_dotplot, filename = 'figures/Annotate_Dotplot_NK_T.pdf',
       height = 1000, width = 2200, units = 'px', dpi = 300)


# UMAP_Subtype----
seu$subtype %>% unique() %>% sort() %>% paste0(collapse = "','")
seu$subtype %<>% factor(levels = c('B_Plasma','CD4T','CD8T','NK','NKT','Tn','Treg','DC','mDC','pDC','Endothelial','Epithelial','Fibroblast','Mast','Mono_Mac','Neutrophil'))
umap3 = DimPlot(seu, group.by = 'subtype',
                cols = c("#9E0142","#D53E4F","#F46D43","#FDAE61",'#FEE08B',
                         "#FFFFBF","#E6F598","#ABDDA4","#8AD088",'#66C2A5',
                         '#017D6F','#19663E','#313695',"#5E4FA2",'#A5678E',
                         "#3288BD"));umap3

ggsave(umap3,filename = 'figures/UMAP_Subtype.pdf',height = 2000,width = 2400,units = 'px',dpi = 300)

# Barplot_Subtype ----
df = table(tibble(clusters = seu$Sample, group = seu$subtype) %>% group_by(group)) %>% as_tibble()
barplot1 = ggplot(data = df, aes(x = clusters, y = n, fill = group)) +
  geom_bar(stat = "identity", position = 'fill') +
  scale_fill_manual(values = c("#9E0142","#D53E4F","#F46D43","#FDAE61",'#FEE08B',
                               "#FFFFBF","#E6F598","#ABDDA4","#8AD088",'#66C2A5',
                               '#017D6F','#19663E','#313695',"#5E4FA2",'#A5678E',
                               "#3288BD"), name = 'Cluster') + 
  labs(y = 'Cell Type Ratio', x = 'Datasets') + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank())
barplot1
ggsave(barplot1, filename = 'figures/Barplot_Subtype.pdf',
       height = 1500, width = 1500, units = 'px', dpi = 300)

# Barplot_Cell_SubType ----
md = seu@meta.data
df = table(tibble(clusters = seu$subtype, sample = seu$orig.ident) %>%
             group_by(sample)) %>% as_tibble() %>% arrange(clusters)
df$group = rep(c('NC','NC','NC','TR','TR','TR'),16)
df_sum = tibble(sample = seu$orig.ident) %>% group_by(sample) %>% table() %>% as_tibble() 
df$sum = rep(df_sum$n, 16)
df$percent = df$n / df$sum
df$color = rep(1:32, each = 3) %>% as.factor()
df$clusters %<>% factor(levels = c('Neutrophil','Mono_Mac','CD8T','Epithelial','B_Plasma','CD4T','DC','Fibroblast','NK','Tn','Endothelial','NKT','Treg','mDC','pDC','Mast'))

bp = ggplot(df, aes(x = clusters, y = percent, group = group, fill = color, colour = color)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(0.8), width = 0.7) +
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = 'grey80',
               width = 0.1, position = position_dodge(0.8)) + 
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                     method = "t.test", 
                     method.args = list(alternative = "two.sided")) + 
  scale_fill_manual(values = c(c("#9E0142",'white',"#D53E4F",'white',"#F46D43",'white',"#FDAE61",'white','#FEE08B','white',
                                 "#FFFFBF",'white',"#E6F598",'white',"#ABDDA4",'white',"#8AD088",'white','#66C2A5','white',
                                 '#017D6F','white','#19663E','white','#313695','white',"#5E4FA2",'white','#A5678E','white',
                                 "#3288BD",'white'))) + 
  scale_color_manual(values = c('#9E0142','#9E0142','#D53E4F','#D53E4F','#F46D43','#F46D43','#FDAE61','#FDAE61','#FEE08B','#FEE08B','#FFFFBF','#FFFFBF','#E6F598','#E6F598','#ABDDA4','#ABDDA4','#8AD088','#8AD088','#66C2A5','#66C2A5','#017D6F','#017D6F','#19663E','#19663E','#313695','#313695','#5E4FA2','#5E4FA2','#A5678E','#A5678E','#3288BD','#3288BD')) + 
  scale_x_discrete("Cell Type")+ 
  ylab('count') + NoLegend() +
  theme(panel.background = element_blank(), axis.line.y = element_line()) + 
  geom_hline(yintercept = 0)
bp

ggsave(bp, filename = 'figures/Barplot_Cell_Subtype.pdf',
       width = 4000, height = 1500, units = 'px',dpi = 300)

# Boxplot_ROGUE_Subtype ----
expr = GetAssayData(seu, assay = 'SCT')
rogue.res = rogue(expr, 
                  labels = seu$subtype,
                  samples = seu$orig.ident,
                  platform = "UMI",
                  span = 1.0)
colnames(rogue.res) %>% sort() %>% paste0(collapse = "','") 
df = pivot_longer(rogue.res,
                  cols = c('B_Plasma','CD4T','CD8T','NK','NKT','Tn','Treg','DC','mDC','pDC','Endothelial','Epithelial','Fibroblast','Mast','Mono_Mac','Neutrophil'),
                  names_to = 'Clusters',
                  values_to = 'ROGUE.Purity')
paste0(rownames(rogue.res), collapse = "','")
df$treat = rep(c('NC1','NC2','NC3','TR1','TR2','TR3'), each = 16)
df$Clusters %<>% as.factor()
df$treat %<>% as.factor()
boxplot.rogue = ggplot(data = df, aes(x = Clusters, y = ROGUE.Purity)) +
  geom_boxplot(aes(group = Clusters, colour = Clusters)) + 
  ylim(c(0,1)) +
  scale_color_manual(values = c("#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B",
                                '#F6CB1D',"#FFFFBF","#E6F598","#ABDDA4","#8AD088",
                                '#66C2A5','#017D6F','#19663E',"#3288BD",'#5E4FA2',
                                "#313695",'#A5678E',
                                "#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B",
                                '#F6CB1D',"#FFFFBF","#E6F598","#ABDDA4","#8AD088",
                                '#66C2A5','#017D6F','#19663E',"#3288BD",'#5E4FA2',
                                "#313695",'#A5678E')) +
  geom_jitter(aes(colour = treat), shape = 16, position = position_jitter(0.2)) + 
  theme_classic() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(hjust = 1,angle = 45))
boxplot.rogue

ggsave(boxplot.rogue,
       filename = 'figures/Boxplot_ROGUE_Subtype.pdf',
       height = 1500,
       width = 2000,
       units = 'px',
       dpi = 300)
# Heatmap_ROE_Subtype----
library(Startrac)
Roe = calTissueDist(seu@meta.data,
                    byPatient = F,
                    colname.cluster = "subtype", # 不同细胞亚群
                    colname.patient = "Sample", # 不同样本
                    colname.tissue = "Group", # 不同组织
                    method = "chisq", # "chisq", "fisher", and "freq" 
                    min.rowSum = 0) 
Roe = data.frame(NC = Roe[,1], TR = Roe[,2])

p1 = pheatmap(Roe, 
              cluster_rows = T,
              display_numbers = T,
              treeheight_row = 0,
              treeheight_col = 0,
              border_color = NA)
ggsave(p1,
       filename = 'figures/Heatmap_ROE_subtype.pdf',
       width = 1000,
       height = 1500,
       units = 'px')


# Heatmap_GSEA_Subtype ----
mdb = msigdbr(species = 'Mus musculus', category = 'H')
fgsea_sets = mdb %>% split(x = .$gene_symbol, f= .$gs_name)
gsea.list = list()
gsea.res.list = list()
seu %<>% JoinLayers(assay = 'RNA')
Celltypes = unique(seu$subtype)

for(i in 1:length(Celltypes)){
  cell_type = Celltypes[i]
  print(paste0('Processing: ', cell_type, ' ...'))
  seu.obj = subset(seu, subtype == cell_type %>% as.character())  # tune
  Idents(seu.obj) = 'Group'
  saveRDS(seu.obj, file = paste0('subcluster_Obj/', cell_type, '.rds'))
  markers = seu.obj %>% 
    PrepSCTFindMarkers() %>% 
    FindMarkers(ident.1 = 'TR',
                ident.2 = 'NC')
  write.csv(markers, file = paste0('markers/', cell_type, '_markers.csv'))
  # markers %<>% filter(p_val_adj <= 0.05)
  markers.gsea = markers %>%
    rownames_to_column('gene') %>%
    arrange(desc(avg_log2FC)) %>%
    dplyr::select(gene, avg_log2FC)
  markers.gsea = markers.gsea[!duplicated(markers.gsea$avg_log2FC),]
  markers.gsea %<>% deframe()
  fgseaRes = fgsea(fgsea_sets, stats = markers.gsea)
  fgseaResTidy = fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES)) %>%
    dplyr::select(-leadingEdge, -ES) %>%
    arrange(padj)
  fgseaResTidy$Group = cell_type %>% as.character()
  gsea.list[[i]] = fgseaResTidy
  gsea.res.list[[i]] = fgseaRes
  write.csv(fgseaResTidy, file = paste0('enrich/', cell_type, '_GSEA.csv'))
  saveRDS(fgseaRes, 'fgseaRes.rds')
}
gsea.all = gsea.list[[1]] %>% dplyr::select(Group, pathway, NES) %>% arrange(pathway)
for(i in 2:length(gsea.list)){
  gsea.all = rbind(gsea.all, gsea.list[[i]] %>% dplyr::select(Group, pathway, NES) %>% arrange(pathway))
}
gsea.all %<>% pivot_wider(names_from = Group, values_from = NES)
gsea.all %<>% column_to_rownames('pathway')
df.hp = as.data.frame(scale(gsea.all));df.hp = df.hp[,sort(colnames(df.hp))]
df.hp[df.hp >= 1.7] = 1.7
df.hp[df.hp <= -1.7] = -1.7

p1 = pheatmap(df.hp,
              # gaps_row = 1:(nrow(df.hp)-1),
              # gaps_col = 1:(ncol(df.hp)-1),
              border_color = 'white',
              treeheight_row = 0,
              treeheight_col = 0,
              cluster_cols = T,
              cluster_rows = T,
              angle_col = 45,
              cellwidth = 15,
              cellheight = 10,
              width = 30);p1
ggsave(p1, filename = 'figures/Heatmap_GSEA_Subtype.pdf',
       height = 3000, width = 3000, units = 'px', dpi = 300)


# BubblePlot_CellChat ----
p1 = netVisual_bubble(cellchat, 
                      sources.use = c(1:2,4:16), 
                      targets.use = 3,
                      signaling = c('CCL','CXCL'),
                      comparison = c(2, 1), 
                      angle.x = 45);p1
p2 = netVisual_bubble(cellchat, 
                      sources.use = c(1:11,13:16), 
                      targets.use = 12,
                      signaling = c('CCL','CXCL'),
                      comparison = c(2, 1), 
                      angle.x = 45);p2
p3 = netVisual_bubble(cellchat, 
                      sources.use = c(1:12,14:16), 
                      targets.use = 13,
                      signaling = c('CCL','CXCL'),
                      comparison = c(2, 1), 
                      angle.x = 45);p3

ggsave(p1, filename = 'figures/BubblePlot_CellChat.pdf',
       width = 2500, height = 1500, units = 'px', dpi = 300)

# VlnPLot_CXCL ----
v1 = VlnPlot(seu, 
             c('Cxcl16','Cxcr6','Prf1','Fas','Fasl','Ifng','Tnf','Gzmb','Lag3','Pdcd1','Havcr2'),
             stack = T,
             split.by = 'Group',
             cols = c("#4575B4","#A50026"));v1
ggsave(v1,
       filename = 'figures/VlnPlot_CXCL.pdf',
       width = 1000,
       height = 2000,
       units = 'px',
       dpi = 300)




# HierarchyPlot_CellChat ----
par(mfrow=c(1,2))
pdf('figures/HierarchyPlot_CellChat.pdf',
    height = 5, width = 10)
gg1 = netVisual_aggregate(cellchat.NC, signaling = 'CXCL', sources.use = 7, targets.use = c(2,3,12,13,15,16))
gg2 = netVisual_aggregate(cellchat.TR, signaling = 'CXCL', sources.use = 7, targets.use = c(2,3,12,13,15,16))
dev.off()

# ChordPlot_CellChat ----
par(mfrow=c(1,2))
pdf('figures/ChordPlot_CellChat.pdf',
    height = 5, width = 10)
netVisual_individual(cellchat.TR, signaling = 'CXCL', pairLR.use = LR.show, layout = "chord")
netVisual_individual(cellchat.NC, signaling = 'CXCL', pairLR.use = LR.show, layout = "chord")
dev.off()
