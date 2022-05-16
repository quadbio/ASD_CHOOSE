source('~/scripts/single_cell/celltype.R')
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/cell_density_diff.R')
source('~/scripts/perturbator/de.R')
source('~/scripts/perturbator/enrichment.R')
source('~/scripts/perturbator/ko_inference.R')
source('~/scripts/perturbator/plots.R')

library(Pando)

select <- dplyr::select

setwd('~/projects/ASD')


#### Read required data ####
crop <- read_rds('data/objects/CROP/asd_chong_annot_withguides_v3.3modules_srt.rds')
crop_ge <- read_rds('data/objects/CROP/asd_chong_annot_wg_ge_v3.1noastro_srt.rds')
cellrank_ge_meta <- read_tsv('data/velocity/ge_opc_cellrank_probs.tsv')
crop_muo_wt <- read_rds('data/objects/CROP_MUO/ASD_CROP_MUO_wt_pando_v2net_srt.rds')

sfari_genes <- read_csv('data/SFARI-Gene_genes_09-02-2021release_11-04-2021export.csv')
tfs <- read_tsv('~/resources/DB/animal_tfdb/tf_human.tsv')

sfari_genes_use <- sfari_genes %>% 
    filter(`gene-score`<=3)

ge_meta <- crop_ge@meta.data %>% as_tibble(rownames='cell') %>% 
    inner_join(cellrank_ge_meta, by=c('cell'='CellID')) %>% 
    left_join(as_tibble(crop_ge[['umap']]@cell.embeddings, rownames='cell')) %>% 
    select(!contains('GO'))

crop_ge <- subset(crop_ge, cells=ge_meta$cell)

crop_ge[['umap']]@cell.embeddings[,1] <- -crop_ge[['umap']]@cell.embeddings[,1]

markers_ge <- c('DLX2', 'DLX5', 'OLIG2', 'PDGFRA', 'PLP1', 'OLIG1')
feature_plot(crop_ge, features=markers_ge, order=T) & 
    article_text() &
    theme_void() & no_legend() &
    theme(plot.title=element_text(face='plain', size=6))
ggsave('plots/paper/fig4_ge_feature_plots_umap.png', width=8, height=8, unit='cm')

#### Check OLIG + DLX expression ####
crop_ge_cc <- subset(crop, lineage=='inhibitory')
crop_ge_cc_meta <- crop_ge_cc@meta.data %>% as_tibble(rownames='cell')

expr_mat <- crop_ge_cc[['RNA']]@data[c('OLIG2', 'DLX2'), ]

expr_df <- expr_mat %>% t() %>% 
    as_tibble(rownames='cell') %>% 
    inner_join(crop_ge_cc_meta) 

plot_df <- expr_df %>% 
    mutate(cat=factor(case_when(
        OLIG2 > 0 & DLX2 > 0 ~ 'OLIG2+DLX2',
        OLIG2 == 0 & DLX2 > 0 ~ 'DLX2',
        OLIG2 > 0 & DLX2 == 0 ~ 'OLIG2',
        T ~ 'none'
    ), levels=c('none', 'DLX2', 'OLIG2+DLX2', 'OLIG2'))) %>% 
    filter(gRNA%in%c('ARID1B', 'Control2')) %>% 
    filter(celltype_cl_refined=='vRG')

color_vec <- c('#158B76', '#4888A2', '#276775', 'darkgrey')
names(color_vec) <- c('OLIG2', 'DLX2', 'OLIG2+DLX2', 'none')

ggplot(plot_df, aes(gRNA, fill=cat)) +
    geom_bar(position='fill') +
    scale_fill_manual(values=color_vec) +
    # scale_x_discrete(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    rotate_x_text(40) +
    theme(
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank()
    ) +
    no_legend() +
    labs(y='Proportion of vRG cells')
ggsave('plots/paper/fig4_OLIG_DLX_proportion_bar.pdf', width=2.4, height=3.5, units='cm')



#### Get GRN ####
grn_features <- NetworkFeatures(crop_muo_wt, network='glm_bidir_GREAT')
grn_module_features <- NetworkModules(crop_muo_wt, network='glm_bidir_GREAT')@features
grn_modules_pos <- grn_module_features$genes_pos
grn_modules_neg <- grn_module_features$genes_neg
grn_modules_full <- list_merge(grn_modules_pos, !!!grn_modules_neg)
grn_modules <- NetworkModules(crop_muo_wt, network='glm_bidir_GREAT')@meta


#### ARID1B GRN #####
### -> ToDo coming soon



#### Cellrank viz ####
cr_space <- cellrank_ge_meta %>% 
    column_to_rownames('CellID') %>% as.matrix() 

lin_order <- c('to_lge2', 'to_cge', 'to_lge', 'to_opc')
probs <- cr_space[,lin_order] 
probs <- probs / colMeans(probs)
probs <- probs / rowSums(probs)

angle_vec = seq(0, 2*pi, length.out=5)[1:4]
angle_vec_sin = cos(angle_vec)
angle_vec_cos = sin(angle_vec)

x = rowSums(t(apply(probs, 1, function(x)x*angle_vec_sin)))
y = rowSums(t(apply(probs, 1, function(x)x*angle_vec_cos)))

ge_meta$C1 <- x
ge_meta$C2 <- y


circle_dr <- cbind(x,y)
colnames(circle_dr) <- NULL

crop_ge[['circular']] <- CreateDimReducObject(circle_dr, key='CR_')

p1 <- ggplot(ge_meta, aes(C1, C2, color=celltype_cl_refined)) +
    geom_point(size=1, alpha=0.8) +
    scale_color_manual(values=asd_celltype_colors2) +
    theme_void()

p2 <- ggplot(arrange(ge_meta, gRNA=='ARID1B'), aes(C1, C2, color=gRNA=='ARID1B')) +
    geom_point(size=1, alpha=0.8) +
    scale_color_manual(values=c('grey', 'black')) +
    theme_void()

p3 <- ggplot(arrange(ge_meta, gRNA=='Control2'), aes(C1, C2, color=gRNA=='Control2')) +
    geom_point(size=1, alpha=0.8) +
    scale_color_manual(values=c('grey', 'black')) +
    theme_void()

p1 + p2 + p3
ggsave('plots/paper/fig4_ge_circular.png', width=20, height=6)


#### UMAPs ####
dim_plot(crop_ge, group.by='celltype_cl_refined') +
    scale_color_manual(values=asd_celltype_colors2)

dim_plot(crop_ge, group.by='terminal_states', label=T) 

feature_plot(crop_ge, features=c('PDGFRA', 'PLP1', 'LHFPL3', 'NKX2-2', 'CSPG4'), order=T, pt.size=0.8)

feature_plot(crop_ge, features=c('OLIG2', 'DLX2', 'OLIG2', 'NR2F2',  'MEIS2', 'DLX5', 'PDGFRA', 'PLP1'), reduction='circular', order=T)
ggsave('plots/paper/fig4_ge_features_circular.png', width=10, height=10)


feature_plot(crop_ge, features=markers_ge, order=T, reduction='circular') & 
    article_text() &
    theme_void() & no_legend() &
    theme(plot.title=element_text(face='plain', size=6))
ggsave('plots/paper/fig4_ge_feature_plots_circular.png', width=8, height=12, unit='cm')



#### Density plots on circular for ARID1B ####
knn <- RANN::nn2(Embeddings(crop_ge, 'circular')[,1:2], k = 100)
knngraph <- FindNeighbors(Embeddings(crop_ge, 'circular')[,1:2], k.param = 100)

knn <- RANN::nn2(probs, k = 100)
knngraph <- FindNeighbors(probs, k.param = 100)

grna <- 'ARID1B'

knn_enrich <- get_knn_cmh_enrichment(knn, crop_ge$gRNA == grna, crop_ge$gRNA == 'Control2', stratums = crop_ge$library)
smooth_enrich <- random_walk_with_restart(knn_enrich, knn_mat = as(knngraph$snn, 'Matrix'), num_rounds = 100, alpha = 0.1)

plot_df <- ge_meta %>% 
    inner_join(enframe(smooth_enrich, 'cell', 'enrich')) %>% 
    mutate(enrich=ifelse(abs(enrich)>0.5, sign(enrich)*0.5, enrich))

p <- ggplot(plot_df, aes(C1, C2, z=enrich)) +
    stat_summary_hex(bins=50, fun='mean') +
    geom_point(data=filter(plot_df, gRNA==grna), alpha=0.5, color='black', size=1) +
    ggtitle(grna) +
    theme_void() +
    no_legend()

pb <- ggplot_build(p)
clim <- max(abs(pb$data[[1]]$value))
p + scale_fill_gradientn(colors=rev(bigrad(pals::brewer.rdbu, bias=1.5)), limits=c(-clim,clim)) 

ggsave('plots/paper/fig4_ge_density_circularhex.png', width=4, height=4)




#### Compare lineage probs ####
plot_df <- ge_meta %>% 
    filter(celltype_cl_refined%in%c('vRG', 'OPC')) %>% 
    mutate(gRNA=factor(gRNA, levels=c('Control2', 'ARID1B'))) %>% 
    filter(gRNA%in%c('ARID1B', 'Control2'))

p1 <- ggplot(plot_df, aes(gRNA, to_opc)) +
    geom_boxplot(outlier.size=0.05, fill='grey', size=0.2) +
    labs(y='OPC transition prob.')


plot_df <- ge_meta %>% 
    group_by(gRNA) %>% 
    filter(celltype_cl_refined%in%c('INP')) %>% 
    mutate(gRNA=factor(gRNA, levels=c('Control2', 'ARID1B'))) %>% 
    filter(gRNA%in%c('ARID1B', 'Control2'))

p2 <- ggplot(plot_df, aes(gRNA, to_cge)) +
    geom_boxplot(outlier.size=0.05, fill='grey', size=0.2) +
    labs(y='CGE transition prob.')

p3 <- ggplot(plot_df, aes(gRNA, to_lge2)) +
    geom_boxplot(outlier.size=0.05, fill='grey', size=0.2) +
    labs(y='LGE transition prob.')

p4 <- ggplot(plot_df, aes(gRNA, to_lge)) +
    geom_boxplot(outlier.size=0.05, fill='grey', size=0.2) +
    labs(y='OB transition prob.')

(p1 | p2 | p3 | p4) & 
    article_text() & 
    rotate_x_text(40) &
    scale_axis_rangeframe() &
    theme_rangeframe() &
    theme(
        axis.title.x = element_blank(),
        plot.margin = unit(rep(0.05, 4), 'line')
    )

ggsave('plots/paper/fig4_ge_opc_probs_boxplot.pdf', width=8, height=4.7, units='cm')



meta_test <- ge_meta %>% 
    filter(celltype_cl_refined%in%c('vRG', 'OPC')) %>% 
    filter(gRNA%in%c('Control2', 'ARID1B'))

wilcox.test(to_opc~gRNA, data=meta_test)


meta_test <- ge_meta %>% 
    filter(celltype_cl_refined%in%c('INP')) %>%  
    filter(gRNA%in%c('Control2', 'ARID1B'))

wilcox.test(to_lge2~gRNA, data=meta_test)
wilcox.test(to_lge~gRNA, data=meta_test)
wilcox.test(to_cge~gRNA, data=meta_test)



plot_df <- meta %>% filter(!is.na(terminal_states)) %>% 
    group_by(gRNA) %>% mutate(count=sum(celltype_cl_refined=='OPC')/n()) %>% 
    arrange(desc(count)) %>% mutate(gRNA=factor(gRNA, levels=unique(.$gRNA)))
ggplot(plot_df, aes(gRNA, fill=celltype_cl_refined)) +
    geom_bar(position='fill') +
    scale_fill_manual(values=asd_celltype_colors2) +
    rotate_x_text(40) +
    ggtitle('Distribution of terminal states')





