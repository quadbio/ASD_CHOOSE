source('~/scripts/single_cell/celltype.R')
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/perturbator/de.R')
source('~/scripts/perturbator/enrichment.R')
source('~/scripts/perturbator/ko_inference.R')
source('~/scripts/perturbator/plots.R')

library(Pando)

select <- dplyr::select

setwd('~/projects/ASD')



#### Read required data ####
crop <- read_rds('data/objects/CROP/asd_chong_annot_withguides_v3.3modules_srt.rds')
# crop %>% write_rds('data/objects/CROP/asd_chong_annot_withguides_v3.3modules_srt.rds')

crop_ge <- read_rds('data/objects/CROP/asd_chong_annot_wg_ge_v3.1noastro_srt.rds')

ge_cellrank_meta <- read_tsv('data/velocity/ge_opc_cellrank_probs.tsv')
ctx_cellrank_meta <- read_tsv('data/velocity/ctx_cellrank_probs.tsv')


crop %>% feature_plot(features=c('S100B', 'STMN2', 'CRABP1', 'SOX10', 'PLP1', 'TFAP2B', 'LMO4', 'ZEB2'), order=T)


#### Transfer PT ####
crop$pseudotime <- NA
crop$pseudotime[ge_cellrank_meta$CellID] <- ge_cellrank_meta$velocity_pseudotime
crop$pseudotime[ctx_cellrank_meta$CellID] <- ctx_cellrank_meta$velocity_pseudotime

crop$pseudotime_ranks <- NA
crop$pseudotime_ranks[ge_cellrank_meta$CellID] <- rank(ge_cellrank_meta$velocity_pseudotime) %>% {./max(.)}
crop$pseudotime_ranks[ctx_cellrank_meta$CellID] <- rank(ctx_cellrank_meta$velocity_pseudotime) %>% {./max(.)}


#### Foramt metadata ###
meta <- crop@meta.data %>% 
    as_tibble(rownames='cell') %>% 
    inner_join(as_tibble(crop[['umap']]@cell.embeddings, rownames='cell'))


#### UMAPS ####
ggplot(meta, aes(UMAP_1, UMAP_2, color=celltype_cl_coarse2)) +
    geom_point(size=0.2) +
    scale_color_manual(values=asd_celltype_colors2) +
    theme_void()
ggsave('plots/paper/fig1_celltype_umap.png', width=8, height=6)


ggplot(arrange(meta, !is.na(pseudotime_ranks), pseudotime_ranks), aes(UMAP_1, UMAP_2, color=pseudotime_ranks)) +
    geom_point(size=0.2) +
    scale_color_gradientn(colors=pals::magma(100)) +
    theme_void()
ggsave('plots/paper/fig1_pseudotime_umap.png', width=8, height=6)


feature_plot(crop, features=c('NEUROD6', 'DLX5', 'PAX6', 'VIM', 'APOE', 'BCL11B', 'SATB2', 'STMN2'), order=T, ncol=2)
ggsave('plots/paper/fig1_feature_umap.png', width=8, height=16)

feature_plot(crop, features=c('MEIS2', 'DLX5', 'DLX2', 'TBR1', 'FOXG1'), order=T, ncol=2)



#### Ventral subclustering ####
ge_meta <- crop_ge@meta.data %>% 
    as_tibble(rownames='cell') %>% 
    inner_join(as_tibble(crop_ge[['umap']]@cell.embeddings, rownames='cell'))

ggplot(ge_meta, aes(UMAP_1, UMAP_2, color=celltype_cl_refined)) +
    geom_point(size=0.8) +
    scale_color_manual(values=asd_celltype_colors2) +
    theme_void()
ggsave('plots/paper/fig1_ge_celltype_umap.png', width=8, height=3.7)


feature_plot(crop_ge, features=c('OLIG2', 'PDGFRA', 'SOX10', 'PCDH15', 'OLIG1', 'LHFPL3'), order=T, ncol=2, pt.size=0.8)
ggsave('plots/paper/fig1_ge_feature_umap.png', width=6, height=6)


#### Marker heatmap ####
crop <- aggregate_assay(crop, group_name = 'celltype_cl_coarse2', assay = 'RNA')

markers_plot <- c('VIM', 'PAX6', 'ASPM', 'GFAP', 'ASCL1',
                  'MKI67', 'TOP2A',
                  'HOPX', 'EOMES', 'TAC3',
                  'BCL11B', 'SATB2', 'NEUROD6',
                  'RORB', 'UNC5D', 'NR2F1',
                  'TLE4', 'FOXP2',
                  'OLIG1', 'OLIG2', 'PDGFRA',
                  'ADARB2', 'CALB2',
                  'DLX2', 'DLX5', 'NR2F2', 'MEIS2', 'GAD1',
                  'S100B', 'APOE', 'ALDH1L1',
                  'AUTS2') %>% 
    intersect(colnames(crop@assays$RNA@misc$summary$celltype_cl_coarse))

# markers_plot <- c('VIM','PAX6','MKI67','HOPX','PTPRZ1','EOMES','TAC3','BCL11B','BCL11A','CTIP2',
#     'STMN2','SATB2','NR2F1','TLE4','NEUROD6','TBR1','OLIG1', 'LIX1','SCGN','GAD2','DLX2','DLX5',
#     'COUPTF2','ADARB2','CALB2','S100B','APOE') %>% 
#     intersect(colnames(crop@assays$RNA@misc$summary$celltype_cl_coarse))

expr_mat <- crop@assays$RNA@misc$summary$celltype_cl_coarse2[, markers_plot]
clust_order <- expr_mat %>% dist() %>% hclust(method='ward.D2') %>% {.$label[.$order]}
gene_order <- expr_mat %>% t() %>% dist() %>% hclust() %>% {.$label[.$order]}

ct_order <- names(asd_celltype_colors2)

expr_df <- expr_mat %>% 
    as_tibble(rownames='celltype') %>% 
    pivot_longer(!celltype) %>% 
    filter(celltype!='mesenchyme') %>% 
    mutate(celltype=factor(celltype, levels=ct_order), name=factor(name, levels=rev(markers_plot))) %>% 
    group_by(name) %>% mutate(value=scale01(value))

ggplot(expr_df, aes(celltype,name, fill=value)) +
    geom_tile() +
    scale_fill_gradientn(colors=pals::brewer.greys(100)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    rotate_x_text(40) +
    no_legend() +
    article_text() +
    theme(
        axis.ticks = element_line(size=0.3),
        panel.border = element_rect(size=0.2),
        plot.margin = unit(rep(0.05,4), 'line')
    ) +
    labs(x='Cell type', y='Gene')

ggsave('plots/paper/fig1_feature_heatmap.png', width=6.4, height=9.5, bg='white', unit='cm')
ggsave('plots/paper/fig1_feature_heatmap.pdf', width=6.4, height=9.5, bg='white', unit='cm')



#### Supp plots (QC, composition) ####
p1 <- ggplot(meta, aes(nCount_RNA)) +
    geom_histogram(color='black', fill='grey', size=0.2) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    labs(x='UMI count', y='# cells') +
    theme(
        plot.margin = unit(rep(0.2,4), 'lines')
    )

p2 <- ggplot(meta, aes(nFeature_RNA)) +
    geom_histogram(color='black', fill='grey', size=0.2) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    labs(x='# features', y='# cells') +
    theme(
        plot.margin = unit(rep(0.2,4), 'lines')
    )

p1 / p2
ggsave('plots/paper/fig1_qc_hist.pdf', width=5, height=3.5, unit='cm')


ggplot(meta, aes(library, fill=factor(celltype_cl_coarse2, levels=names(asd_celltype_colors2)))) +
    geom_bar(size=0.2, position='fill') +
    scale_fill_manual(values=asd_celltype_colors2) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    no_x_text() +
    labs(x='Library', y='Fraction of cells') +
    theme(
        plot.margin = unit(rep(0.2,4), 'lines'),
        axis.line.x = element_blank()
    ) +
    no_legend()

ggsave('plots/paper/fig1_ct_copm_bar.pdf', width=3, height=3.5, unit='cm')



ggplot(meta, aes(gRNA, fill=library)) +
    geom_bar(size=0.1, color='black') +
    scale_axis_rangeframe() + theme_rangeframe() +
    scale_fill_manual(values=pals::brewer.greys(8)) +
    article_text() +
    rotate_x_text(90) +
    scale_y_continuous(expand=c(0,0)) +
    labs(x='gRNA', y='# cells') +
    theme(
        plot.margin = unit(rep(0.2,4), 'lines'),
        axis.line.x = element_blank()
    ) +
    no_legend()

ggsave('plots/paper/fig1_grna_comp_bar.pdf', width=7, height=3.5, unit='cm')


grna_names <- sort(unique(meta$gRNA))
grna_names <- grna_names[grna_names!='Control2']
grna_colors <- many[1:length(grna_names)]
names(grna_colors) <- grna_names
grna_colors['Control2'] <- 'white'


ggplot(meta, aes(library, fill=gRNA)) +
    geom_bar(size=0.1, color='black') +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_manual(values=grna_colors) +
    labs(x='Library', y='# cells') +
    theme(
        plot.margin = unit(rep(0.2,4), 'lines'),
        axis.line.x = element_blank()
    ) +
    no_legend() +
    no_x_text()

ggsave('plots/paper/fig1_grna_library_comp_bar.pdf', width=4, height=3.5, unit='cm')


markers_plot <- c('VIM', 'PAX6', 'MKI67', 'EOMES', 'TAC3','BCL11B', 
                  'SATB2', 'NR2F1','OLIG1', 'PDGFRA','DLX2', 'DLX5', 
                  'NR2F2', 'MEIS2','S100B', 'APOE', 'ALDH1L1','AUTS2')

feature_plot(crop, features=markers_plot, order=T, ncol=6, pt.size=0.001) & 
    article_text() &
    theme_void() & no_legend() &
    theme(plot.title=element_text(face='plain', size=6))
ggsave('plots/paper/fig1_feature_plots_umap.png', width=18, height=8, unit='cm')




