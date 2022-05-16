source('~/scripts/single_cell/celltype.R')
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/single_cell/cell_density_diff.R')
source('~/scripts/perturbator/de.R')
source('~/scripts/perturbator/enrichment.R')
source('~/scripts/perturbator/ko_inference.R')
source('~/scripts/perturbator/plots.R')

library(Pando)
library(clusterProfiler)
library(org.Hs.eg.db)
library(doParallel)
library(eulerr)
library(ggraph)
library(tidygraph)
library(ggpointdensity)

select <- dplyr::select
filter <- dplyr::filter

setwd('~/projects/ASD')


#### Read required data ####
crop <- read_rds('data/objects/CROP/asd_chong_annot_withguides_v3.3modules_srt.rds')
crop_muo_wt <- read_rds('data/objects/CROP_MUO/ASD_CROP_MUO_wt_pando_v2net_srt.rds')

sfari_genes <- readr::read_csv('data/SFARI-Gene_genes_09-02-2021release_11-04-2021export.csv')
sysid_genes <- readr::read_csv('data/SysID_17.03.2022.csv') %>% 
    {colnames(.) <- str_replace(colnames(.), ' ', '_');.}
tfs <- read_tsv('~/resources/DB/animal_tfdb/tf_human.tsv')

perc_expr <- rowMeans(crop[['RNA']]@data > 0)
expr_genes <- names(perc_expr[perc_expr>0.05])

expr_genes %>% write('data/results/diff_expression/DE_bg_genes.txt')

all_features <- bitr(expr_genes, fromType = 'SYMBOL', toType = c('ENSEMBL', 'ENTREZID'), OrgDb = org.Hs.eg.db) %>%
    as_tibble()

de_lin_lr <- read_tsv('data/results/diff_expression/lineage_lr_cl_refined_de.tsv')

deg_label_ex <- readxl::read_xlsx('data/DEG_sig_selection_Plot Fig3a_23.03.2022 shared.xlsx', sheet = 1) %>% 
    select(feature, gRNA, lineage, Function) %>% mutate(label=T)

deg_label_in <- readxl::read_xlsx('data/DEG_sig_selection_Plot Fig3a_23.03.2022 shared.xlsx', sheet = 2) %>% 
    select(feature, gRNA, lineage, Function) %>% mutate(label=T)

deg_label <- bind_rows(deg_label_ex, deg_label_in)


#### Lineage DE ####
de_lin <- de_lin_lr %>% 
    mutate(p_val_adj=p.adjust(p_val, method='fdr')) %>% 
    mutate(spval = sign(avg_log2FC)*-log10(p_val)) %>% 
    mutate(spval_clip = ifelse(abs(spval)>10, sign(spval)*10, spval)) %>% 
    filter(!str_detect(feature, 'masked|Cas9|SEPTIN6')) 

de_sig <- de_lin %>% 
    # filter(abs(avg_log2FC)>0.3) %>%
    filter(p_val_adj<0.05) %>%
    group_by(gRNA, lineage) 

de_sig_top <- de_sig %>% 
    arrange(desc(abs(avg_log2FC))) %>% 
    filter(feature%in%all_features$SYMBOL) %>% 
    filter(row_number()<30)

de_sig_top %>% write_tsv('data/results/diff_expression/lineage_top_DEG_lr.tsv')
de_sig %>% write_tsv('data/results/diff_expression/lineage_DEG_lr.tsv')


de_lin_sig <- de_lin %>% 
    filter(p_val_adj<0.4)

ggplot(de_lin_sig, aes(gRNA, -log10(p_val_adj), color=p_val_adj<0.05)) +
    geom_jitter(size=0.01, width=0.2) +
    geom_hline(yintercept = -log10(0.05), color='darkgrey', size=0.3, linetype='dashed') +
    scale_color_manual(values=c('grey', orange)) +
    facet_grid(lineage~.) +
    article_text() + rotate_x_text(40) + no_legend() +
    scale_axis_rangeframe() + theme_rangeframe() +
    labs(x='Target gene', y='-log10(FDR)')

ggsave('plots/paper/fig3_de_ex_in_manhatten.png', width=9, height=5, units = 'cm', bg='white')
ggsave('plots/paper/fig3_de_ex_in_manhatten.pdf', width=9, height=5, units = 'cm')



plot_df <- de_lin %>% 
    filter(p_val_adj<0.4) %>% 
    left_join(deg_label) 

ggplot(plot_df, aes(gRNA, -log10(p_val_adj), color=p_val_adj<0.05, label=feature)) +
    geom_jitter(size=0.01, width=0.2) +
    geom_hline(yintercept = -log10(0.05), color='darkgrey', size=0.3, linetype='dashed') +
    geom_text_repel(
        data=filter(plot_df, !is.na(Function)), 
        nudge_y=1000,
        angle=90,
        size=5/ggplot2::.pt,
        color='black',
        segment.size=0.2,
        force=2,
        min.segment.length=0
    ) +
    scale_color_manual(values=c('grey', orange)) +
    scale_y_continuous(limits=c(0,20), na.value = 20) +
    facet_grid(lineage~.) +
    article_text() + rotate_x_text(40) + no_legend() +
    scale_axis_rangeframe() + theme_rangeframe() +
    theme(strip.text = element_blank()) +
    labs(x='Target gene', y='-log10(FDR)')

ggsave('plots/paper/fig3_de_ex_in_manhatten_labelled.png', width=10, height=6, units = 'cm', bg='white')
ggsave('plots/paper/fig3_de_ex_in_manhatten_labelled.pdf', width=10, height=6, units = 'cm')


plot_df <- de_sig %>% 
    group_by(feature, lineage) %>% 
    mutate(
        mean_fc = mean(avg_log2FC),
        std_fc = sd(avg_log2FC),
        min_fc = min(avg_log2FC),
        max_fc = max(avg_log2FC),
        n_kos = n()
    )

dodge_width <- 0.8
ggplot(plot_df, aes(n_kos, avg_log2FC, label=feature, color=lineage, group=feature)) +
    geom_hline(yintercept = 0, color='darkgrey', size=0.5) +
    geom_linerange(aes(ymin=min_fc, ymax=max_fc), position=position_dodge(dodge_width), color='grey', alpha=0.5, size=0.2) +
    geom_point(position=position_dodge(dodge_width), size=0.15) +
    # geom_point(data=filter(plot_df, n_kos>1), aes(y=mean_fc), position=position_dodge(dodge_width), size=1.3) +
    geom_text_repel(
        data=filter(plot_df, n_kos>3, avg_log2FC==max_fc),
        size=5/ggplot2::.pt,
        color='black',
        segment.size=0.2,
        force=2,
        min.segment.length=0,
        position=position_dodge(dodge_width)
    ) +
    article_text() +
    scale_axis_rangeframe() + theme_rangeframe() +
    scale_color_manual(values=asd_lineage_colors) +
    no_legend() +
    no_x_text() +
    theme(
        plot.margin = unit(rep(0,4), 'lines'),
        axis.line.x = element_blank(),
        axis.title.x = element_blank()
    ) +
    labs(y='log2(fold change)')

ggsave('plots/paper/fig3_de_nkos_jitter.png', width=8.3, height=3.5, units = 'cm', bg='white')
ggsave('plots/paper/fig3_de_nkos_jitter.pdf', width=8.3, height=3.5, units = 'cm')



#### Do pathway enrichment ####

de_genes_all <- de_sig %>%
    group_by(lineage) %>%
    group_split() %>%
    {names(.) <- map_chr(., ~.x$lineage[1]); .}

go_results_all <- map(de_genes_all, function(x){
    gene_ids <- filter(all_features, SYMBOL%in%x$feature)
    ego <- enrichGO(
        gene = gene_ids$ENTREZID,
        universe = all_features$ENTREZID,
        OrgDb = org.Hs.eg.db,
        pvalueCutoff = 0.05,
        # qvalueCutoff = 0.05,
        pAdjustMethod = 'fdr',
        ont = 'ALL',
        readable = T
    )
    print(ego)
    return(ego)
})

go_all_df <- map_dfr(go_results_all, ~if (!is.null(.x)){as_tibble(.x@result)} , .id='KO')

go_all_df <- go_all_df %>%
    mutate(ngenes = as.numeric(str_replace(GeneRatio, '\\d+/(\\d+)', '\\1'))) %>%
    mutate(ngenes_enrich = as.numeric(str_replace(GeneRatio, '(\\d+)/\\d+', '\\1'))) %>%
    filter(ngenes>5, ngenes_enrich>1)
go_all_df %>% write_tsv('data/results/diff_expression/lineage_de_all_GO_enrich.tsv')


kegg_results_all <- map(de_genes_all, function(x){
    gene_ids <- filter(all_features, SYMBOL%in%x$feature)
    ego <- enrichKEGG(
        gene = gene_ids$ENTREZID,
        universe = all_features$ENTREZID,
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        pAdjustMethod = 'fdr'
    )
    print(ego)
    return(ego)
})

kegg_all_df <- map_dfr(kegg_results_all, ~if (!is.null(.x)){as_tibble(.x@result)} , .id='KO')

kegg_all_df <- kegg_all_df %>% 
    mutate(ngenes = as.numeric(str_replace(GeneRatio, '\\d+/(\\d+)', '\\1'))) %>% 
    mutate(ngenes_enrich = as.numeric(str_replace(GeneRatio, '(\\d+)/\\d+', '\\1'))) %>% 
    filter(ngenes>5, ngenes_enrich>1) 

kegg_all_df %>% write_tsv('data/results/diff_expression/lineage_de_all_KEGG_enrich.tsv')




#### DEG enrichment in SFARI/sysID ####
all_sig_genes <- de_sig$feature %>% unique()

# enrichtest_df <- map_dfr(set_names(unique(de_sig$gRNA)), function(x){
#     print(x)
#     
#     all_use <- de_lin %>% 
#         filter(gRNA==x, feature%in%expr_genes) %>% 
#         arrange(p_val_adj) %>% pull(feature) %>% unique()
#     
#     sig_use <- de_sig %>% 
#         filter(gRNA==x, feature%in%expr_genes) %>% 
#         arrange(p_val_adj) %>% pull(feature) %>% unique()
#     
#     kstest <- ks.test(which(all_use%in%sfari_genes$`gene-symbol`), 1:length(all_use), alternative='less')
#     
#     test_table <- rbind(
#         table(expr_genes%in%sfari_genes$`gene-symbol`),
#         table(sig_use%in%sfari_genes$`gene-symbol`)
#     )
#     
#     ftest <- fisher.test(test_table)
#     
#     return(tibble(
#         kstest = kstest$p.value,
#         fishtest = ftest$p.value,
#         logodds = log(ftest$estimate),
#         ngenes = length(sig_use),
#         total_ratio = test_table[2,1]/test_table[1,1],
#         in_ratio = test_table[2,2]/test_table[1,2]
#     ))
#     
# }, .id='target')


test_table <- rbind(
    table(expr_genes%in%sfari_genes$`gene-symbol`),
    table(all_sig_genes%in%sfari_genes$`gene-symbol`)
)

ftest_sfari <- fisher.test(test_table)

primid_genes <- sysid_genes %>% 
    filter(Gene_group=='Current primary ID genes')

test_table <- rbind(
    table(expr_genes%in%primid_genes$Gene_symbol),
    table(all_sig_genes%in%primid_genes$Gene_symbol)
)

ftest_sysid <- fisher.test(test_table)

ftest_df <- tibble(
    pval=c(ftest_sfari$p.value,ftest_sysid$p.value),
    odds_ratio=c(ftest_sfari$estimate,ftest_sysid$estimate),
    origin=c('SFARI', 'sysID'),
)
    
p1 <- ggplot(ftest_df, aes(-log10(pval), origin)) +
    geom_bar(stat='identity') +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text()

p2 <- ggplot(ftest_df, aes(log2(odds_ratio), origin)) +
    geom_bar(stat='identity') +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text()

p1 / p2


ggplot(ftest_df, aes(log2(odds_ratio), origin)) +
    geom_bar(stat='identity', color='black', fill='grey', size=0.2) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    theme(
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(rep(0,4), 'lines'),
    ) +
    labs(x='log2(fold enrichment)')

ggsave('plots/paper/fig3_deg_enrichment_sysid_sfari_bar.png', width=4, height=1.5, units='cm', bg='white')
ggsave('plots/paper/fig3_deg_enrichment_sysid_sfari_bar.pdf', width=4, height=1.5, units='cm')


intersect(sfari_genes$`gene-symbol`, primid_genes$Gene_symbol) %>% length()

ex_sig_genes <- de_sig %>% filter(lineage=='excitatory') %>% pull(feature) %>% unique()
in_sig_genes <- de_sig %>% filter(lineage=='inhibitory') %>% pull(feature) %>% unique()

sfari_intersect <- intersect(all_sig_genes, sfari_genes$`gene-symbol`)

euler_plot <- euler(c(
    'a' = length(ex_sig_genes),
    'b' = length(in_sig_genes),
    'c' = length(sfari_genes$`gene-symbol`),
    'a&b' = length(intersect(ex_sig_genes, in_sig_genes)),
    'b&c' = length(intersect(in_sig_genes, sfari_genes$`gene-symbol`)),
    'a&c' = length(intersect(ex_sig_genes, sfari_genes$`gene-symbol`)),
    'a&b&c' = length(intersect(intersect(ex_sig_genes, sfari_genes$`gene-symbol`), in_sig_genes))
))

pdf('plots/paper/fig3_SFARI_DEG_venn.pdf')
plot(euler_plot)
dev.off()


kriegstein_deg <- readxl::read_xls('data/DEG_Kriegstein_Science.xls') %>% 
    {colnames(.) <- str_replace(colnames(.), ' ', '_');.}

krig_excit_deg <- kriegstein_deg %>% 
    filter(q_value<0.05, Cell_type%in%c('L2/3', 'L4', 'L5/6', 'L5/6-CC', 'Neu-mat')) %>% 
    pull(Gene_name) %>% unique()

euler_plot <- euler(c(
    'a' = length(ex_sig_genes),
    'b' = length(krig_excit_deg),
    'a&b' = length(intersect(ex_sig_genes, krig_excit_deg))
))

pdf('plots/paper/fig3_KS_ex_DEG_venn.pdf')
plot(euler_plot)
dev.off()



all_genes <- union(ex_sig_genes, in_sig_genes) %>% union(sfari_genes$`gene-symbol`) %>% union(primid_genes$Gene_symbol)
euler_df <- data.frame(
    'excitatory'=all_genes%in%ex_sig_genes,
    'inhibitory'=all_genes%in%in_sig_genes,
    'SAFRI'=all_genes%in%sfari_genes$`gene-symbol`,
    'sysID'=all_genes%in%primid_genes$Gene_symbol
)

euler_plot <- venn(euler_df)

pdf('plots/paper/fig3_SFARI_DEG_venn.pdf')
plot(euler_plot)
dev.off()


#### GO enrichment of DEG ####
gene_ids <- filter(all_features, SYMBOL%in%sfari_intersect)
ego <- enrichGO(
    gene = gene_ids$ENTREZID,
    universe = all_features$ENTREZID,
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    pAdjustMethod = 'fdr',
    ont = 'ALL',
    readable = T
)

gene_ids <- filter(all_features, SYMBOL%in%sfari_genes$`gene-symbol`)
ego_sfari <- enrichGO(
    gene = gene_ids$ENTREZID,
    universe = all_features$ENTREZID,
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    pAdjustMethod = 'fdr',
    ont = 'ALL',
    readable = T
)

gene_ids <- filter(all_features, SYMBOL%in%all_sig_genes)
ego_choose <- enrichGO(
    gene = gene_ids$ENTREZID,
    universe = all_features$ENTREZID,
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    pAdjustMethod = 'fdr',
    ont = 'ALL',
    readable = T
)

as_tibble(ego_choose@result) %>% write_tsv('data/results/diff_expression/all_CHOOSE_DE_GO_enrich.tsv')

gene_ids <- filter(all_features, SYMBOL%in%in_sig_genes)
ego_in <- enrichGO(
    gene = gene_ids$ENTREZID,
    universe = all_features$ENTREZID,
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    pAdjustMethod = 'fdr',
    ont = 'ALL',
    readable = T
)

as_tibble(ego_in@result) %>% write_tsv('data/results/diff_expression/all_CHOOSE_inhib_DE_GO_enrich.tsv')

gene_ids <- filter(all_features, SYMBOL%in%ex_sig_genes)
ego_ex <- enrichGO(
    gene = gene_ids$ENTREZID,
    universe = all_features$ENTREZID,
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    pAdjustMethod = 'fdr',
    ont = 'ALL',
    readable = T
)

as_tibble(ego_ex@result) %>% write_tsv('data/results/diff_expression/all_CHOOSE_excit_DE_GO_enrich.tsv')


gene_ids <- filter(all_features, SYMBOL%in%intersect(intersect(ex_sig_genes, sfari_genes$`gene-symbol`), in_sig_genes))
ego_intersect <- enrichGO(
    gene = gene_ids$ENTREZID,
    universe = all_features$ENTREZID,
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    pAdjustMethod = 'fdr',
    ont = 'ALL',
    readable = T
)

as_tibble(ego_intersect@result) %>% write_tsv('data/results/diff_expression/all_CHOOSE_intersect_DE_GO_enrich.tsv')


#### Plot enrichment ####
ego_choose <- read_tsv('data/results/diff_expression/all_CHOOSE_DE_GO_enrich.tsv')

plot_df <- ego_choose %>% 
    filter(ONTOLOGY=='BP', p.adjust<0.01) %>% 
    mutate(ngenes = as.numeric(str_replace(GeneRatio, '\\d+/(\\d+)', '\\1'))) %>% 
    mutate(ngenes_enrich = as.numeric(str_replace(GeneRatio, '(\\d+)/\\d+', '\\1'))) %>% 
    mutate(bggenes = as.numeric(str_replace(BgRatio, '\\d+/(\\d+)', '\\1'))) %>% 
    mutate(bggenes_enrich = as.numeric(str_replace(BgRatio, '(\\d+)/\\d+', '\\1'))) %>% 
    filter(ngenes_enrich>20) %>% 
    mutate(oddsratio=(ngenes_enrich/ngenes)/(bggenes_enrich/bggenes)) %>% 
    arrange(-p.adjust) %>% 
    mutate(Description=factor(Description, levels=unique(.$Description))) %>% 
    arrange(p.adjust) %>% 
    filter(row_number()<=20)
    

ggplot(plot_df, aes(-log10(p.adjust), Description, label=Description, fill=oddsratio)) +
    geom_bar(stat='identity', color='black', size=0.2) +
    geom_text(hjust=0, size=5/ggplot2::.pt) +
    scale_fill_gradientn(colors=pals::brewer.greys(30), limits=c(1,4), na.value='red') +
    theme_rangeframe() + scale_axis_rangeframe() +
    article_text() +
    no_y_text() + no_legend()

ggsave('plots/paper/fig3_DEG_GO_enrich_bar.pdf', width=3.8, height=5, units='cm')



#### Plot multiome/GRN parmeters ####
muo_meta <- crop_muo_wt@meta.data %>% 
    as_tibble() %>% 
    filter(dataset=='MUO')


p1 <- ggplot(muo_meta, aes(nCount_RNA)) +
    geom_histogram(color='black', fill='grey', size=0.2) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    labs(x='UMI count', y='# cells') +
    theme(
        plot.margin = unit(rep(0.2,4), 'lines')
    )

p2 <- ggplot(muo_meta, aes(nFeature_RNA)) +
    geom_histogram(color='black', fill='grey', size=0.2) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    labs(x='# features', y='# cells') +
    theme(
        plot.margin = unit(rep(0.2,4), 'lines')
    )

p3 <- ggplot(muo_meta, aes(nCount_ATAC)) +
    geom_histogram(color='black', fill='grey', size=0.2) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    labs(x='Fragment count', y='# cells') +
    theme(
        plot.margin = unit(rep(0.2,4), 'lines')
    )

p4 <- ggplot(muo_meta, aes(nFeature_ATAC)) +
    geom_histogram(color='black', fill='grey', size=0.2) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    labs(x='# peaks', y='# cells') +
    theme(
        plot.margin = unit(rep(0.2,4), 'lines')
    )

(p1 / p2) | (p3 / p4)

ggsave('plots/paper/fig3_muo_qc_hist.pdf', width=10, height=3.5, unit='cm')


muo <- read_rds('data/objects/MUO/H9_joined_srt.rds')

DefaultAssay(muo) <- 'RNA'
p1 <- dim_plot(muo, label=T)
p2 <- feature_plot(muo, features=c('NEUROD6', 'DLX5', 'SATB2', 'GAD2'), order=T)

muo$lineage <- case_when(
    muo$seurat_clusters%in%c(5,9,13) ~ 'inhibitory',
    T ~ 'excitatory',
)

p1 <- CoveragePlot(
    object = muo,
    assay = 'ATAC',
    group.by = 'lineage',
    region = 'NEUROD6', 
    extend.downstream = 2000
) & article_text()
p1
ggsave('plots/paper/fig3_muo_NEUROD6_tracks.pdf', width=10, height=10, unit='cm')

p2 <- CoveragePlot(
    object = muo,
    assay = 'ATAC',
    group.by = 'lineage',
    region = 'GAD2', 
    extend.downstream = 2000,
    extend.upstream = 2000
) & article_text()
p2
ggsave('plots/paper/fig3_muo_GAD2_tracks.pdf', width=10, height=10, unit='cm')




crop_muo_all <- read_rds('data/objects/CROP_MUO/ASD_CROP_MUO_all_srt.rds')
dim_plot(crop_muo_all, group.by='dataset', order=T, pt.size=1)


crop_muo_all <- AddMetaData(crop_muo_all, crop@meta.data['celltype_cl_coarse2'])
dim_plot(crop_muo_all, group.by='celltype_cl_coarse2', order=T, pt.size=1)

crop_muo_meta <- crop_muo_all@meta.data %>% 
    as_tibble(rownames='cell') %>% 
    inner_join(as_tibble(crop_muo_all[['umap']]@cell.embeddings, rownames='cell')) %>% 
    mutate(
        alpha = (gRNA == 'Control2') | (dataset=='MUO'),
        color = ifelse(dataset=='MUO', 'MUO', 'CTRL')
    ) %>% 
    filter(UMAP_1<10)
    
ggplot(crop_muo_meta, aes(UMAP_1, UMAP_2, color=celltype_cl_coarse2)) +
    geom_point(data=filter(crop_muo_meta, !alpha), size=1, alpha=0.2) +
    geom_point(data=filter(crop_muo_meta, alpha & color=='CTRL'), size=2, shape=21, fill='#BBBBBB', color='black', stroke=0.5) +
    geom_point(data=filter(crop_muo_meta, alpha & color=='MUO'), size=2, shape=21, fill='#444444', color='black', stroke=0.5) +
    scale_color_manual(values=asd_celltype_colors2) +
    theme_void() +
    no_legend()
    

ggsave('plots/paper/fig3_muo_cell_dist_umap.png', width=15, height=15, unit='cm')


feature_plot(crop_muo_all, features=c('SATB2', 'BCL11B', 'VIM', 'PAX6', 'DLX2', 'GAD2', 'OLIG2', 'NEUROD6', 'DLX5'), order=T, pt.size=0.001, ncol=3) & 
    article_text() &
    theme_void() & no_legend() &
    theme(plot.title=element_text(face='plain', size=6))
ggsave('plots/paper/fig3_muo_feature_umap.png', width=13, height=13, unit='cm')





#### Get GRN ####
grn_features <- NetworkFeatures(crop_muo_wt, network='glm_bidir_GREAT')
grn_module_features <- NetworkModules(crop_muo_wt, network='glm_bidir_GREAT')@features
grn_modules_pos <- grn_module_features$genes_pos
grn_modules_neg <- grn_module_features$genes_neg
grn_modules_full <- list_merge(grn_modules_pos, !!!grn_modules_neg)
grn_modules <- NetworkModules(crop_muo_wt, network='glm_bidir_GREAT')@meta %>% 
    filter(n_genes>=5)


#### General GRN metrics ####
grn_coefs <- coef(crop_muo_wt, network='glm_bidir_GREAT')
grn_gof <- gof(crop_muo_wt, network='glm_bidir_GREAT') %>% 
    dplyr::filter(rsq<=1, rsq>=0) %>% 
    mutate(nice=rsq>0.10&nvariables>10) 

## GOF metrics
p1 <- ggplot(grn_gof, aes(rsq, nvariables, alpha=nice)) +
    geom_pointdensity(size=0.01, shape=16) +
    geom_hline(yintercept=10, size=0.2, color='darkgrey', linetype='dashed') +
    geom_vline(xintercept=0.1, size=0.2, color='darkgrey', linetype='dashed') +
    scale_color_gradientn(colors=rev(reds(1.3))) +
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(0, 1, 10, 100, 1000, 10000)) +
    scale_x_continuous(breaks=seq(0,1,0.2)) +
    article_text() +
    no_legend() +
    labs(x=expression('Explained variance'~(R**2)), y='# variables in model') +
    theme(
        plot.margin = unit(c(0,0,0,0), 'line'),
        strip.text = element_blank()
    )

p2 <- ggplot(grn_gof, aes(rsq)) +
    geom_histogram(fill='darkgray', bins=20, color='black', size=0.2) +
    theme_void() +
    no_legend()

p3 <- ggplot(grn_gof, aes(nvariables)) +
    geom_histogram(fill='darkgray', bins=20, color='black', size=0.2) +
    scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(0, 1, 10, 100, 1000, 10000)) +
    theme_void() +
    coord_flip() +
    no_legend() 


layout <- '
AAAA#
BBBBC
BBBBC
'
p2 + p1 + p3 + plot_layout(design = layout) & no_margin()
ggsave('plots/paper/fig3_nvar_vs_r2_scatter.pdf', width=4.1, height=3, units='cm')


## Module sizes
plot_df <- grn_modules %>% 
    distinct(target, n_regions)

p1 <- ggplot(plot_df, aes(1, n_regions)) +
    theme_rangeframe() + scale_axis_rangeframe() +
    article_text() +
    no_x_text() +
    theme(
        axis.line.x = element_blank(),
        axis.title.x = element_blank()
    ) +
    scale_x_discrete(expand=c(0.1,0.1)) +
    geom_boxplot(width=0.2, size=0.2, fill='darkgrey', outlier.size = 0.1) +
    labs(y='# peaks')

plot_df <- grn_modules %>% 
    distinct(target, n_tfs)

p2 <- ggplot(plot_df, aes(1, n_tfs)) +
    theme_rangeframe() + scale_axis_rangeframe() +
    article_text() +
    no_x_text() +
    scale_x_discrete(expand=c(0.1,0.1)) +
    theme(
        axis.line.x = element_blank(),
        axis.title.x = element_blank()
    ) +
    geom_boxplot(width=0.2, size=0.2, fill='darkgrey', outlier.size = 0.1) +
    labs(y='# TFs')


plot_df <- grn_modules %>% 
    distinct(tf, n_genes)

p3 <- ggplot(plot_df, aes(1, n_genes)) +
    theme_rangeframe() + scale_axis_rangeframe() +
    article_text() +
    no_x_text() +
    scale_x_discrete(expand=c(0.1,0.1)) +
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(0, 5, 50, 500)) +
    theme(
        axis.line.x = element_blank(),
        axis.title.x = element_blank()
    ) +
    geom_boxplot(width=0.2, size=0.2, fill='darkgrey', outlier.size = 0.1) +
    labs(y=expression('# genes'))

p1 | p2 | p3 & no_margin()
ggsave('plots/paper/fig3_model_module_specs_violin.pdf', width=5, height=4.5, units='cm')





#### Test each module for enrichment ####
module_sfari_enrich <- map_dfr(grn_modules_full, function(x){
    bg_tab <- table(grn_features%in%sfari_genes$`gene-symbol`)
    x_tab <- table(x%in%sfari_genes$`gene-symbol`)
    test_table <- rbind(
        c(bg_tab['FALSE'], bg_tab['TRUE']),
        c(x_tab['FALSE'], x_tab['TRUE'])
    )
    test_table[is.na(test_table)] <- 0
    ftest <- fisher.test(test_table)
    return(tibble(
        pval = ftest$p.value,
        logodds = log2(ftest$estimate),
        ngenes = length(x),
        ngenes_sfari = test_table[[2,'TRUE']]
    ))    
}, .id='tf') %>% mutate(padj=p.adjust(pval, method='fdr'))

pos_module_sfari_enrich <- map_dfr(grn_modules_pos, function(x){
    bg_tab <- table(grn_features%in%sfari_genes$`gene-symbol`)
    x_tab <- table(x%in%sfari_genes$`gene-symbol`)
    test_table <- rbind(
        c(bg_tab['FALSE'], bg_tab['TRUE']),
        c(x_tab['FALSE'], x_tab['TRUE'])
    )
    test_table[is.na(test_table)] <- 0
    ftest <- fisher.test(test_table)
    return(tibble(
        pval = ftest$p.value,
        logodds = log2(ftest$estimate),
        ngenes = length(x),
        ngenes_sfari = test_table[[2,'TRUE']]
    ))    
}, .id='tf') %>% mutate(padj=p.adjust(pval, method='fdr'))

neg_module_sfari_enrich <- map_dfr(grn_modules_neg, function(x){
    bg_tab <- table(grn_features%in%sfari_genes$`gene-symbol`)
    x_tab <- table(x%in%sfari_genes$`gene-symbol`)
    test_table <- rbind(
        c(bg_tab['FALSE'], bg_tab['TRUE']),
        c(x_tab['FALSE'], x_tab['TRUE'])
    )
    test_table[is.na(test_table)] <- 0
    ftest <- fisher.test(test_table)
    return(tibble(
        pval = ftest$p.value,
        logodds = log2(ftest$estimate),
        ngenes = length(x),
        ngenes_sfari = test_table[[2,'TRUE']]
    ))    
}, .id='tf') %>% mutate(padj=p.adjust(pval, method='fdr'))


both_df <- module_sfari_enrich %>% 
    filter(ngenes_sfari>0) %>% 
    arrange(desc(sign(logodds)*-log10(padj))) %>% 
    mutate(tf=factor(tf, levels=unique(.$tf)))

neg_df <- neg_module_sfari_enrich %>% 
    filter(ngenes_sfari>0) %>% 
    arrange(desc(sign(logodds)*-log10(padj))) %>% 
    mutate(tf=factor(tf, levels=unique(.$tf)))

pos_df <- pos_module_sfari_enrich %>% 
    filter(ngenes_sfari>0) %>% 
    arrange(desc(sign(logodds)*-log10(padj))) %>% 
    mutate(tf=factor(tf, levels=unique(.$tf)))

p1 <- ggplot(both_df, aes(x=tf, y=sign(logodds)*-log10(padj), fill=padj<0.01)) +
    geom_bar(stat='identity') +  
    scale_fill_manual(values=c('darkgrey', '#AB6472')) +
    rotate_x_text(40) +
    ggtitle('combined')

p2 <- ggplot(neg_df, aes(x=tf, y=sign(logodds)*-log10(padj), fill=padj<0.01)) +
    geom_bar(stat='identity') +  
    scale_fill_manual(values=c('darkgrey', '#AB6472')) +
    rotate_x_text(40) +
    ggtitle('negative')

p3 <- ggplot(pos_df, aes(x=tf, y=sign(logodds)*-log10(padj), fill=padj<0.01)) +
    geom_bar(stat='identity') +  
    scale_fill_manual(values=c('darkgrey', '#AB6472')) +
    rotate_x_text(40) +
    ggtitle('positive')

p1 / p2 / p3 

ggsave('plots/MUO/grn_sfari_score1_enrich_logodds_bar.png', width=25, height=15, bg='white')



tf_module_enrich <- bind_rows('combined'=both_df, 'negative'=neg_df, 'positive'=pos_df, .id='module_type')
tf_module_enrich %>% write_tsv('data/results/grn/grn_sfari_module_enrichment.tsv')

#### Plots for TF module / SFARI enrichment ####
tf_module_enrich <- read_tsv('data/results/grn/grn_sfari_module_enrichment.tsv')


plot_df <- tf_module_enrich %>% 
    filter(module_type=='combined') %>% 
    arrange(desc(sign(logodds)*-log10(padj))) %>% 
    mutate(tf=factor(tf, levels=unique(.$tf)))


ggplot(plot_df, aes(x=tf, y=sign(logodds)*-log10(padj), fill=padj<0.01)) +
    geom_bar(stat='identity') +  
    scale_fill_manual(values=c('darkgrey', '#AB6472')) +
    rotate_x_text(40) +
    ggtitle('combined')


plot_df <- tf_module_enrich %>% 
    filter(module_type=='combined') %>% 
    arrange(logodds) %>% 
    mutate(tf=factor(tf, levels=unique(.$tf)))

ggplot(plot_df, aes(x=logodds, y=tf, color=padj<0.01&logodds>1, size=ngenes, fill=padj<0.01&logodds>1)) +
    geom_bar(stat='identity', color='black', width=0.1, size=0.1, fill='black') +  
    geom_point(shape=21, color='black') +  
    scale_size_continuous(range=c(0.5,2.5)) +
    scale_fill_manual(values=c('darkgrey', '#AB6472')) +
    article_text() +
    no_legend() +
    theme(
        axis.title.y = element_blank(),
        panel.grid.minor.x = element_line(size=0.1, color='lightgrey'),
        panel.grid.major.x = element_line(size=0.1, color='lightgrey')
    ) +
    labs(x='log2(fold enrichment)')

ggsave('plots/paper/fig3_module_enrich_SFARI_lolli.png', width=4, height=20, bg='white', units='cm')
ggsave('plots/paper/fig3_module_enrich_SFARI_lolli.pdf', width=4, height=20, units='cm')


plot_df %>% filter(padj<0.01&logodds>1) %>% pull(tf) %>% unique() %>% length()
plot_df %>% filter(padj<0.01&logodds>1&tf%in%sfari_genes$`gene-symbol`) %>% pull(tf) %>% unique() %>% length()
plot_df %>% filter(padj<0.01&logodds>1&!tf%in%sfari_genes$`gene-symbol`) %>% pull(tf) %>% unique() %>% length()


ggplot(plot_df, aes(x=logodds, y=tf, color=padj<0.01&logodds>1, size=ngenes, fill=tf%in%sfari_genes$`gene-symbol`)) +
    geom_bar(stat='identity', color='black', width=0.1, size=0.1, fill='black') +  
    geom_point(shape=21, color='black') +  
    scale_size_continuous(range=c(0.5,2.5)) +
    scale_fill_manual(values=c('darkgrey', '#AB6472')) +
    article_text() +
    no_legend() +
    theme(
        axis.title.y = element_blank(),
        panel.grid.minor.x = element_line(size=0.1, color='lightgrey'),
        panel.grid.major.x = element_line(size=0.1, color='lightgrey')
    ) +
    labs(x='log2(fold enrichment)')



#### Full GRN with ASD modules highlighted ####
#### Compute tf coexpression ####
rna_expr <- t(GetAssayData(crop_muo_wt, assay='RNA', slot='data')[unique(grn_modules$tf), ])
gene_cov <- sparse_cov(rna_expr)
gene_cor_sparse <- gene_cov$cor %>%
    # {.[abs(.)<0.1] <- 0; .} %>%
    Matrix::Matrix(sparse=T)

gene_cor_df <- gene_cor_sparse %>%
    as_tibble(rownames='source') %>%
    pivot_longer(!source, names_to='target', values_to='corr')


#### Get adjacency df and matrix ####
coefs_modules_var <- grn_modules %>%
    filter(target%in%.$tf) %>%
    group_by(target)

tf_net <- coefs_modules_var %>%
    select(tf, target, everything()) %>%
    group_by(target) %>%
    left_join(gene_cor_df, by=c('tf'='source', 'target')) %>% {.$corr[is.na(.$corr)] <- 0; .}

reg_mat <- tf_net %>%
    select(target, tf, estimate) %>%
    pivot_wider(names_from=tf, values_from=estimate, values_fill=0) %>%
    as_matrix() %>% Matrix::Matrix(sparse=T)


#### Layout with UMAP on adjacency matrix ####
reg_factor_mat <- abs(reg_mat) + 1

weighted_coex_mat <- gene_cor_sparse[rownames(reg_factor_mat), colnames(reg_factor_mat)] * sqrt(reg_factor_mat)

set.seed(111)
weight_coex_umap <- umap(weighted_coex_mat, n_pcs=30)

weight_coex_meta <- weight_coex_umap %>%
    dplyr::rename('gene'='cell')

tf_graph <- as_tbl_graph(tf_net) %>%
    activate(edges) %>%
    mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
    activate(nodes) %>%
    mutate(
        # dir=sign(estimate),
        central_pr=centrality_pagerank(weights = estimate),
        central_betw=centrality_betweenness(),
        central_eig=centrality_eigen(),
        central_deg=centrality_degree(),
        outdegree=centrality_degree(mode='out'),
        indegree=centrality_degree(mode='in')
    ) %>%
    inner_join(weight_coex_meta, by=c('name'='gene')) %>%
    # inner_join(module_sfari_enrich, by=c('name'='tf')) %>%
    activate(edges) %>%
    filter(padj<1e-2)


ggraph(tf_graph, x=UMAP1, y=UMAP2) +
    geom_edge_diagonal(width=0.1, color='lightgrey') +
    geom_node_point(aes(filter=!name%in%enrich_tfs), size=1, color='darkgrey') +
    geom_node_point(aes(filter=name%in%sfari_genes$`gene-symbol` & name%in%enrich_tfs), size=1, color=orange, fill='black', shape=21, stroke=1) +
    geom_node_point(aes(filter=name%in%sfari_genes$`gene-symbol` & !name%in%enrich_tfs), size=1, color=orange) +
    geom_node_point(aes(filter=name%in%enrich_tfs & !name%in%sfari_genes$`gene-symbol`), size=1, color='black') +
    geom_node_text(
        aes(label=name, size=central_pr, filter=name%in%enrich_tfs | name%in%sfari_genes$`gene-symbol`),
        repel=T, size=5/ggplot2::.pt, max.overlaps=99999
    ) +
    # geom_node_text(aes(label=name), size=3) +
    scale_edge_color_gradientn(colors=rev(rdbu()), limits=c(-6,6)) +
    scale_edge_alpha_continuous(range=c(0.2,0.8), limits=c(2,20)) +
    scale_fill_viridis(option='magma', direction = -1) +
    article_text() +
    theme_void() + no_legend()

ggsave('plots/paper/fig3_sfari_enrich_grn_graph.png', width=7.6, height=6.6, units = 'cm', bg='white')
ggsave('plots/paper/fig3_sfari_enrich_grn_graph.pdf', width=7.6, height=6.6, units = 'cm')



#### Enrichment of individual KO DEG in modules ####

tf_module_enrich <- read_tsv('data/results/grn/grn_sfari_module_enrichment.tsv')
enrich_tfs <- tf_module_enrich %>% 
    filter(logodds>1, padj<0.05) %>% 
    filter(module_type=='combined') %>% 
    pull(tf) %>% unique()


module_enrich <- map_dfr(set_names(c('excitatory', 'inhibitory')), function(lin){
    print(lin)
    enrich_df <- map_dfr(grn_modules_full, function(x){
        enrich_df <- map_dfr(set_names(unique(de_sig$gRNA)), function(grna){
            deg_use <- de_sig %>% 
                filter(lineage==lin, gRNA==grna) %>% 
                pull(feature) %>% unique() %>% intersect(grn_features)
            
            if (length(deg_use)==0){
                return(tibble())
            }
            
            bg_tab <- table(grn_features%in%deg_use)
            x_tab <- table(x%in%deg_use)
            
            test_table <- rbind(
                c(bg_tab['FALSE'], bg_tab['TRUE']),
                c(x_tab['FALSE'], x_tab['TRUE'])
            )
            test_table[is.na(test_table)] <- 0
            ftest <- fisher.test(test_table)
            
            return(tibble(
                pval = ftest$p.value,
                logodds = log(ftest$estimate),
                ngenes = length(x),
                ngenes_de = test_table[[2,'TRUE']]
            ))  
        }, .id='gRNA')
        return(enrich_df)
    }, .id='tf')
    return(enrich_df)
}, .id='lineage')

module_enrich %>% write_tsv('data/results/grn/grn_DEG_module_enrichment.tsv')
module_enrich <- read_tsv('data/results/grn/grn_DEG_module_enrichment.tsv')


plot_df <- module_enrich %>% 
    mutate(padj=p.adjust(pval, method='fdr')) %>% 
    filter(padj<1)

ggplot(plot_df, aes(gRNA, -log10(pval), color=padj<0.05, label=tf)) +
    geom_jitter() +
    geom_text()




#### All guides combined per trajectory ####
module_enrich_all <- map_dfr(set_names(c('excitatory', 'inhibitory')), function(lin){
    print(lin)
    enrich_df <- map_dfr(grn_modules_full, function(x){
        deg_use <- de_sig_top %>% 
            filter(lineage==lin) %>% 
            pull(feature) %>% unique() %>% intersect(grn_features)
        
        if (length(deg_use)==0){
            return(tibble())
        }
        
        bg_tab <- table(grn_features%in%deg_use)
        x_tab <- table(x%in%deg_use)
        
        test_table <- rbind(
            c(bg_tab['FALSE'], bg_tab['TRUE']),
            c(x_tab['FALSE'], x_tab['TRUE'])
        )
        test_table[is.na(test_table)] <- 0
        ftest <- fisher.test(test_table)
        
        return(tibble(
            pval = ftest$p.value,
            logodds = log(ftest$estimate),
            ngenes = length(x),
            ngenes_de = test_table[[2,'TRUE']]
        ))  
    }, .id='tf')
    return(enrich_df)
}, .id='lineage')

module_enrich_all %>% write_tsv('data/results/grn/grn_DEG_all_module_enrichment.tsv')
module_enrich_all <- read_tsv('data/results/grn/grn_DEG_all_module_enrichment.tsv')

module_enrich_all <- module_enrich_all %>% 
    mutate(padj=p.adjust(pval, method='fdr')) %>% 
    mutate(signed_p=sign(logodds)*-log10(padj)) 

module_enrich_split <- module_enrich_all %>% 
    group_by(lineage) %>% group_split() %>% 
    map(
        ~ .x %>% arrange(desc(signed_p)) %>% 
            # filter(signed_p>1.3) %>% 
            mutate(tf=factor(tf, levels=rev(unique(.$tf))))
    )


plot_df <- module_enrich_all %>% 
    mutate(cat=factor(case_when(
        tf%in%sfari_genes$`gene-symbol` & !tf%in%enrich_tfs ~ 'SFARI',
        tf%in%sfari_genes$`gene-symbol` & tf%in%enrich_tfs ~ 'SFARI+enrich',
        tf%in%enrich_tfs & !tf%in%sfari_genes$`gene-symbol` ~ 'enrich',
        T ~ 'none'
    ), levels=c('none', 'SFARI', 'enrich', 'SFARI+enrich'))) %>% 
    mutate(ngenes_bins=cut(ngenes, c(0,50,100,1000)))
    

ggplot(plot_df, aes(cat, signed_p)) +
    geom_boxplot(alpha=0.5, outlier.shape = NA, color='white') +
    geom_quasirandom(data=filter(plot_df, cat=='none'), color='grey', size=0.2) +
    geom_quasirandom(data=filter(plot_df, tf%in%enrich_tfs & !tf%in%sfari_genes$`gene-symbol`), color='black', size=0.2) +
    geom_quasirandom(data=filter(plot_df, tf%in%sfari_genes$`gene-symbol` & !tf%in%enrich_tfs), color=yellow, size=0.2) +
    geom_quasirandom(data=filter(plot_df, tf%in%sfari_genes$`gene-symbol` & tf%in%enrich_tfs), color=yellow, fill='black', shape=21, stroke=0.2, size=0.4) +
    geom_boxplot(alpha=0.7, outlier.shape = NA, size=0.3) +
    theme_rangeframe() + scale_axis_rangeframe() +
    article_text() +
    labs(y='Signed FDR', x='TF group')
ggsave('plots/paper/fig3_DEG_enrich_group_boxplot.pdf', width=3, height=4, units = 'cm')





test_df <- filter(plot_df, cat%in%c('none', 'SFARI+enrich')) %>% 
    mutate(cat=ifelse(cat=='none', cat, 'enrich'))

null_model <- lm(formula = signed_p~ngenes+cat, data=test_df)
summary(null_model)


label_tfs_ex <- c('NEUROD4','EOMES','SOX4','NEUROD6','ETV1','MEF2C')
label_tfs_in <- c('OLIG1','NFIA','DACH1','MEIS2','ARX','SOX6')

pt_size <- 0.2
yellow <- '#FBC02D'
plots <- map(module_enrich_split, 
    ~ggplot(.x, aes(signed_p, tf, label=tf)) +
        # geom_bar(stat='identity', width=0.2, size=0) +
        # geom_point(aes(color=abs(signed_p)>2, size=ngenes_de)) +
        geom_vline(xintercept=1.3, size=0.5, color='grey', linetype='dashed') +
        geom_vline(xintercept=0, size=0.3, color='darkgrey') +
        geom_point(size=pt_size, color='darkgrey') +
        geom_point(data=filter(.x, tf%in%sfari_genes$`gene-symbol` & !tf%in%enrich_tfs), size=0.1, color=yellow) +
        geom_point(data=filter(.x, tf%in%enrich_tfs & !tf%in%sfari_genes$`gene-symbol`), size=pt_size, color='black') +
        geom_point(data=filter(.x, tf%in%sfari_genes$`gene-symbol` & tf%in%enrich_tfs), size=0.4, color=yellow, fill='black', shape=21, stroke=0.4) +
        geom_text_repel(data=filter(.x, tf%in%label_tfs_in, lineage=='inhibitory'), nudge_x = 999, min.segment.length = 0, size=5/ggplot2::.pt, segment.size=0.2,) +
        geom_text_repel(data=filter(.x, tf%in%label_tfs_ex, lineage=='excitatory'), nudge_x = 999, min.segment.length = 0, size=5/ggplot2::.pt, segment.size=0.2,) +
        scale_color_manual(values=c('black', yellow)) +
        scale_axis_rangeframe() + theme_rangeframe() +
        facet_grid(~lineage) +
        article_text() +
        no_legend() + no_y_text() +
        theme(
            plot.margin = unit(rep(0,4), 'lines'),
            axis.line.y = element_blank(),
            strip.text = element_blank()
        ) +
        labs(x='Signed FDR', y='TF module')
)
wrap_plots(plots)

ggsave('plots/paper/fig3_module_de_enrichment_lolli.png', width=6, height=7.6, bg='white', units='cm')
ggsave('plots/paper/fig3_module_de_enrichment_lolli.pdf', width=6, height=7.6, units='cm')


plots <- map(module_enrich_split, function(x){
    plot_df <- filter(x, signed_p>1.3)
    p <- ggplot(plot_df, aes(signed_p, tf)) +
        # geom_bar(stat='identity', width=0.2, size=0) +
        # geom_point(aes(color=abs(signed_p)>2, size=ngenes_de)) +
        geom_point(size=0.2, color='grey') +
        geom_point(data=filter(plot_df, tf%in%sfari_genes$`gene-symbol` & !tf%in%enrich_tfs), size=0.1, color=yellow) +
        geom_point(data=filter(plot_df, tf%in%enrich_tfs & !tf%in%sfari_genes$`gene-symbol`), size=pt_size, color='black') +
        geom_point(data=filter(plot_df, tf%in%sfari_genes$`gene-symbol` & tf%in%enrich_tfs), size=0.4, color=yellow, fill='black', shape=21, stroke=0.4) +
        geom_vline(xintercept = 1.3, size=0.5, color='darkgrey') +
        geom_vline(xintercept=0, size=0.3, color='lightgrey') +
        scale_color_manual(values=c('black', orange)) +
        scale_axis_rangeframe() + theme_rangeframe() +
        facet_grid(~lineage) +
        # article_text() +
        no_legend() +
        labs(x='Signed FDR', y='TF module')
})
wrap_plots(plots)





#### Get gene-centered GRN graphs ####
grn_graph <- grn_modules %>% select(tf, target, everything()) %>% as_tbl_graph()
all_genes <- grn_graph %N>% as_tibble() %>% pull(name)

genes_plot <- c('EOMES', 'OLIG1', 'NEUROD2')

registerDoParallel(16)
gene_graphs <- map(genes_plot, function(gene){
    
    spaths <- all_shortest_paths(grn_graph, gene, all_genes, mode='out')$res
    spath_list <- Pando::map_par(spaths, function(p){
        edg <- names(p)
        edg_graph <- grn_graph %N>%
            filter(name%in%edg) %>%
            convert(to_shortest_path, from=which(.N()$name==edg[1]), to=which(.N()$name==edg[length(edg)])) %E>%
            mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
            as_tibble()
        edg_dir <- edg_graph %>% pull(estimate) %>% sign() %>% prod()
        edg_p <- edg_graph %>% pull(padj) %>% {-log10(.)} %>% mean()
        return(
            list(
                path = tibble(
                    start_node = edg[1],
                    end_node = edg[length(edg)],
                    dir = edg_dir,
                    path = paste(edg, collapse=';'),
                    path_regions = paste(edg_graph$regions, collapse=';'),
                    order = length(edg)-1,
                    mean_padj = edg_p
                ),
                graph = mutate(edg_graph, path=paste(edg, collapse=';'), end_node=edg[length(edg)], comb_dir=edg_dir)
            )
        )
    }, parallel = T)
    
    spath_dir <- map_dfr(spath_list, function(x) x$path) %>%
        mutate(
            path_genes=str_split(path, ';'),
            path_regions=str_split(path_regions, ';')
        )
    spath_graph <- map_dfr(spath_list, function(x) x$graph)
    
    grn_pruned <- spath_dir %>%
        select(start_node, end_node, everything()) %>%
        group_by(end_node) %>% filter(order<=2) %>% 
        filter(order==1 | mean_padj==max(mean_padj))
    
    spath_graph_pruned <- spath_graph %>%
        filter(path%in%grn_pruned$path) %>%
        select(from_node, to_node, end_node, comb_dir) %>% distinct()
    
    grn_graph_pruned <- grn_graph %E>%
        mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
        as_tibble() %>% distinct() %>%
        inner_join(spath_graph_pruned) %>%
        select(from_node, to_node, everything(), -from, -to) %>% arrange(comb_dir) %>% as_tbl_graph()
    
    return(grn_graph_pruned)
})


ggraph(gene_graphs[[1]], layout='tree', circular=T) +
    geom_edge_diagonal(aes(color=sign(estimate)), width=0.1) +
    # geom_node_label(aes(label=name, size=central_pr, filter=name%in%grn_modules$tf), size=3) +
    geom_node_point(aes(filter=!name%in%enrich_tfs), size=1, color='darkgrey') +
    geom_node_point(aes(filter=name%in%enrich_tfs & !name%in%sfari_genes$`gene-symbol`), size=0.7, color='black') +
    geom_node_point(aes(filter=name%in%sfari_genes$`gene-symbol` & !name%in%enrich_tfs), size=0.7, color='#fbc02d') +
    geom_node_point(aes(filter=name%in%sfari_genes$`gene-symbol` & name%in%enrich_tfs), size=0.7, color='#fbc02d', fill='black', shape=21, stroke=1) +
    geom_node_label(aes(label=name, size=central_pr, filter=name%in%c('EOMES', 'BCL11A', 'BCL11B', 'DACH1', 'MEIS2', 'NFIA', 'NEUROD4', 'FOXP2', 'MEF2C', 'RFX4', 'ZBTB20')), size=5/ggplot2::.pt, label.padding=unit(0.1, 'line')) +
    scale_edge_color_gradientn(colors=c('#aed6f1', '#fadbd8')) +
    scale_x_continuous(expand=c(0.1, 0)) +
    scale_y_continuous(expand=c(0.1, 0)) +
    article_text() + theme_void() + no_legend()

ggsave('plots/paper/fig3_EOMES_grn_cirular.png', width=4, height=4, bg='white', units='cm')
ggsave('plots/paper/fig3_EOMES_grn_cirular.pdf', width=4, height=4, units='cm')



ggraph(gene_graphs[[2]], layout='tree', circular=T) +
    geom_edge_diagonal(aes(color=sign(estimate)), width=0.1) +
    # geom_node_label(aes(label=name, size=central_pr, filter=name%in%grn_modules$tf), size=3) +
    geom_node_point(aes(filter=!name%in%enrich_tfs), size=1, color='darkgrey') +
    geom_node_point(aes(filter=name%in%enrich_tfs & !name%in%sfari_genes$`gene-symbol`), size=0.7, color='black') +
    geom_node_point(aes(filter=name%in%sfari_genes$`gene-symbol` & !name%in%enrich_tfs), size=0.7, color='#fbc02d') +
    geom_node_point(aes(filter=name%in%sfari_genes$`gene-symbol` & name%in%enrich_tfs), size=0.7, color='#fbc02d', fill='black', shape=21, stroke=1) +
    geom_node_label(aes(label=name, size=central_pr, filter=name%in%c('RFX4', 'OLIG1', 'SOX6', 'ETV1', 'DACH1', 'NPAS3', 'ZNF521')), size=5/ggplot2::.pt, label.padding=unit(0.1, 'line')) +
    scale_edge_color_gradientn(colors=c('#aed6f1', '#fadbd8')) +
    scale_x_continuous(expand=c(0.1, 0)) +
    scale_y_continuous(expand=c(0.1, 0)) +
    article_text() + theme_void() + no_legend()

ggsave('plots/paper/fig3_OLIG1_grn_cirular.png', width=4, height=4, bg='white', units='cm')
ggsave('plots/paper/fig3_OLIG1_grn_cirular.pdf', width=4, height=4, units='cm')



ggraph(gene_graphs[[3]], layout='tree', circular=T) +
    geom_edge_diagonal(aes(color=sign(estimate)), width=0.3) +
    # geom_node_label(aes(label=name, size=central_pr, filter=name%in%grn_modules$tf), size=3) +
    geom_node_point(aes(filter=!name%in%enrich_tfs), size=1, color='darkgrey') +
    geom_node_point(aes(filter=name%in%enrich_tfs & !name%in%sfari_genes$`gene-symbol`), size=0.7, color='black') +
    geom_node_point(aes(filter=name%in%sfari_genes$`gene-symbol` & !name%in%enrich_tfs), size=0.7, color='#fbc02d') +
    geom_node_point(aes(filter=name%in%sfari_genes$`gene-symbol` & name%in%enrich_tfs), size=0.7, color='#fbc02d', fill='black', shape=21, stroke=1) +
    geom_node_label(aes(label=name, size=central_pr, filter=name%in%c('NEUROD2', 'NFIA', 'SATB2', 'SOX11', 'SOX4', 'NEUROD6', 'BCL11A', 'BHLHE22')), size=5/ggplot2::.pt, label.padding=unit(0.1, 'line')) +
    scale_edge_color_gradientn(colors=c('#aed6f1', '#fadbd8')) +
    scale_x_continuous(expand=c(0.1, 0)) +
    scale_y_continuous(expand=c(0.1, 0)) +
    article_text() + theme_void() + no_legend()

ggsave('plots/paper/fig3_NEUROD2_grn_cirular.png', width=4, height=4, bg='white', units='cm')
ggsave('plots/paper/fig3_NEUROD2_grn_cirular.pdf', width=4, height=4, units='cm')




#### ARID1B enriched GRN ####
arid_tfs <- module_enrich %>% 
    mutate(padj=p.adjust(pval, method='fdr')) %>% 
    filter(padj<0.05, gRNA=='ARID1B', lineage=='inhibitory') %>% top_n(4, -padj) %>% pull(tf) 

arid_modules <- grn_modules %>% filter(tf%in%arid_tfs) %>% 
    group_by(tf) %>% group_split()

x <- arid_modules[[1]]
this_tf <- x$tf[1]
ggraph(x, layout='fr') +
    geom_edge_diagonal(width=0.1, color='grey') +
    # geom_node_label(aes(label=name, size=central_pr), size=3) +
    geom_node_point(aes(filter=!name%in%enrich_tfs), size=1, color='darkgrey') +
    geom_node_point(aes(filter=name%in%enrich_tfs & !name%in%sfari_genes$`gene-symbol`), size=1, color='black') +
    geom_node_point(aes(filter=name%in%sfari_genes$`gene-symbol` & !name%in%enrich_tfs), size=1, color='#fbc02d') +
    geom_node_point(aes(filter=name%in%sfari_genes$`gene-symbol` & name%in%enrich_tfs), size=1, color='#fbc02d', fill='black', shape=21, stroke=1) +
    geom_node_label(aes(label=name, size=central_pr, filter=name==this_tf), size=5/ggplot2::.pt, label.padding=unit(0.1, 'line')) +
    theme_void()
ggsave('plots/paper/fig3_BCL11A_module_fr.pdf', width=2, height=2, units='cm')


x <- arid_modules[[4]]
this_tf <- x$tf[1]
ggraph(x, layout='fr') +
    geom_edge_diagonal(width=0.1, color='grey') +
    # geom_node_label(aes(label=name, size=central_pr), size=3) +
    geom_node_point(aes(filter=!name%in%enrich_tfs), size=1, color='darkgrey') +
    geom_node_point(aes(filter=name%in%enrich_tfs & !name%in%sfari_genes$`gene-symbol`), size=1, color='black') +
    geom_node_point(aes(filter=name%in%sfari_genes$`gene-symbol` & !name%in%enrich_tfs), size=1, color='#fbc02d') +
    geom_node_point(aes(filter=name%in%sfari_genes$`gene-symbol` & name%in%enrich_tfs), size=1, color='#fbc02d', fill='black', shape=21, stroke=1) +
    geom_node_label(aes(label=name, size=central_pr, filter=name==this_tf), size=5/ggplot2::.pt, label.padding=unit(0.1, 'line')) +
    theme_void()
ggsave('plots/paper/fig3_SATB2_module_fr.pdf', width=3.5, height=3.5, units='cm')




















