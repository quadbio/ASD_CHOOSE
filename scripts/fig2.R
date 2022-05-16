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
# crop %>% write_rds('data/objects/CROP/asd_chong_annot_withguides_v3.3modules_srt.rds')


#### Format meta ####
crop$state <- case_when(
    crop$celltype_cl_coarse2 %in% c('mesenchyme', 'Astrocytes') ~ 'other',
    str_detect(crop$celltype_cl_coarse2, 'RG|INP|IP|OPC') ~ 'progenitor',
    T ~ 'neuron'
)

crop$lineage <- case_when(
    crop$celltype_cl_coarse2 %in% asd_celltype_names2[1:8] ~ 'excitatory',
    crop$celltype_cl_coarse2 %in% asd_celltype_names2[9:16] ~ 'inhibitory',
    T ~ 'other'
)

crop$celltype_cl_coarse2 <- case_when(
    crop$clusters %in% c(35) ~ 'L6_CThPN',
    T ~ crop$celltype_cl_coarse
)




#### Detemine level to tast at ####
crop$test_celltype <- crop$celltype_cl_coarse2


meta <- crop@meta.data %>% 
    as_tibble(rownames='cell') %>% 
    inner_join(as_tibble(crop[['umap']]@cell.embeddings, rownames='cell'))


#### Do enrichment tests ####
dim_plot(crop, group.by=c('test_celltype', 'clusters_ge'), label=T)
feature_plot(crop, features=c('DLX5', 'SCGN', 'DLX2', 'GAD2'), order=T)

fisher_enrich <- test_guide_enrichment(
    crop, 
    test_groups = 'test_celltype', 
    guide_assay = 'guide_assignments_ctrl',
    nt_name = 'Control2',
    groups = 'library'
)

constistency <- fisher_enrich %>% 
    mutate(padj=p.adjust(pval, method='fdr')) %>% 
    filter(padj<0.05) %>% 
    group_by(x, y) %>% 
    mutate(dir=sign(log_odds_ratio)) %>% 
    mutate(const=max(dir)-min(dir)<2) %>% 
    distinct(x, y, const)


cmh_enrich <- test_guide_enrichment(
    crop,
    method = 'cmh',
    test_groups = 'test_celltype', 
    guide_assay = 'guide_assignments_ctrl',
    nt_name = 'Control2',
    groups = 'library'
)

lin_cmh_enrich <- test_guide_enrichment(
    crop,
    method = 'cmh',
    test_groups = 'lineage', 
    guide_assay = 'guide_assignments_ctrl',
    nt_name = 'Control2',
    groups = 'library'
)



# cmh_enrich %>% write_tsv('data/results/enrichment/celltype_coarse_cmh_enrich.tsv')
# fisher_enrich %>% write_tsv('data/results/enrichment/celltype_coarse_fisher_enrich.tsv')

state_meta <- crop@meta.data %>% 
    as_tibble() %>% distinct(test_celltype, lineage, state)
    
plot_df <- cmh_enrich %>% 
    mutate(padj=p.adjust(pval, method = 'fdr')) %>% 
    left_join(constistency) %>% 
    inner_join(state_meta, by=c('y'='test_celltype')) %>% 
    filter(y!='mesenchyme') %>% 
    mutate(log_padj=pmin(-log10(padj), 5)) %>% 
    mutate(signed_log_padj=sign(log_odds_ratio)*(log_padj)) %>% 
    mutate(y=factor(y, levels=names(asd_celltype_colors2))) 

lodds_mat <- plot_df %>% 
    distinct(x,y,signed_log_padj) %>% 
    pivot_wider(names_from = x, values_from=signed_log_padj, values_fill=0) %>% 
    column_to_rownames('y') %>% as.matrix()

gene_clust <- lodds_mat %>% t() %>% dist() %>% hclust(method='ward.D')
gene_order <- gene_clust %>% {.$labels[.$order]}

p_ct <- ggplot(plot_df, aes(y, factor(x, levels=gene_order), fill=signed_log_padj)) +
    geom_tile() +
    # geom_point(data=filter(plot_df, const), size=1, color='darkgrey') +
    geom_point(data=filter(plot_df, padj<5e-2), size=0.8, fill='black', shape=21, color='white', stroke=0.4) +
    scale_fill_gradientn(colors=rev(pals::brewer.rdbu(100)), limits=c(-5,5)) +
    rotate_x_text(45) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    no_legend() +
    article_text() +
    theme(
        axis.ticks = element_line(size=0.3),
        panel.border = element_rect(size=0.2),
        plot.margin = unit(rep(0.05,4), 'line'),
        axis.title.y = element_blank()
    ) +
    facet_grid(~lineage, space='free', scales='free') +
    labs(x='Celltype', y='KO gene') +
    no_y_text()

p_ct



#### Heatmap strip with only DV ####

# cmh_enrich %>% write_tsv('data/results/enrichment/celltype_coarse_cmh_enrich.tsv')
# fisher_enrich %>% write_tsv('data/results/enrichment/celltype_coarse_fisher_enrich.tsv')

plot_df <- lin_cmh_enrich %>% 
    mutate(padj=p.adjust(pval, method = 'fdr')) %>% 
    filter(y!='other', y!='inhibitory') %>% 
    mutate(log_padj=pmin(-log10(padj), 5)) %>% 
    mutate(signed_log_padj=sign(log_odds_ratio)*(log_padj)) 

p_dv <- ggplot(plot_df, aes(y, factor(x, levels=gene_order), fill=signed_log_padj)) +
    geom_tile() +
    # geom_point(data=filter(plot_df, const), size=1, color='darkgrey') +
    geom_point(data=filter(plot_df, padj<5e-2), size=0.8, fill='black', shape=21, color='white', stroke=0.4) +
    scale_fill_gradientn(colors=rev(pals::brewer.puor(100))[10:90], limits=c(-5,5)) +
    rotate_x_text(45) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    no_legend() +
    article_text() +
    theme(
        axis.ticks = element_line(size=0.3),
        panel.border = element_rect(size=0.2),
        plot.margin = unit(rep(0.05,4), 'line'),
        axis.title.x = element_blank()
    ) +
    labs(y='Perturbed gene') +
    no_x_text() 

(p_dv | p_ct) + plot_layout(widths=c(1,15))

ggsave('plots/paper/fig2_enrichment_heatmap.pdf', width=6, height=9.4, unit='cm')
ggsave('plots/paper/fig2_enrichment_heatmap.png', width=6, height=9.4, bg='white',unit='cm')


print_scale(rev(pals::brewer.rdbu(100)))
ggsave('plots/paper/fig2_enrichment_heatmap_scale.pdf', width=2, height=2)

print_scale(rev(pals::brewer.puor(100)))
ggsave('plots/paper/fig2_enrichment_heatmap_scale2.pdf', width=2, height=2)



#### Pseudotime enrichment ####
#### Inhib trajectory ####
#### Over Inhib PT ####
ge_meta <- crop@meta.data %>% 
    filter(lineage=='inhibitory', !is.na(pseudotime_ranks)) %>% 
    as_tibble(rownames='cell') 

ge_density_meta <- ge_meta %>% 
    group_by(gRNA) %>%
    nest() %>% 
    ungroup() %>%  
    mutate(density =  map(data, ~density(.x$pseudotime_ranks))) %>% 
    hoist(density, "x", "y") %>% 
    select(-density, -data) %>% 
    unnest(c(x,y)) %>% filter(x>-0.1, x<1.1) %>% 
    mutate(pt_bins=cut(x,seq(-0.1,1.1,0.01))) %>% 
    group_by(gRNA, pt_bins) %>% 
    summarise(density=mean(y), x=mean(x)) %>% 
    group_by(gRNA) 

ctrl_dens <- ge_density_meta %>% filter(gRNA=='Control2') %>% pull(density)

ge_diff_meta <- ge_density_meta %>% 
    mutate(diff_density=density-ctrl_dens) %>% 
    mutate(diff_sign=sign(diff_density)) %>% 
    filter(gRNA!='Control2')

p <- ggplot(ge_diff_meta, aes(pt_bins, diff_density, fill=gRNA)) +
    geom_bar(color='black', stat='identity', size=0.1, alpha=0.7, width=1) +
    # scale_fill_manual(values=grna_scale) +
    article_text() +
    no_x_text() +
    theme_rangeframe() +
    scale_axis_rangeframe() +
    labs(x='Pseudotime bins', y='Differential density vs control')

p
ggsave('plots/paper/fig2_grna_legend.pdf', width=4, height=4)
p + no_legend()
ggsave('plots/paper/fig2_inhib_diff_density.pdf', width=6, height=2.2, unit='cm')
ggsave('plots/paper/fig2_inhib_diff_density.png', width=6, height=2.2, unit='cm')




ggplot(ge_diff_meta, aes(pt_bins, diff_density)) +
    geom_bar(color='black', stat='identity', size=0.1, alpha=0.7, width=1) +
    facet_wrap(~gRNA) +
    article_text() +
    no_x_text() +
    theme_rangeframe() +
    scale_axis_rangeframe() +
    theme(
        axis.line.x = element_blank(),
        panel.spacing.y = unit(0, 'cm')
    ) +
    labs(x='Pseudotime bins', y='Differential density vs control')

ggsave('plots/paper/fig2_inhib_diff_density_split.pdf', width=9, height=7, unit='cm')
ggsave('plots/paper/fig2_inhib_diff_density_split.png', width=9, height=7, unit='cm', bg='white')



#### Over Excit PT ####
ctx_meta <- crop@meta.data %>% 
    filter(lineage=='excitatory', !is.na(pseudotime_ranks)) %>% 
    as_tibble(rownames='cell') 

ctx_density_meta <- ctx_meta %>% 
    group_by(gRNA) %>%
    nest() %>% 
    ungroup() %>%  
    mutate(density =  map(data, ~density(.x$pseudotime_ranks))) %>% 
    hoist(density, "x", "y") %>% 
    select(-density, -data) %>% 
    unnest(c(x,y)) %>% filter(x>-0.1, x<1.1) %>% 
    mutate(pt_bins=cut(x,seq(-0.1,1.1,0.01))) %>% 
    group_by(gRNA, pt_bins) %>% 
    summarise(density=mean(y), x=mean(x)) %>% 
    group_by(gRNA) 

ctrl_dens <- ctx_density_meta %>% filter(gRNA=='Control2') %>% pull(density)

ctx_diff_meta <- ctx_density_meta %>% 
    mutate(diff_density=density-ctrl_dens) %>% 
    mutate(diff_sign=sign(diff_density)) %>% 
    filter(gRNA!='Control2')

ggplot(ctx_diff_meta, aes(pt_bins, diff_density, fill=gRNA)) +
    geom_bar(color='black', stat='identity', size=0.1, alpha=0.7, width=1) +
    # scale_fill_manual(values=grna_scale) +
    article_text() +
    no_x_text() +
    theme_rangeframe() +
    scale_axis_rangeframe() +
    no_legend() +
    labs(x='Pseudotime bins', y='Differential density vs control')

ggsave('plots/paper/fig2_excit_diff_density.pdf', width=6, height=2.2, unit='cm')
ggsave('plots/paper/fig2_excit_diff_density.png', width=6, height=2.2, unit='cm')



ggplot(ctx_diff_meta, aes(pt_bins, diff_density)) +
    geom_bar(color='black', stat='identity', size=0.1, alpha=0.7, width=1) +
    facet_wrap(~gRNA) +
    article_text() +
    no_x_text() +
    theme_rangeframe() +
    scale_axis_rangeframe() +
    theme(
        axis.line.x = element_blank(),
        panel.spacing.y = unit(0, 'cm')
    ) +
    labs(x='Pseudotime bins', y='Differential density vs control')

ggsave('plots/paper/fig2_excit_diff_density_split.pdf', width=9, height=7, unit='cm')
ggsave('plots/paper/fig2_excit_diff_density_split.png', width=9, height=7, unit='cm', bg='white')




#### Expression over pseudotime ####
genes_plot <- c('PAX6', 'VIM', 'GFAP', 'ASCL1', 'SOX2',
    'NEUROG1', 'NEUROG2', 'EOMES', 'CTIP2', 'SATB2',
    'NEROD2', 'NEUROD6', 'NRXN3', 'CUX1','PROX1',
    'DLX2', 'DL5', 'GAD1', 'GAD2', 'OLIG1', 'OLIG2', 'MEIS2','NR2F2', 'BCL11B', 'NEUROD2') %>% intersect(rownames(crop))

gene_expr <- crop[['RNA']]@data[genes_plot, ] 
gene_expr_df <- gene_expr %>% as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='cell', values_to='expr') %>% 
    group_by(gene) %>% mutate(scale_expr=scale01(expr))

ge_meta_expr <- ge_meta %>% 
    inner_join(gene_expr_df)

ge_meta_expr$pt_bins_10 <- cut(ge_meta_expr$pseudotime_ranks, 10, labels=c(1:10))

p1 <- ggplot(ge_meta_expr, aes(pseudotime_ranks, scale_expr)) +
    geom_point(color='lightgrey', alpha=0.05) +
    geom_smooth(color='black') +
    facet_wrap(~gene) +
    scale_x_continuous(breaks=seq(0,1,0.1)) +
    theme(
        panel.grid.major.x = element_line(color='grey'),
        panel.grid.minor.x = element_line(color='grey')
    )

p1
ggsave('plots/CROP/asd_ge_pt_expression_smooth.png', width=12, height=10, bg='white')

ctx_meta_expr <- ctx_meta %>% 
    inner_join(gene_expr_df)

ctx_meta_expr$pt_bins_10 <- cut(ctx_meta_expr$pseudotime_ranks, 10, labels=c(1:10))

p1 <- ggplot(filter(ctx_meta_expr, gene=='NEUROD2'), aes(pseudotime_ranks, expr)) +
    geom_point(color='lightgrey', alpha=0.05) +
    geom_smooth(color='black') +
    facet_wrap(~gene) +
    scale_x_continuous(breaks=seq(0,1,0.1)) +
    theme(
        panel.grid.major.x = element_line(color='grey'),
        panel.grid.minor.x = element_line(color='grey')
    )
p1 
ggsave('plots/CROP/asd_ctx_pt_expression_smooth.png', width=12, height=10, bg='white')


#### Cluster pseudotime bins ####
ge_meta$pt_bins_10 <- cut(ge_meta$pseudotime_ranks, 10, labels=c(1:10))
ge_gene_expr <- crop[['RNA']]@data[VariableFeatures(crop), crop$lineage=='inhibitory' & !is.na(crop$pseudotime_ranks)]
ge_bin_gene_expr <- t(ge_gene_expr) %>% 
    aggregate_matrix(groups=ge_meta$pt_bins_10)

ge_bin_clust <- ge_bin_gene_expr %>% dist() %>% hclust()
ge_bin_clust_tbl <- ge_bin_clust %>% as.dendrogram() %>% as_tbl_graph()

ggraph(ge_bin_clust_tbl) +
    geom_edge_elbow() +
    geom_node_label(aes(label=label, filter=leaf)) +
    theme_void()

ge_hclusters <- cutree(ge_bin_clust, h=10) %>% 
    enframe('pt_bins_10', 'bin_cluster')

ge_meta <- ge_meta %>% 
    inner_join(ge_hclusters)

p1 <- ggplot(ge_meta, aes(factor(pt_bins_10, levels=sort(unique(as.numeric(pt_bins_10)))), '1',  fill=factor(bin_cluster))) +
    geom_tile() +
    scale_fill_manual(values=pals::brewer.greys(5)[2:5]) +
    theme_void() + no_label() + no_legend()

p2 <- ggplot(ge_meta, aes(factor(pt_bins_10, levels=sort(unique(as.numeric(pt_bins_10)))), fill=celltype_cl_coarse2)) +
    geom_bar() +
    scale_fill_manual(values=asd_celltype_colors2) +
    theme_minimal() + no_legend() + no_label() + no_y_text()

p_ge <- p1 / p2



ctx_meta$pt_bins_10 <- cut(ctx_meta$pseudotime_ranks, 10, labels=c(1:10))
ctx_gene_expr <- crop[['RNA']]@data[VariableFeatures(crop), crop$lineage=='excitatory' & !is.na(crop$pseudotime_ranks)]
ctx_bin_gene_expr <- t(ctx_gene_expr) %>% 
    aggregate_matrix(groups=ctx_meta$pt_bins_10)

ctx_bin_clust <- ctx_bin_gene_expr %>% dist() %>% hclust()
ctx_bin_clust_tbl <- ctx_bin_clust %>% as.dendrogram() %>% as_tbl_graph()

ggraph(ctx_bin_clust_tbl) +
    geom_edge_elbow() +
    geom_node_label(aes(label=label, filter=leaf)) +
    theme_void()

ctx_hclusters <- cutree(ctx_bin_clust, h=8) %>% 
    enframe('pt_bins_10', 'bin_cluster')

ctx_meta <- ctx_meta %>% 
    inner_join(ctx_hclusters)

p1 <- ggplot(ctx_meta, aes(factor(pt_bins_10, levels=sort(unique(as.numeric(pt_bins_10)))), '1',  fill=factor(bin_cluster))) +
    geom_tile() +
    scale_fill_manual(values=pals::brewer.greys(6)[2:6]) +
    theme_void() + no_label() + no_legend()

p2 <- ggplot(ctx_meta, aes(factor(pt_bins_10, levels=sort(unique(as.numeric(pt_bins_10)))), fill=celltype_cl_coarse2)) +
    geom_bar() +
    scale_fill_manual(values=asd_celltype_colors2) +
    theme_minimal() + no_legend() + no_label() + no_y_text()

p_ctx <- p1 / p2


(p_ctx) / (p_ge) + plot_layout(heights=c(1,1,2))



ctx_meta$pt_stages_manual <- case_when(
    ctx_meta$pt_bins_10 %in% c(1) ~ 'prolif',
    ctx_meta$pt_bins_10 %in% c(2,3,4) ~ 'diff',
    ctx_meta$pt_bins_10 %in% c(5,6,7) ~ 'deep_spec',
    ctx_meta$pt_bins_10 %in% c(8,9,10) ~ 'upper_spec'
)

ge_meta$pt_stages_manual <- case_when(
    ge_meta$pt_bins_10 %in% c(1,2) ~ 'prolif',
    ge_meta$pt_bins_10 %in% c(3,4,5) ~ 'diff',
    ge_meta$pt_bins_10 %in% c(6,7) ~ 'LGE_spec',
    ge_meta$pt_bins_10 %in% c(8,9,10) ~ 'CGE_spec'
)


#### Expression heatmap over bins ####

ge_genes_plot <- c('GFAP', 'VIM', 'OLIG2', 'ASCL1', 'DLX2', 'PROX1', 'NRXN3','DLX5', 'GAD1','MEIS2', 'NR2F2') %>% intersect(rownames(crop))

ge_gene_expr <- crop[['RNA']]@data[ge_genes_plot, ge_meta$cell] 
ge_gene_expr_bins <- t(ge_gene_expr) %>% aggregate_matrix(groups=ge_meta$pt_bins_10)

ge_gene_expr_df <- ge_gene_expr_bins %>% as_tibble(rownames='pt_bin') %>% 
    pivot_longer(!pt_bin, names_to='gene', values_to='expr') %>% 
    group_by(gene) %>% mutate(scale_expr=scale01(expr))

gene_order <- ge_gene_expr_bins %>% t() %>% dist() %>% hclust() %>% {.$label[.$order]}

plot_df <- ge_gene_expr_df %>% 
    mutate(gene=factor(gene, levels=rev(ge_genes_plot)))

ggplot(plot_df, aes(as.numeric(pt_bin), gene, fill=scale_expr)) +
    geom_tile() +
    scale_fill_gradientn(colors=pals::brewer.greys(100)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    article_text() +
    no_x_text() + no_legend() +
    theme(axis.title.x = element_blank()) +
    labs(y='Gene')

ggsave('plots/paper/fig2_ge_pt_bin_expr_heatmap.pdf', width=6, height=2.2, unit='cm')
ggsave('plots/paper/fig2_ge_pt_bin_expr_heatmap.png', width=6, height=2.2, unit='cm')


ctx_genes_plot <- c('VIM', 'SOX2', 'PAX6', 'NEUROG2', 'EOMES', 'BCL11B', 'NEUROD6', 'NEUROD2','SATB2') %>% intersect(rownames(crop))

ctx_gene_expr <- crop[['RNA']]@data[ctx_genes_plot, ctx_meta$cell] 
ctx_gene_expr_bins <- t(ctx_gene_expr) %>% aggregate_matrix(groups=ctx_meta$pt_bins_10)

ctx_gene_expr_df <- ctx_gene_expr_bins %>% as_tibble(rownames='pt_bin') %>% 
    pivot_longer(!pt_bin, names_to='gene', values_to='expr') %>% 
    group_by(gene) %>% mutate(scale_expr=scale01(expr))

gene_order <- ctx_gene_expr_bins %>% t() %>% dist() %>% hclust() %>% {.$label[.$order]}

plot_df <- ctx_gene_expr_df %>% 
    mutate(gene=factor(gene, levels=rev(ctx_genes_plot)))

ggplot(plot_df, aes(as.numeric(pt_bin), gene, fill=scale_expr)) +
    geom_tile() +
    scale_fill_gradientn(colors=pals::brewer.greys(100)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    article_text() +
    no_x_text() + no_legend() +
    theme(axis.title.x = element_blank()) +
    labs(y='Gene')

ggsave('plots/paper/fig2_ctx_pt_bin_expr_heatmap.pdf', width=6.4, height=2.2, unit='cm')
ggsave('plots/paper/fig2_ctx_pt_bin_expr_heatmap.png', width=6.4, height=2.2, unit='cm')



#### Sig test for (coarser) pt bins ####
ge_meta_use <- filter(ge_meta, !is.na(pt_stages_manual))
ctx_meta_use <- filter(ctx_meta, !is.na(pt_stages_manual))

crop_ge <- crop %>% subset(cells=ge_meta_use$cell)
crop_ctx <- crop %>% subset(cells=ctx_meta_use$cell)

crop_ge$pt_stages_manual <- ge_meta_use$pt_stages_manual
crop_ctx$pt_stages_manual <- ctx_meta_use$pt_stages_manual

ge_pt_cmh_enrich <- test_guide_enrichment(
    crop_ge,
    method = 'cmh',
    test_groups = 'pt_stages_manual', 
    guide_assay = 'guide_assignments_ctrl',
    nt_name = 'Control2',
    groups = 'library'
)

ctx_pt_cmh_enrich <- test_guide_enrichment(
    crop_ctx,
    method = 'cmh',
    test_groups = 'pt_stages_manual', 
    guide_assay = 'guide_assignments_ctrl',
    nt_name = 'Control2',
    groups = 'library'
)

ge_pt_fisher_enrich <- test_guide_enrichment(
    crop_ge, 
    test_groups = 'pt_stages_manual', 
    guide_assay = 'guide_assignments_ctrl',
    nt_name = 'Control2',
    groups = 'library'
)

ctx_pt_fisher_enrich <- test_guide_enrichment(
    crop_ctx, 
    test_groups = 'pt_stages_manual', 
    guide_assay = 'guide_assignments_ctrl',
    nt_name = 'Control2',
    groups = 'library'
)

ge_constistency <- ge_pt_fisher_enrich %>% 
    filter(pval<0.05) %>%
    group_by(x, y) %>% 
    mutate(dir=sign(log_odds_ratio)) %>% 
    mutate(const=max(dir)-min(dir)<2) %>% 
    distinct(x, y, const)

ctx_constistency <- ctx_pt_fisher_enrich %>% 
    filter(pval<0.05) %>%
    group_by(x, y) %>% 
    mutate(dir=sign(log_odds_ratio)) %>% 
    mutate(const=max(dir)-min(dir)<2) %>% 
    distinct(x, y, const)


plot_df <- ge_pt_cmh_enrich %>% 
    left_join(ge_constistency) %>%
    # mutate(y=factor(y, levels=names(asd_celltype_colors2))) %>%
    mutate(y=factor(y, levels=c('prolif', 'diff', 'LGE_spec', 'CGE_spec'))) %>%
    # mutate(padj=p.adjust(pval, method='fdr')) %>% 
    mutate(padj=pval) %>%
    mutate(log_padj=pmin(-log10(padj), 5)) %>% 
    mutate(signed_log_padj=sign(log_odds_ratio)*(log_padj)) %>% 
    group_by(x) %>% 
    filter(any(pval<0.05))

lodds_mat <- plot_df %>% 
    select(x,y,signed_log_padj) %>% 
    pivot_wider(names_from = x, values_from=signed_log_padj, values_fill=0) %>% 
    column_to_rownames('y') %>% as.matrix()

gene_clust <- lodds_mat %>% t() %>% dist() %>% hclust(method='ward.D2')
gene_order <- gene_clust %>% {.$labels[.$order]}

p1 <- ggplot(plot_df, aes(y, factor(x, levels=gene_order), fill=signed_log_padj)) +
    geom_tile() +
    # geom_point(data=filter(plot_df, const), size=3, color='darkgrey') +
    geom_point(data=filter(plot_df, padj<5e-2), size=0.8, fill='black', shape=21, color='white', stroke=0.4) +
    scale_fill_gradientn(colors=rev(pals::brewer.rdbu(100)), limits=c(-5,5)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    labs(x='Pseudotime bin', y='Target gene') +
    article_text() +
    no_label() + no_legend() +
    theme(plot.margin = unit(rep(0,4), 'lines')) 

p2 <- ggplot(as_tibble(crop_ge@meta.data), aes(pt_bins_10, fill=celltype_cl_coarse)) +
    geom_bar(position='fill') +
    scale_fill_manual(values=asd_celltype_colors2) +
    theme_void() + no_legend() +
    theme(plot.margin = unit(rep(0,4), 'lines'))

p1 
ggsave('plots/paper/fig2_inhib_enrich_bins.pdf', width=5.5, height=4, unit='cm')
ggsave('plots/paper/fig2_inhib_enrich_bins.png', width=5.5, height=4, unit='cm')






plot_df <- ctx_pt_cmh_enrich %>% 
    left_join(ctx_constistency) %>%
    # mutate(y=factor(y, levels=names(asd_celltype_colors2))) %>%
    mutate(y=factor(y, levels=c('prolif', 'diff', 'deep_spec', 'upper_spec'))) %>%
    # mutate(padj=p.adjust(pval, method='fdr')) %>% 
    mutate(padj=pval) %>%
    mutate(log_padj=pmin(-log10(padj), 5)) %>% 
    mutate(signed_log_padj=sign(log_odds_ratio)*(log_padj)) %>% 
    group_by(x) %>% 
    filter(any(pval<0.05))

lodds_mat <- plot_df %>% 
    select(x,y,signed_log_padj) %>% 
    pivot_wider(names_from = x, values_from=signed_log_padj, values_fill=0) %>% 
    column_to_rownames('y') %>% as.matrix()

gene_clust <- lodds_mat %>% t() %>% dist() %>% hclust(method='ward.D2')
gene_order <- gene_clust %>% {.$labels[.$order]}

p1 <- ggplot(plot_df, aes(y, factor(x, levels=gene_order), fill=signed_log_padj)) +
    geom_tile() +
    # geom_point(data=filter(plot_df, const), size=3, color='darkgrey') +
    geom_point(data=filter(plot_df, padj<5e-2), size=0.8, fill='black', shape=21, color='white', stroke=0.4) +
    scale_fill_gradientn(colors=rev(pals::brewer.rdbu(100)), limits=c(-5,5)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    labs(x='Pseudotime bin', y='Target gene') +
    article_text() +
    no_label() + no_legend() +
    theme(plot.margin = unit(rep(0,4), 'lines')) 

p2 <- ggplot(as_tibble(crop_ctx@meta.data), aes(pt_bins_6, fill=celltype_cl_coarse)) +
    geom_bar(position='fill') +
    scale_fill_manual(values=asd_celltype_colors2) +
    theme_void() + no_legend() +
    theme(plot.margin = unit(rep(0,4), 'lines'))

p1 
ggsave('plots/paper/fig2_excit_enrich_bins.pdf', width=5.5, height=2, unit='cm')
ggsave('plots/paper/fig2_excit_enrich_bins.png', width=5.5, height=2, unit='cm')





#### Density plots for examples ####
knn <- RANN::nn2(Embeddings(crop, 'pca')[,1:10], k = 200)
knngraph <- FindNeighbors(Embeddings(crop, 'pca')[,1:10], k.param = 200)

grnas_plot <- sort(unique(crop$gRNA))

plots <- map(grnas_plot, function(grna){
    knn_enrich <- get_knn_cmh_enrichment(knn, crop$gRNA == grna, crop$gRNA == 'Control2', stratums = crop$library)
    smooth_enrich <- random_walk_with_restart(knn_enrich, knn_mat = as(knngraph$snn, 'Matrix'), num_rounds = 100, alpha = 0.1)
    
    plot_df <- meta %>% 
        inner_join(enframe(smooth_enrich, 'cell', 'enrich')) %>%
        filter(UMAP_1>-7)
    
    p <- ggplot(plot_df, aes(UMAP_1, UMAP_2, z=enrich)) +
        stat_summary_hex(bins=50, fun='mean') +
        # geom_point(data=filter(plot_df, gRNA==grna), alpha=0.5, color='black', size=0.7) +
        ggtitle(grna) +
        theme_void() +
        no_legend()
    
    pb <- ggplot_build(p)
    clim <- max(abs(pb$data[[1]]$value))
    p + scale_fill_gradientn(colors=rev(bigrad(pals::brewer.rdbu, bias=1.5)), limits=c(-clim,clim)) 
    
})


wrap_plots(plots)
ggsave('plots/paper/fig2_examples_hexumap.png', width=25, height=20)
ggsave('plots/paper/fig2_examples_hexumap.pdf', width=25, height=20)


enrich_list <- map(grnas_plot, function(grna){
    print(grna)
    knn_enrich <- get_knn_cmh_enrichment(knn, crop$gRNA == grna, crop$gRNA == 'Control2', stratums = crop$library)
    smooth_enrich <- random_walk_with_restart(knn_enrich, knn_mat = as(knngraph$snn, 'Matrix'), num_rounds = 100, alpha = 0.1)
    return(list(knn=knn_enrich, smooth=smooth_enrich))
})

names(enrich_list) <- grnas_plot

plots <- map(grnas_plot, function(grna){
    smooth_enrich <- enrich_list[[grna]]$knn
    knn_enrich <- enrich_list[grna]$smooth
    
    plot_df <- meta %>% 
        inner_join(enframe(smooth_enrich, 'cell', 'enrich')) %>%
        filter(UMAP_1>-7)
    
    clim <- max(abs(plot_df$enrich))
    
    p <- ggplot(plot_df, aes(UMAP_1, UMAP_2, z=enrich)) +
        stat_summary_hex(bins=50) +
        geom_point(data=filter(plot_df, gRNA==grna), alpha=0.5, color='black', size=0.01, shape=16) +
        scale_fill_gradientn(colors=rev(bigrad(pals::brewer.rdbu, bias=2)), limits=c(-clim,clim)) +
        ggtitle(grna) +
        theme_void() +
        no_legend()
})

wrap_plots(plots) & theme(title = element_text(size=7))
ggsave('plots/paper/fig2_examples_hexumap2.png', width=18, height=14, units = 'cm')
ggsave('plots/paper/fig2_examples_hexumap2.pdf', width=18, height=14, units = 'cm')






