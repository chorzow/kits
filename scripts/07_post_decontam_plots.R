library(ggplot2)
library(ggsci)
library(ggh4x)
library(vegan)
library(patchwork)
library(scales)
library(microshades)
library(cowplot)

workdir <- '/Users/asobolev/skoltech/articles/kits'
source(file.path(workdir, 'scripts/01_preload.R'))
method <- 'either'  # used for decontam
threshold <- 0.2  # used for decontam
selected_groups_control <- c("Bacteroidota", "Cyanobacteria", "Firmicutes", 
                             "Actinobacteriota", "Proteobacteria")
decontam_otu_path <- file.path(workdir, 
                               paste0('results/decontam/OTUs_decontam_', method,
                                      '_', threshold, 
                                      '_grouped_all_vs_all_2.tsv'))

# rarefaction for controls ------------------------------------------------
control_s <- rbind(kitome_s, splashome_s)
cur_samples <- control_s

sdata <- sample_data(cur_samples)
ps <- phyloseq(OTU, TAX, sdata)
tab <- as.matrix(as.data.frame(otu_table(ps)))
rare_c <- rarecurve(t(tab), step=500, cex=0.5, col = 'steelblue4', label = T)


# pcoa for controls -------------------------------------------------------

carbom_abund <- get_relabund(OTU, TAX, sdata, abund = 1)

set.seed(42)
pcoa_c12 <- pcoa(carbom_abund, color = 'Biome', shapes = 'Group', alpha = 0.8)
pcoa_c13 <- pcoa(carbom_abund, nmds = c(1, 3),
                 color = 'Biome', shapes = 'Group', alpha = 0.8)
pcoa_c12 <- pcoa_c12 + coord_fixed() + 
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20))
pcoa_c13 <- pcoa_c13 + coord_fixed() + 
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20))

pdf('plots/pcoa_control.pdf')
pcoa_c12 / pcoa_c13
dev.off()

pdf('plots/pcoa_control_12.pdf')
pcoa_c12 + coord_fixed() + theme(axis.text = element_text(size = 18),
                               axis.title = element_text(size = 20))
dev.off()

pdf('plots/pcoa_control_13.pdf')
pcoa_c13 + coord_fixed() + theme(axis.text = element_text(size = 18),
                                 axis.title = element_text(size = 20))
dev.off()

car_merge <- tax_glom(carbom_abund, 'Genus')

# contaminants' abundance -------------------------------------------------
mdf <- prep_mdf(car_merge, subgroup_level = 'Genus')
color_objs <- create_color_dfs(mdf, selected_groups = selected_groups_control, 
                               group_level = 'Phylum',
                               subgroup_level = 'Genus',
                               cvd = TRUE)
mdf <- color_objs$mdf
cdf <- color_objs$cdf
new_groups <- extend_group(mdf, cdf, 'Phylum', 'Genus', 'Proteobacteria', 
                           existing_palette = 'micro_purple',
                           new_palette = 'micro_cvd_orange', n_add = 5)

order <- c('21243','21244','21245','21213','21214','21215','21228','21229',
           '21230','21198','21199','21200','21258','21259','21260','21273',
           '21274','21275','21303','21304','21305','21288','21289','21290',
           '21318','21319','21320','21321','21322','21323','21240','21241',
           '21242','21210','21211','21212','21225','21226','21227','21195',
           '21196','21197','21255','21256','21257','21270','21271','21272',
           '21300','21301','21302','21285','21286','21287','21315','21316',
           '21317')

reordered <- reorder_samples_by(new_groups$mdf, new_groups$cdf,
                                sample_variable = 'Sample',
                                sample_ordering = order)

control_bars <- plot_microshades(reordered$mdf, reordered$cdf) + 
  scale_y_continuous(labels = scales::percent, expand = expansion(0)) + 
  facet_wrap(~ Biome, scales = "free_x", ncol = 1) + 
  theme_bars +
  theme(legend.position = "none")  +
  theme(axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 16),
        
        panel.spacing = unit(1, 'lines'),
        strip.text = element_text(size = 20, face = "bold"),
        strip.background = element_blank(),
        plot.margin = margin(6, 20, 6, 6)
        )


pdf("plots/control_bars.pdf")
control_bars
dev.off()


legend <- custom_legend(new_groups$mdf, new_groups$cdf, legend_key_size = 0.8,
                        legend_text_size = 12, legend_orientation = 'horizontal',
                        group_level = 'Phylum', subgroup_level = 'Genus')
pdf("plots/control_bars_legend.pdf", width = 12)
legend
dev.off()


# final plots after decontam ----------------------------------------------

cur_samples <- env_s

order <- c('21237', '21238', '21239', '21207', '21208', '21209', '21222', '21223', '21224',
           '21192', '21193', '21194', '21252', '21253', '21254', '21267', '21268', '21269',
           '21297', '21298', '21299', '21282', '21283', '21284', '21312', '21313', '21314', # w
           '21246', '21247', '21248', '21216', '21217', '21218', '21231', '21232', '21233',
           '21201', '21202', '21203', '21261', '21262', '21263', '21276', '21277', '21278',
           '21306', '21307', '21308', '21291', '21292', '21293', '21324', '21325', '21326', # s
           '21249', '21250', '21251', '21219', '21220', '21221', '21234', '21235', '21236',
           '21204', '21205', '21206', '21264', '21265', '21266', '21279', '21280', '21281',
           '21309', '21310', '21311', '21294', '21295', '21296', '21327', '21328', '21329', # g
           '25759', '25760', '25761', '25747', '25748', '25749', '25753', '25754', '25755',
           '25741', '25742', '25743', '25765', '25766', '25767', '25771', '25772', '25773',
           '25783', '25784', '25785', '25777', '25778', '25779', '25789', '25790', '25791',
           '25795', '25796', '25797' # f
)

group_breaks = seq(3.5, 30.5, 3)


otu_mat <- load_otu_tax(otu_path = decontam_otu_path)$otu
OTU <- otu_table(otu_mat, taxa_are_rows = T)
abund <- get_relabund(OTU, TAX, cur_samples)


# microshades abundance plots ---------------------------------------------

mdf <- prep_mdf(abund, subgroup_level = 'Order')
# Create abundance dataframe & color mapping dataframe
color_objs <- create_color_dfs(mdf, selected_groups = 
                                 c("Verrucomicrobiota", 
                                   "Cyanobacteria", 
                                   "Actinobacteriota",
                                   "Firmicutes",
                                   "Proteobacteria"
                                   ), 
                               group_level = 'Phylum',
                               subgroup_level = 'Order',
                               cvd = TRUE)

mdf <- color_objs$mdf
cdf <- color_objs$cdf
new_groups <- extend_group(mdf, cdf, 'Phylum', 'Order', 'Other', 
                           existing_palette = 'micro_cvd_gray',
                           new_palette = 'micro_purple', n_add = 5)

reordered <- reorder_samples_by(new_groups$mdf, new_groups$cdf,
                                sample_variable = 'Sample',
                                sample_ordering = order)

legend <- custom_legend(new_groups$mdf, new_groups$cdf, legend_key_size = 0.8,
                        legend_text_size = 12, legend_orientation = 'horizontal',
                        group_level = 'Phylum', subgroup_level = 'Order')
center_legend <- plot_grid(NULL, legend, nrow=2, rel_widths = c(0.04, 0.96))
plot_diff <- plot_microshades(new_groups$mdf, new_groups$cdf) + 
  scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme_bars + 
  theme(legend.position = "none",
        panel.spacing = unit(0, 'lines'),
        strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_blank()) +
  facet_wrap(~Sample.type, scales = "free_x", ncol = 1) +
  theme(axis.text.x = element_text(size= 8, angle = 65)) +
  theme(plot.margin = margin(6,20,6,6)) +
  ylab('Relative abundance, %') + 
  xlab('Kit') + 
  # scale_x_discrete("", breaks = as.character(breaks), labels = labels) +
  geom_vline(xintercept = group_breaks)

plot_diff
pdf("plots/microshades.pdf")
plot_diff
dev.off()

legend
pdf("plots/microshades_legend.pdf", width = 10, height = 10)
legend
dev.off()

plot_grid(plot_diff, legend, nrow = 1, rel_widths = c(1, .45),
          rel_heights = c(1, .4)) 


# PCoA, Shannon index and rarefaction after decontam -----------------------
cur_samples <- water_s
carbom_abund <- get_relabund(OTU, TAX, cur_samples)
ps <- phyloseq(OTU, TAX, cur_samples)
tab <- as.matrix(as.data.frame(otu_table(ps)))
rare_w_d <- rarecurve(t(tab), step=500, cex=0.5, col = 'steelblue4', label = T)
order <- c('MagBac', 'MagMic', 'MagSoi', 'MagSto', 'SilMet', 'SilSoi', 'SkySto',
           'SkySoi', 'B&T')
# _w_d == water_decontam
# _12 and _13 denote principal components
pcoa_w_d_12 <- pcoa(carbom_abund, alpha = 1, color = 'Group') + 
  coord_fixed(ratio = 2) +
  theme(
    aspect.ratio = 1,
    legend.position = 'none',
    axis.title = element_text(size = 20),
    axis.title.x = element_blank()
    ) + 
  scale_color_igv(limits = order)

pcoa_w_d_13 <- pcoa(carbom_abund, nmds = c(1, 3), alpha = 1, color = 'Group') + 
  coord_fixed(ylim = c(-0.65, 0.4)) +
  theme(
    aspect.ratio = 1,
    # legend.position = 'none',
    axis.title = element_text(size = 20)
    ) + 
  scale_color_igv(limits = order)

pdf('plots/pcoa_water_decontam.pdf')
pcoa_w_d_12 / pcoa_w_d_13
dev.off()

cur_samples <- soil_s
carbom_abund <- get_relabund(OTU, TAX, cur_samples)
ps <- phyloseq(OTU, TAX, cur_samples)
tab <- as.matrix(as.data.frame(otu_table(ps)))
rare_s_d <- rarecurve(t(tab), step=500, cex=0.5, col = 'steelblue4', label = T)
order <- c('MagBac', 'MagMic', 'MagSoi', 'MagSto', 'SilMet', 'SilSoi', 'SkySto',
           'SkySoi', 'PowSoi')
# _s_d == _soil_decontam
pcoa_s_d_12 <- pcoa(carbom_abund, alpha = 1, color = 'Group') + 
  coord_fixed(ylim = c(-0.3, 0.5)) +
  theme(
    aspect.ratio = 1,
    legend.position = 'none',
    axis.title = element_text(size = 20),
    axis.title.x = element_blank()
    ) +
  scale_color_igv(limits = order)

pcoa_s_d_13 <- pcoa(carbom_abund, nmds = c(1, 3), alpha = 1, color = 'Group') + 
  coord_fixed(ylim = c(-0.3, 0.5)) +
  theme(
    aspect.ratio = 1,
    # legend.position = 'none',
    axis.title = element_text(size = 20)) +
  scale_color_igv(limits = order)

pdf('plots/pcoa_soil_decontam.pdf')
pcoa_s_d_12 / pcoa_s_d_13
dev.off()

cur_samples <- animals_s
carbom_abund <- get_relabund(OTU, TAX, cur_samples)
ps <- phyloseq(OTU, TAX, cur_samples)
tab <- as.matrix(as.data.frame(otu_table(ps)))
rare_a_d <- rarecurve(t(tab), step=500, cex=0.5, col = 'steelblue4', label = T)
order <- c('MagBac', 'MagMic', 'MagSoi', 'MagSto', 'SilMet', 'SilSoi', 'SkySto',
           'SkySoi', 'PowSoi')
# _a_d == _animals_decontam

pcoa_a_d_12 <- pcoa(carbom_abund, alpha = 1, color = 'Group') + 
  coord_fixed() +
  theme(
    aspect.ratio = 1,
    legend.position = 'none',
    axis.title = element_text(size = 20),
    axis.title.x = element_blank()
    ) + 
  scale_color_igv(limits = order)

pcoa_a_d_13 <- pcoa(carbom_abund, nmds = c(1, 3), alpha = 1, color = 'Group') + 
  coord_fixed(ylim = c(-0.3, 0.65)) +
  theme(
    aspect.ratio = 1,
    # legend.position = 'none',
    axis.title = element_text(size = 20)
    ) +
  scale_color_igv(limits = order)

pdf('plots/pcoa_animals_decontam.pdf')
pcoa_a_d_12 / pcoa_a_d_13
dev.off()

cur_samples <- feces_s
carbom_abund <- get_relabund(OTU, TAX, cur_samples)
ps <- phyloseq(OTU, TAX, cur_samples)
tab <- as.matrix(as.data.frame(otu_table(ps)))
rare_a_d <- rarecurve(t(tab), step=500, cex=0.5, col = 'steelblue4', label = T)
order <- c('MagBac', 'MagMic', 'MagSoi', 'MagSto', 'SilMet', 'SilSoi', 'SkySto',
           'SkySoi', 'PowSoi', 'PowFec')
# _a_d == _animals_decontam
pcoa_f_d_12 <- pcoa(carbom_abund, alpha = 1, color = 'Group') + 
  coord_fixed() +
  theme(
    aspect.ratio = 1,
    legend.position = 'none', 
    axis.title = element_text(size = 20),
    axis.title.x = element_blank()
    ) + 
  scale_color_igv(limits = order)

pcoa_f_d_13 <- pcoa(carbom_abund, nmds = c(1, 3), alpha = 1, color = 'Group') + 
  coord_fixed(ylim = c(-0.4, 0.4)) +
  theme(
    aspect.ratio = 1,
    # legend.position = 'none',
    axis.title = element_text(size = 20)
    ) + 
  scale_color_igv(limits = order)

pdf('plots/pcoa_feces_decontam.pdf')
pcoa_f_d_12 / pcoa_f_d_13
dev.off()



# barplots of universal OTUs ----------------------------------------------

universal.otu.df <- read.csv2(
  file.path(
    workdir, 'results/reproducibility/Universal_OTUs_by_sample_type.csv'),
  sep = ',')

cur_samples <- water_s
universal.otu <- universal.otu.df$water

uni_ps <- phyloseq(OTU, TAX, cur_samples)
uni_ps <- prune_taxa(universal.otu, uni_ps)
uni_ps <- merge_samples(uni_ps, 'Sample.type')
uni_ps <- transform_sample_counts(uni_ps, percentFraction)
uni_ps.abund <- filter_taxa(uni_ps, function(x) sum(x > 1) > 0, T)
uni_ps.abund <- transform_sample_counts(uni_ps.abund, percentFraction)

bars_order_w = plot_bar(uni_ps.abund, fill = "Order") + 
  geom_bar(aes(color = Family, fill = Order), 
           stat = "identity", position = "stack", color = NA) + 
  labs(x = '', y = '') + theme_classic() + 
  theme(legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8), legend.position = "right",
        axis.text = element_text(size = 10), axis.title = element_blank(),
        axis.text.x = element_text(angle = 65, hjust = 1), 
        axis.text.y = element_text(size = 8)
  ) +
  theme(
    strip.text = element_text(size = 12, face = 'bold'),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 8),
    strip.placement = "outside",
    panel.spacing.y = unit(2, "lines")) +
  scale_fill_manual(values = order_colors) + 
  guides(fill = guide_legend(ncol = 1))

### soil
universal.otu <- universal.otu.df$soil
cur_samples <- soil_s
uni_ps <- phyloseq(OTU, TAX, cur_samples)
uni_ps <- prune_taxa(universal.otu, uni_ps)
uni_ps <- merge_samples(uni_ps, 'Sample.type')
uni_ps <- transform_sample_counts(uni_ps, percentFraction)
uni_ps.abund <- filter_taxa(uni_ps, function(x) sum(x > 1) > 0, T)
uni_ps.abund <- transform_sample_counts(uni_ps.abund, percentFraction)

bars_order_s = plot_bar(uni_ps.abund, fill = "Order") + 
  geom_bar(aes(color = Family, fill = Order), 
           stat = "identity", position = "stack", color = NA) + 
  labs(x = '', y = '') + theme_classic() + 
  theme(legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8), legend.position = "right",
        axis.text = element_text(size = 10), axis.title = element_blank(),
        axis.text.x = element_text(angle = 65, hjust = 1), 
        axis.text.y = element_text(size = 8)
  ) +
  theme(
    strip.text = element_text(size = 12, face = 'bold'),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 8),
    strip.placement = "outside",
    panel.spacing.y = unit(2, "lines")) +
  scale_fill_manual(values = order_colors) + 
  guides(fill = guide_legend(ncol = 1))


### gut flora
universal.otu <- universal.otu.df$organism
cur_samples <- animals_s
uni_ps <- phyloseq(OTU, TAX, cur_samples)
uni_ps <- prune_taxa(universal.otu, uni_ps)
uni_ps <- merge_samples(uni_ps, 'Sample.type')
uni_ps <- transform_sample_counts(uni_ps, percentFraction)
uni_ps.abund <- filter_taxa(uni_ps, function(x) sum(x > 1) > 0, T)
uni_ps.abund <- transform_sample_counts(uni_ps.abund, percentFraction)

bars_order_a = plot_bar(uni_ps.abund, fill = "Order") + 
  geom_bar(aes(color = Family, fill = Order), 
           stat = "identity", position = "stack", color = NA) + 
  labs(x = '', y = '') + theme_classic() + 
  theme(legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8), legend.position = "right",
        axis.text = element_text(size = 10), axis.title = element_blank(),
        axis.text.x = element_text(angle = 65, hjust = 1), 
        axis.text.y = element_text(size = 8)
  ) +
  theme(
    strip.text = element_text(size = 12, face = 'bold'),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 8),
    strip.placement = "outside",
    panel.spacing.y = unit(2, "lines")) +
  scale_fill_manual(values = order_colors) + 
  guides(fill = guide_legend(ncol = 1))

### feces
universal.otu <- universal.otu.df$feces
cur_samples <- feces_s
uni_ps <- phyloseq(OTU, TAX, cur_samples)
uni_ps <- prune_taxa(universal.otu, uni_ps)
uni_ps <- merge_samples(uni_ps, 'Sample.type')
uni_ps <- transform_sample_counts(uni_ps, percentFraction)
uni_ps.abund <- filter_taxa(uni_ps, function(x) sum(x > 1) > 0, T)
uni_ps.abund <- transform_sample_counts(uni_ps.abund, percentFraction)

bars_order_f = plot_bar(uni_ps.abund, fill = "Order") + 
  geom_bar(aes(color = Family, fill = Order), 
           stat = "identity", position = "stack", color = NA) + 
  labs(x = '', y = '') + theme_classic() + 
  theme(legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8), legend.position = "right",
        axis.text = element_text(size = 10), axis.title = element_blank(),
        axis.text.x = element_text(angle = 65, hjust = 1), 
        axis.text.y = element_text(size = 8)
  ) +
  theme(
    strip.text = element_text(size = 12, face = 'bold'),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 8),
    strip.placement = "outside",
    panel.spacing.y = unit(2, "lines")) +
  scale_fill_manual(values = order_colors) + 
  guides(fill = guide_legend(ncol = 1))

pdf('plots/universal_otu_barplots.pdf', width = 12)
bars_order_w | bars_order_s | bars_order_a | bars_order_f
dev.off()

cur_samples <- env_s
universal.otu <- unlist(universal.otu.df)
uni_ps <- phyloseq(OTU, TAX, cur_samples)
uni_ps <- prune_taxa(universal.otu, uni_ps)
uni_ps <- merge_samples(uni_ps, 'Sample.type')
uni_ps <- transform_sample_counts(uni_ps, percentFraction)
uni_ps.abund <- filter_taxa(uni_ps, function(x) sum(x > 1) > 0, T)
uni_ps.abund <- transform_sample_counts(uni_ps.abund, percentFraction)

mdf <- prep_mdf(uni_ps.abund, subgroup_level = 'Order')
# Create abundance dataframe & color mapping dataframe
color_objs <- create_color_dfs(mdf, selected_groups = 
                                 c("Verrucomicrobiota", 
                                   "Cyanobacteria", 
                                   "Actinobacteriota",
                                   "Firmicutes",
                                   "Proteobacteria"
                                 ), 
                               group_level = 'Phylum',
                               subgroup_level = 'Order',
                               cvd = TRUE)

mdf <- color_objs$mdf
cdf <- color_objs$cdf
new_groups <- extend_group(mdf, cdf, 'Phylum', 'Order', 'Other', 
                           existing_palette = 'micro_cvd_green',
                           new_palette = 'micro_purple', n_add = 5)

plot_diff <- plot_microshades(new_groups$mdf, new_groups$cdf) + 
  scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme_bars + 
  theme(legend.position = "none",
        panel.spacing = unit(0, 'lines'),
        strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_blank()) +
  facet_wrap(~Sample.type, scales = "free_x", ncol = 4) +
  theme(axis.text.x = element_text(size= 8, angle = 65)) +
  theme(plot.margin = margin(6,20,6,6)) +
  ylab('Relative abundance, %') + 
  xlab('Kit')

legend <- custom_legend(new_groups$mdf, new_groups$cdf, legend_key_size = 0.8,
                        legend_text_size = 12, legend_orientation = 'horizontal',
                        group_level = 'Phylum', subgroup_level = 'Order')
pdf('plots/microshades_universal.pdf')
plot_diff
dev.off()

pdf('plots/microshades_universal_legend.pdf')
legend
dev.off()
