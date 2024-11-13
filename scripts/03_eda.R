library(dplyr)
library(reshape2)
library(vegan)
library(cowplot)

workdir <- "/Users/asobolev/skoltech/articles/kits/"
source(file.path(workdir, "scripts/01_preload.R"))

### merge stats
merge_stats <- read.csv2(file.path(workdir, "results/track_stats_DADA2.tsv"),
                         header = TRUE, check.names = FALSE, sep = " ",
                         row.names = NULL)
merge_stats$Sample_ID <- lapply(merge_stats$row.names,
                                function(x) (strsplit(x, "_")[[1]][1]))
merge_stats <- merge_stats %>% select(-c(passed_frac, row.names))
merge_stats <- mutate_all(merge_stats, function(x) as.numeric(as.character(x)))
merge_stats$Sample_ID <- as.integer(as.vector(merge_stats$Sample_ID))

env <- subset(samples_all, samples_all$`Sample type` %in% c("MQ", "Kitome"))
cols <- c(colnames(merge_stats), "Group", "Sample type")
stats_df <- right_join(merge_stats, samples_all, by = c("Sample_ID")) %>%
  select(all_of(cols))
stats_df$Sample_ID <- as.character(stats_df$Sample_ID)
stats_df$isneg <- ifelse(stats_df$`Sample type` %in% c("MQ", "Kitome"),
                         "Control",
                         "Environment")

stat_mean <- melt(stats_df, id.vars = c(7, 8, 9, 10)) %>%
  group_by(Group, isneg, variable) %>%
  summarise(mean = mean(value), median = median(value), sd = sd(value))

stat_meandf <- subset(stat_mean,
                      stat_mean$variable == "nonchim" &
                        stat_mean$Group != "PowFec") %>%
  group_by(Group) %>%
  mutate(
    ratio_median = median[isneg == "Environment"] / median[isneg == "Control"],
    ratio_mean = mean[isneg == "Environment"] / mean[isneg == "Control"]
  )

kits_order <- c("MagBac",
                "MagMic",
                "MagSoi",
                "MagSto",
                "SilMet",
                "SilSoi",
                "SkySto",
                "SkySoi",
                "B&T",
                "PowSoi",
                "PowFec")
stats_df$Group <- factor(stats_df$Group, levels = kits_order)

p_stat <- ggplot(stats_df, aes(x = Group, y = nonchim, color = isneg,
                               shape = `isneg`)) +
  facet_wrap(~ `Group`, scales = "free_x") +
  geom_boxplot(outlier.shape = NA, position = "dodge") +
  geom_jitter(position = position_dodge(width = 0.75)) +
  scale_color_igv() + ylab("Read counts") +
  theme_bw() + theme_bars + theme(axis.title = element_text(size = 12))

pdf("plots/nonchim_readcounts.pdf")
p_stat
dev.off()



input_stat <- ggplot(stats_df, aes(x = Sample_ID, y = input, fill = Group)) +
  geom_bar(stat = "identity") +
  facet_wrap(~`Sample type`, scales = "free_x") +
  scale_fill_igv() +
  theme_classic() + theme_bars
ggsave(file.path(workdir, "plots/prep_input_stats.png"),
       input_stat, width = 1309 / 90, height = 818 / 90, dpi = 300)

mean_p <- ggplot(stat_mean, aes(x = variable, y = mean, fill = Group)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean, ymax = mean + sd)) +
  facet_wrap(~ Group, nrow = 2) +
  scale_fill_igv() +
  theme_light() +
  theme(
    strip.text = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(angle = 90, size = 10),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 8),
    panel.spacing.y = unit(2, "lines")
  ) +
  labs(x = "", y = "")
ggsave(file.path(workdir, "plots/prep_merge_stats.png"),
       mean_p, width = 1309 / 90, height = 818 / 90, dpi = 300)

### abundance
abund <- get_relabund(otu, tax, samples)
pcoa_prep <- pcoa(abund, color = "Sample.type")
ggsave(file.path(workdir, "plots/pcoa_prep.png"),
       pcoa_prep, width = 1309 / 90, height = 818 / 90, dpi = 300)

abund <- get_relabund(otu, tax, env_s)

mdf <- prep_mdf(abund, subgroup_level = "Order", remove_na = TRUE)
# Create abundance dataframe & color mapping dataframe
color_objs_updated <- create_color_dfs(mdf,
                                       selected_groups =
                                         c("Verrucomicrobiota",
                                           "Cyanobacteria",
                                           "Actinobacteriota",
                                           "Firmicutes",
                                           "Proteobacteria"
                                         ),
                                       group_level = "Phylum",
                                       subgroup_level = "Order",
                                       cvd = TRUE)

mdf_updated <- color_objs_updated$mdf
cdf_updated <- color_objs_updated$cdf
new_groups <- extend_group(mdf_updated, cdf_updated,
                           "Phylum", "Order", "Other",
                           existing_palette = "micro_cvd_gray",
                           new_palette = "micro_purple", n_add = 5)
legend_new <- custom_legend(new_groups$mdf, new_groups$cdf,
                            group_level = "Phylum",
                            subgroup_level = "Order")
center_legend <- plot_grid(NULL, legend_new, nrow = 1,
                           rel_widths = c(0.04, 0.96))
plot_diff <- plot_microshades(new_groups$mdf, new_groups$cdf) +
  scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 6, angle = 90),
        panel.spacing = unit(0, "lines"),
        strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_blank()) +
  facet_wrap(~Sample.type, scales = "free_x") +
  theme(axis.text.x = element_text(size = 8)) +
  theme(plot.margin = margin(6, 20, 6, 6)) +
  ylab("Relative abundance, %") + xlab("Kit")
plot_grid(plot_diff, legend_new, rel_widths = c(1, .25))

bars_order <- plot_bar(abund, fill = "Order") +
  geom_bar(aes(fill = Order), color = "black", linewidth = 0.05,
           stat = "identity", position = "stack") +
  theme_classic() + theme_bars + labs(x = "", y = "") +
  facet_wrap(Sample.type ~ ., nrow = 4, scales = "free_x",
             strip.position = "left") +
  scale_fill_manual(values = order_colors) +
  guides(fill = guide_legend(ncol = 2))

### rarefaction

### control samples
control_s <- rbind(kitome_s, splashome_s)
cur_samples <- control_s
sdata <- sample_data(cur_samples)
ps <- phyloseq(otu, tax, sdata)

tab <- as.matrix(as.data.frame(otu_table(ps)))

rare_c <- rarecurve(t(tab), step = 500, cex = 0.5, col = "steelblue4",
                    label = TRUE)

cur_samples <- soil_s
sdata <- sample_data(cur_samples)
cur_abund <- get_relabund(otu, tax, sdata)
ps <- phyloseq(otu, tax, sdata)
tab <- as.matrix(as.data.frame(otu_table(ps)))
rare_c <- rarecurve(t(tab), step = 500, cex = 0.5, col = "steelblue4",
                    label = TRUE)
pcoa_soil <- pcoa(cur_abund, color = "Group")

cur_samples <- water_s
sdata <- sample_data(cur_samples)
cur_abund <- get_relabund(otu, tax, sdata)
ps <- phyloseq(otu, tax, sdata)
tab <- as.matrix(as.data.frame(otu_table(ps)))
rare_c <- rarecurve(t(tab), step = 500, cex = 0.5, col = "steelblue4",
                    label = TRUE)
pcoa_water <- pcoa(cur_abund, color = "Group")

cur_samples <- feces_s
sdata <- sample_data(cur_samples)
cur_abund <- get_relabund(otu, tax, sdata)
ps <- phyloseq(otu, tax, sdata)
tab <- as.matrix(as.data.frame(otu_table(ps)))
rare_c <- rarecurve(t(tab), step = 500, cex = 0.5, col = "steelblue4",
                    label = TRUE)
pcoa_feces <- pcoa(cur_abund, color = "Group")

cur_samples <- animals_s
sdata <- sample_data(cur_samples)
cur_abund <- get_relabund(otu, tax, sdata)
ps <- phyloseq(otu, tax, sdata)
tab <- as.matrix(as.data.frame(otu_table(ps)))
rare_c <- rarecurve(t(tab), step = 500, cex = 0.5, col = "steelblue4",
                    label = TRUE)
pcoa_animals <- pcoa(cur_abund, color = "Group")
