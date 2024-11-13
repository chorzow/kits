library(dplyr)

workdir <- "/Users/asobolev/skoltech/articles/kits/"

source(file.path(workdir, "scripts/01_preload.R"))

method <- "either"  # used for decontam. TODO: add to config
threshold <- 0.2  # used for decontam. TODO: add to config

### read otu table with contaminants' relative abundance

otu_mat_contam <- read.csv2(file.path(workdir,
                                      paste0("results/decontam/otu_contam_",
                                             threshold, "_relative.tsv")),
                            header = TRUE, check.names = FALSE, sep = "\t")

# rename otu_table
sample_id <- colnames(otu_mat_contam)
id_names <- c()
for (i in sample_id) {
  id_names <- c(id_names, substr(i, 1, 5))
}
colnames(otu_mat_contam) <- id_names
otu_mat_contam <- mutate_all(otu_mat_contam,
                             function(x) as.numeric(as.character(x)))

otu_mat_contam <- as.matrix(otu_mat_contam)

# taxonomy table
tax_mat <- read.csv2(file.path(workdir, "results/all_phylogeny.tsv"),
                    header = TRUE, check.names = FALSE, sep = "\t")
tax_mat[1] <- paste0("OTU_", seq(nrow(tax_mat)))


# filter metadata with IDs
sample_id <- colnames(otu_mat_contam)
samples_df <- data.frame()
for (i in sample_id){
  print(i)
  row <- samples_all[samples_all$Sample_ID == i, ]
  samples_df <- rbind(samples_df, row)
}

rm(row)

OTU_name <- tax_mat[1]
colnames(OTU_name) <- "OTU"

rownames(tax_mat) <- OTU_name$OTU
tax_mat[1] <- NULL
rownames(otu_mat) <- OTU_name$OTU
# otu_mat[1] = NULL

# create phyloseq object
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU <- otu_table(otu_mat_contam, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat)
samples <- sample_data(env_s)
rownames(samples) <- samples$Sample_ID


#### All contaminants
ps_contam <- phyloseq(OTU, TAX, samples)

# ps_contam = transform_sample_counts(ps_contam, percentFraction)
sample_sums(ps_contam)

ps_contam_abund <- filter_taxa(ps_contam, function(x) sum(x > 1) > 0, TRUE)

### barplots
contam_bars_genus <- plot_bar(ps_contam_abund, fill = "Genus") +
  geom_bar(aes(fill = Genus),
           stat = "identity", position = "stack", color = NA) +
  labs(x = "", y = "") + theme_classic() + theme_bars + scale_fill_frontiers() +
  facet_wrap(~Sample.type, nrow = 4, scales = "free",
             strip.position = "top") +
  guides(fill = guide_legend(ncol = 1))
pdf("plots/contam_levels_bysample.pdf", width = 10, height = 10)
contam_bars_genus
dev.off()


rel_abund_df <- data.frame(sample_sums(ps_contam))
colnames(rel_abund_df) <- "Contamination level, %"
rel_abund_df <- cbind(rel_abund_df, samples)

env_df <- read.csv2(paste0(workdir, "/results/decontam/contam_levels_",
                           method, "_", threshold, "_relative.tsv"),
                    sep = "\t", check.names = FALSE)
colnames(env_df) <- "Contamination level, %"
env_df$`Contamination level, %` <- as.numeric(env_df$`Contamination level, %`)

mean_contam_df <- cbind(env_df, env_s) %>%
  group_by(Group, Sample.type) %>%
  summarise(mean = mean(`Contamination level, %`),
            sd = sd(`Contamination level, %`))

mean_contam_df_bysample <- cbind(env_df, env_s) %>%
  group_by(Sample.type) %>%
  summarise(mean = mean(`Contamination level, %`),
            sd = sd(`Contamination level, %`))

# env_df <- rel_abund_df %>% 
#   group_by(Group, Host.organism, Step) %>% 
#   summarise(mean = mean(`Contamination level, %`),
#             sd = sd(`Contamination level, %`))

# write.table(env_df, file = file.path(workdir, 'results/'))

order <- c("MagBac", "MagMic", "MagSoi", "MagSto", "SilMet", "SilSoi", "SkySto",
           "SkySoi", "B&T", "PowSoi", "PowFec")
mean_contam_df$Group <- factor(mean_contam_df$Group, levels = order)

### plot contamination level bars
mean_contam_levels <- ggplot(mean_contam_df, aes(x = Group, y = mean,
                                                 fill = Group)) +
  geom_bar(stat = "identity") + theme_classic(base_size = 10) +
  facet_wrap(~Sample.type, scales = "free_x", nrow = 2) + theme_bars +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .2) +
  labs(y = "Contamination level, %") +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.text.x = element_text(size = 20, angle = 90, hjust = 1,
                               vjust = 0.5, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    strip.text.y = element_text(size = 20)
  ) +
  scale_fill_igv()

pdf("plots/contam_levels.pdf")
mean_contam_levels
dev.off()


# shannon and contam vs shannon -------------------------------------------

# load decontaminated data
otu_mat <- load_otu_tax(otu_path = paste0("results/decontam/OTUs_decontam_",
                                          method, "_", threshold,
                                          "_grouped_all_vs_all_2.tsv"))$otu

OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
cur_samples <- water_s
# carbom_abund <- get_relabund(OTU, TAX, cur_samples)
ps <- phyloseq(OTU, TAX, cur_samples)
water_shannon_data <- estimate_richness(ps, measures = c("Chao1", "Shannon"))

cur_samples <- soil_s
# carbom_abund <- get_relabund(OTU, TAX, cur_samples)
ps <- phyloseq(OTU, TAX, cur_samples)
soil_shannon_data <- estimate_richness(ps, measures = c("Chao1", "Shannon"))

cur_samples <- animals_s
# carbom_abund <- get_relabund(OTU, TAX, cur_samples)
ps <- phyloseq(OTU, TAX, cur_samples)
animals_shannon_data <- estimate_richness(ps, measures = c("Chao1", "Shannon"))

cur_samples <- feces_s
# carbom_abund <- get_relabund(OTU, TAX, cur_samples)
ps <- phyloseq(OTU, TAX, cur_samples)
feces_shannon_data <- estimate_richness(ps, measures = c("Chao1", "Shannon"))

shannon_data <- rbind(water_shannon_data, soil_shannon_data, 
                      animals_shannon_data, feces_shannon_data)
rownames(shannon_data) <- gsub("X", "", rownames(shannon_data))

rel_abund_df <- cbind(shannon_data, env_df)
rel_abund_df <- cbind(rel_abund_df, env_s)

order <- c("MagBac", "MagMic", "MagSoi", "MagSto", "SilMet", "SilSoi", "SkySto",
           "SkySoi", "B&T", "PowSoi", "PowFec")
rel_abund_df$Group <- factor(rel_abund_df$Group, levels = order)

shannon_vs_contam <- ggplot(rel_abund_df,
                            aes(x = `Contamination level, %`, y = Shannon,
                                color = Group, shape = Group)) +
  geom_point(size = 6, alpha = 0.6) + scale_color_igv() +
  xlab("Contamination level, %") + ylab("Shannon index") +
  scale_shape_manual(values = c(15, 16, 17, 18, 4, 10, 20, 9, 13, 22, 25)) +
  theme_classic(base_size = 10)

pdf("plots/shannon_vs_contam.pdf")
shannon_vs_contam
dev.off()

shannon <- ggplot(rel_abund_df, aes(x = Group, y = Shannon,
                                    color = Group)) +
  geom_point(size = 5, show.legend = FALSE) +
  facet_wrap(~ Sample.type, scales = "free_x", nrow = 2) +
  theme_bars +
  theme_classic(base_size = 10) +
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 20, angle = 90, hjust = 1,
                                   vjust = 0.5, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 20, face = "bold"),
        panel.spacing = unit(1, "lines")) +
  labs(y = "Shannon index") +
  scale_color_igv()

pdf("plots/shannon.pdf")
shannon
dev.off()
