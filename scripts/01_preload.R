library(phyloseq)
library(Biostrings)
library(stringr)
library(grid)
library(ggsci)
library(ggplot2)

workdir <- "/Users/asobolev/skoltech/articles/kits/"

all_otu_path <- "results/all_OTU_frequency.tsv"
all_tax_path <- "results/all_phylogeny.tsv"

# metadata
samples_all <- read.csv2(file.path(workdir, "metadata_nofc.csv"),
                         header = TRUE, check.names = FALSE, sep = ",",
                         row.names = NULL)

samples_all$Group <- as.factor(samples_all$Group)
samples_all$`Sample type` <- as.factor(samples_all$`Sample type`)

levels(samples_all$`Sample type`) <- list(
  "Water" = "water",
  "MQ" = "control",
  "Kitome" = "kitome",
  "Sediment" = "soil",
  "Gut flora" = "organism",
  "Feces" = "feces"
)

load_otu_tax <- function(otu_path = all_otu_path, tax_path = all_tax_path) {
  otu_mat <- read.csv2(otu_path, header = TRUE, check.names = FALSE, sep = "\t")
  tax_mat <- read.csv2(tax_path, header = TRUE, check.names = FALSE, sep = "\t")

  sample_id <- colnames(otu_mat)
  id_names <- c()
  for (i in sample_id){
    id_names <- c(id_names, substr(i, 1, 5))
  }
  colnames(otu_mat) <- id_names

  tax_mat[1] <- paste0("OTU_", seq_len(nrow(tax_mat)))
  otu_name <- tax_mat[1]
  colnames(otu_name) <- "OTU"

  rownames(otu_mat) <- otu_name$OTU
  rownames(tax_mat) <- otu_name$OTU
  otu_mat[1] <- NULL
  tax_mat[1] <- NULL

  otu_mat <- as.matrix(otu_mat)
  tax_mat <- as.matrix(tax_mat)
  return(list(otu = otu_mat, tax = tax_mat))
}

percent_fraction <- function(x) (x / sum(x) * 100)

get_relabund <- function(otu, tax, cur_samples, abund = 1) {
  ps <- phyloseq(otu, tax, cur_samples)
  print("Transforming into percent abundance...")
  ps <- transform_sample_counts(ps, percent_fraction)
  print("Filtering taxa...")
  ps_abund <- filter_taxa(ps, function(x) sum(x > abund) > 0, TRUE)
  print("Transforming into relative percent...")
  ps_abund <- transform_sample_counts(ps_abund, percent_fraction)
  return(ps_abund)
}

### pcoa plotting
pcoa <- function(abund_ps, nmds = c(1, 2), shapes, alpha = 1, color) {
  if (missing(shapes)) {
    shapes <- NULL
  }
  set.seed(42)
  print("Ordinating...")
  ps_ord <- ordinate(abund_ps, "PCoA", "bray", k = 3)
  pcoa_plot <- plot_ordination(abund_ps, ps_ord, axes = nmds,
                               type = "samples", color = color,
                               shape = shapes) +
    geom_point(size = 5, alpha = alpha) + scale_color_igv() +
    scale_shape_manual(values = c(15, 16, 17, 18, 3, 4, 8, 9, 10, 11)) +
    theme_light() +
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12))
  return(pcoa_plot)
}

# create phyloseq object --------------------------------------------------

otu_mat <- load_otu_tax(otu_path = file.path(workdir, all_otu_path),
                        tax_path = file.path(workdir, all_tax_path))$otu
tax_mat <- load_otu_tax(otu_path = file.path(workdir, all_otu_path),
                        tax_path = file.path(workdir, all_tax_path))$tax

otu <- otu_table(otu_mat, taxa_are_rows = TRUE)
tax <- tax_table(tax_mat)
samples <- sample_data(samples_all)
rownames(samples) <- samples$Sample_ID


ps <- phyloseq(otu, tax, samples)


# transform data to relative ----------------------------------------------

ps <- transform_sample_counts(ps, percent_fraction)
ps_abund <- filter_taxa(ps, function(x) sum(x > 1) > 0, TRUE)
ps_abund <- transform_sample_counts(ps_abund, percent_fraction)

tax_table_df <- as.data.frame(tax_table(ps_abund))

unique_phyla <- unique(tax_table_df[["Phylum"]])
unique_orders <- unique(tax_table_df[["Order"]])

# Print the list of unique orders
print(unique_orders)
print(paste("Number of unique orders in the table:", length(unique_orders)))

custom_colors_100 <- c("#962c59", "#60e54e", "#7b38d6", "#a0e53f", "#3a2aa4",
                       "#e3e939", "#4e63e6", "#60b937", "#cd48d9", "#5be389",
                       "#e135ac", "#3fa147", "#702695", "#afc538", "#8f5fd9",
                       "#b5eb7c", "#31176c", "#dec441", "#2d367f", "#e0a12d",
                       "#7385e5", "#6d9c32", "#9f44a1", "#66e1ae", "#e1352b",
                       "#64e9d9", "#c93081", "#88c777", "#cd81e0", "#3a762c",
                       "#e26bba", "#e2e582", "#37144c", "#c9eaac", "#6f276b",
                       "#938920", "#4f59ad", "#d87a29", "#6797d8", "#da5425",
                       "#53cad0", "#d63c52", "#3da67d", "#e35680", "#2a5424",
                       "#dca8e1", "#1f3318", "#c5eeda", "#261332", "#aab364",
                       "#3b3462", "#ce9c54", "#2e699c", "#dd7659", "#71c7ea",
                       "#91331b", "#9ccbcd", "#7c272f", "#8ec199", "#5a213d",
                       "#d7d5b5", "#301b26", "#e6c28c", "#122630", "#e3bac4",
                       "#37301c", "#b1b8db", "#4c251a", "#5d98b3", "#8f5e28",
                       "#9674b3", "#687730", "#c5709d", "#558660", "#774c73",
                       "#4b928c", "#ac5f61", "#224c4f", "#db9484", "#3e445b",
                       "#918054", "#65628c", "#594f20", "#a68093", "#486555",
                       "#9e9f88", "#537085", "#765850", "#0d7cff", "#208600",
                       "#00c3aa", "#99ff70", "#004c8c", "#655156", "#fff640",
                       "#dfadba", "#001044", "#ff6560", "#001e25", "#ffc12b")

order_colors <- setNames(custom_colors_100[seq_along(unique_orders)],
                         unique_orders)

theme_bars <- theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10),
  legend.position = "right",
  axis.text = element_text(size = 12),
  axis.title = element_blank(),
  axis.text.x = element_text(angle = 90, hjust = 1),
  axis.text.y = element_text(size = 12),
  strip.text = element_text(size = 12, face = "bold"),
  strip.background = element_blank(),
  strip.text.x = element_text(size = 12),
  strip.placement = "outside",
  panel.spacing.y = unit(1, "lines")
)

# subset data -------------------------------------------------------------

all_data <- subset(samples, (samples$`Sample type` != ""))
env_s <- subset(samples, (!samples$`Sample type` %in% c("MQ", "Kitome")))
soil_s <- subset(samples, (samples$`Sample type` == "Sediment"))
water_s <- subset(samples, (samples$`Sample type` == "Water"))
animals_s <- subset(samples, (samples$`Sample type` == "Gut flora"))
feces_s <- subset(samples, (samples$`Sample type` == "Feces"))
kitome_s <- subset(samples, (samples$`Sample type` == "Kitome"))
splashome_s <- subset(samples, (samples$`Sample type` == "MQ"))