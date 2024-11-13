library(phyloseq)
library(dplyr)
library(ggplot2)
library(vegan)
library(ape)
library(ggpubr)
library(scales)
library(ggradar)
library(ggbeeswarm)
library(patchwork)
library(ggsci)
library(lemon)

workdir <- "/Users/asobolev/skoltech/articles/kits"
dna_stats <- read.csv2(file.path(workdir, "wetlab_metrics_upd.csv"),
                       header = TRUE, sep = ",", check.names = FALSE)

dna_stats$Kit <- as.factor(dna_stats$Kit)
levels(dna_stats$Kit) <- list(
  "MagBac" = "Magen Bacterial",
  "MagMic" = "Magen Microbiome",
  "MagSoi" = "Magen Soil",
  "MagSto" = "Magen Stool",
  "SilMet" = "Sileks Metagenomic",
  "SilSoi" = "Sileks Soil",
  "SkyFec" = "Skygen Fecal",
  "SkySoi" = "Skygen Soil",
  "B&T" = "Blood&Tissue",
  "PowSoi" = "PowerSoil",
  "PowFec" = "PowerFecal"
)

numcol <- c("DNA, ng/mcl", "DNA, ng/mcl (log scale)",
            "DIN", "260/280", "260/230", "18S-16S")

dna_stats[numcol] <- lapply(dna_stats[numcol],
                            function(x) (as.numeric(as.vector(x))))

dna_stats$Sample <- as.factor(dna_stats$Sample)
dna_stats$Sample <- factor(dna_stats$Sample, levels = c("Water", "Sediment",
                                                        "Gut flora", "Feces"))

### facet_wrap custom function
kits_order <- c("MagBac", "MagMic", "MagSoi", "MagSto", "SilMet", "SilSoi",
                "SkySto", "SkySoi", "B&T", "PowSoi", "PowFec")

facet_kit <- function(metric) {
  p <- ggplot(dna_stats, aes(x = Kit, y = {{metric}}, color = Kit)) +
    geom_point(size = 2, show.legend = FALSE) + facet_wrap(~ Sample, ncol = 4) +
    scale_x_discrete() + scale_color_igv() +
    theme_classic(base_size = 10) +
    theme(axis.text.y = element_text(size = 9, color = "black"),
          axis.text.x = element_text(size = 9, angle = 90, hjust = 1,
                                     vjust = 0.5, color = "black"),
          strip.background = element_blank(),
          strip.text = element_text(size = 12, face = "bold"),
          panel.spacing = unit(1, "lines"))
  return(p)
}

facet_kit(`DIN`)

p_dna <- facet_kit(metric = `DNA, ng/mcl`)
p_logdna <- facet_kit(metric = `DNA, ng/mcl (log scale)`)
p_din <- facet_kit(metric = `DIN`) + geom_hline(yintercept = 3, col = "red")
p_260280 <- facet_kit(metric = `260/280`) + geom_hline(yintercept = 1.9,
                                                       col = "red")
p_260230 <- facet_kit(metric = `260/230`) + geom_hline(yintercept = 1.9,
                                                       col = "red")
p_18s16s <- facet_kit(metric = `18S-16S`) + ylab("Ct(18S) - Ct(16S)") +
  geom_hline(yintercept = 0, col = "red")

pdf("plots/logdna_din.pdf")
p_logdna / p_din + plot_annotation(tag_levels = "A")
dev.off()

pdf("plots/absorbance.pdf")
p_260280 / p_260230 + plot_annotation(tag_levels = "A")
dev.off()

vs <- ggplot(dna_stats, aes(x = DIN, y = `DNA, ng/mcl`, color = Kit)) +
  facet_rep_wrap(~ Sample, scales = "free", repeat.tick.labels = "all") +
  geom_point(size = 3, alpha = 0.7) + scale_color_igv() +
  xlab("DIN") + ylab("DNA concentration, ng/mcl") +
  scale_shape_manual(values = c(15, 16, 17, 18, 4, 10, 20, 9, 13, 22, 25)) +
  scale_x_continuous(breaks = seq(0, 7, 1), limits = c(0, 7)) +
  theme_bw(base_size = 10) +
  theme(aspect.ratio = 1,
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, hjust = 1,
                                   vjust = 0.5, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        panel.spacing = unit(1, "lines"))

pdf("plots/din_vs_yield.pdf")
vs
dev.off()

pdf("plots/18s-16s.pdf")
p_18s16s / p_logdna
dev.off()
