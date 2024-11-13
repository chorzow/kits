library(phyloseq)
library(decontam)
library(jsonlite)

workdir <- "/Users/asobolev/skoltech/articles/kits"

source(file.path(workdir, "scripts/01_preload.R"))

# for animal, soil and feces, kitome samples are used as negative control
asf <- subset(samples,
              samples$`Sample type` %in%
                c("Gut flora", "Sediment", "Feces", "Kitome"))
asf$is.neg <- asf$Sample.type == "Kitome"

# for water, miliq samples are used as negative control
water_all <- subset(samples, samples$`Sample type`  %in% c("Water", "MQ"))
water_all$is.neg <- water_all$Sample.type == "MQ"


samples_df <- rbind(asf, water_all)
samples <- sample_data(samples_df) # metadata for phyloseq object
rownames(samples) <- samples$Sample_ID

methods <- c("either")
threshold_array <- seq(0.1, 0.5, 0.1)

experiment_list <- list()
experiment_list[[1]] <- asf$Sample_ID
experiment_list[[2]] <- water_all$Sample_ID

output_dir <- file.path(workdir, "results/decontam")

json_out <- data.frame()


for (method in methods) {
  print(paste("Start method:", method))
  for (threshold in threshold_array){
    print(paste("Use threshold:", threshold))

    json_out <- data.frame()  # empty dataframe for contaminants

    for (experiment in experiment_list) {
      print(experiment)

      current_set <- subset(samples, samples$Sample_ID %in% experiment)

      current_set$DNA.concentration..ng.mkl <-
        as.numeric(as.vector(current_set$DNA.concentration..ng.mkl))

      sample_set <- subset(current_set, current_set$is.neg != TRUE)$Sample_ID


      ## make phyloseq object
      current_samples <- current_set
      rownames(current_samples) <- current_samples$Sample_ID

      ps <- phyloseq(otu, tax, current_samples)  # ps == phyloseq object


      print("decontaminating...")
      ## run decontam
      contamdf_either <- isContaminant(ps, method = method, neg = "is.neg",
                                       conc = "DNA.concentration..ng.mkl",
                                       threshold = threshold, normalize = TRUE)

      print(table(contamdf_either$contaminant))

      json_columns <- data.frame(sample = sample_set)
      json_columns$contams <- list(c(which(contamdf_either$contaminant)))
      json_out <- rbind(json_out, json_columns)
    }

    json_string <- toJSON(json_out, auto_unbox = TRUE, pretty = TRUE)
    print("Save JSON...")
    write(json_string,
          paste0(output_dir,
                 "/contamin_list_",
                 method,
                 "_",
                 as.character(threshold),
                 "_grouped_all_vs_all_2.json"))

    ### filter the contaminant OTUs
    all_samples <- sample_data(samples_df)
    rownames(all_samples) <- all_samples$Sample_ID
    ps <- phyloseq(otu, tax, all_samples)

    print(sample_names(ps@otu_table))

    otu_df <- as.data.frame(otu_table(ps))

    otu_df[json_out[1, 2][[1]],
           json_out$sample %in% colnames(otu_df)] <- 0  # set contam. freq. to 0

    ### save otu table
    write.table(otu_df,
                paste0(output_dir,
                       "/OTUs_decontam_",
                       method,
                       "_",
                       as.character(threshold),
                       "_grouped_all_vs_all_2.tsv"),
                sep = "\t", row.names = TRUE, col.names = NA)
    print("Saved filtered OTU table...")

    ### save contaminants table

    carbom_pos <- prune_samples(sample_data(ps)$is.neg == FALSE, ps)

    otu_sums <- sample_sums(carbom_pos)
    carbom_contam <- prune_taxa(contamdf_either$contaminant == TRUE, carbom_pos)

    otu_contam <- sample_sums(carbom_contam)
    otu_res <- otu_contam / otu_sums * 100
    write.table(otu_res,
                file.path(workdir,
                          paste0("results/decontam/contam_levels_",
                                 method, "_", threshold,
                                 "_relative.tsv")),
                sep = "\t", row.names = TRUE, col.names = TRUE)

    otu_contam <- data.frame(otu_table(carbom_contam), check.names = FALSE)
    new_otu <- otu_contam

    for (i in seq(1, length(otu_sums))) {
      for (j in seq(1, length(otu_contam[, 1]))) {
        new_otu[j, i] <- otu_contam[j, i] / as.numeric(otu_sums[i]) * 100
      }
    }
    write.table(new_otu,
                file.path(workdir,
                          paste0("results/decontam/otu_contam_",
                                 threshold, "_relative.tsv")),
                sep = "\t", row.names = TRUE, col.names = TRUE)

    otu_df <- as.data.frame(otu_table(ps))
    contam_otu <- data.frame(matrix(0,
                                    nrow = nrow(otu_table(ps)),
                                    ncol = ncol(otu_table(ps))))
    colnames(contam_otu) <- colnames(otu_table(ps))
    rownames(contam_otu) <- rownames(otu_table(ps))

    set1 <- as.vector(json_out[1, 2][[1]])
    set2 <- seq_len(nrow(otu_df))
    diff <- setdiff(set2, set1)

    contam_otu[set1, ] <- otu_df[set1, ]  # set non-contaminant freq. to 0

    write.table(contam_otu,
                paste0(output_dir,
                       "/OTUs_contamin_",
                       method,
                       "_",
                       as.character(threshold),
                       "_grouped_all_vs_all.tsv"),
                sep = "\t", row.names = TRUE, col.names = TRUE)
    print("Saved contaminated OTU table...")
  }
}
