# Processing the Seurat object to get UMAP/clusters
prefilter_atac_processing <- function(){
  # extract gene annotations from EnsDb
  library(EnsDb.Hsapiens.v86)
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  # change to UCSC style since the data was mapped to hg19
  seqlevelsStyle(annotations) <- 'UCSC'
  Annotation(signac_object) <- annotations
  signac_object <- NucleosomeSignal(object = signac_object)
  signac_object$nucleosome_group <- ifelse(signac_object$nucleosome_signal > 2, 'NS > 2', 'NS < 2')
  signac_object <- TSSEnrichment(object = signac_object, fast = FALSE)
  signac_object$high.tss <- ifelse(signac_object$TSS.enrichment > 2, 'High', 'Low')
  signac_object <- subset(
    x = signac_object,
    subset = nCount_peaks > 2000 &
      nCount_peaks < 20000 &
      nFeature_peaks < 10000 &
      nucleosome_signal < 2 &
      TSS.enrichment > 2)
}
postfilter_atac_processing <- function(){
  signac_object <- FindTopFeatures(signac_object, min.cutoff = 10)
  signac_object <- RunTFIDF(signac_object, method = 3)
  signac_object <- RunSVD(signac_object)
  signac_object <- RunUMAP(signac_object, reduction = "lsi", dims = 2:30)
  signac_object <- FindNeighbors(object = signac_object, reduction = 'lsi', dims = 2:30)
  signac_object <- FindClusters(object = signac_object, verbose = FALSE, algorithm = 3)
}
process_atac <- function(signac_object) {
  # extract gene annotations from EnsDb
  library(EnsDb.Hsapiens.v86)
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  # change to UCSC style since the data was mapped to hg19
  seqlevelsStyle(annotations) <- 'UCSC'
  Annotation(signac_object) <- annotations
  signac_object <- NucleosomeSignal(object = signac_object)
  signac_object$nucleosome_group <- ifelse(signac_object$nucleosome_signal > 2, 'NS > 2', 'NS < 2')
  signac_object <- TSSEnrichment(object = signac_object, fast = FALSE)
  signac_object$high.tss <- ifelse(signac_object$TSS.enrichment > 2, 'High', 'Low')
  signac_object <- subset(
    x = signac_object,
    subset = nCount_peaks > 2000 &
      nCount_peaks < 20000 &
      nFeature_peaks < 10000 &
      nucleosome_signal < 2 &
      TSS.enrichment > 2)
  signac_object <- FindTopFeatures(signac_object, min.cutoff = 10)
  signac_object <- RunTFIDF(signac_object, method = 3)
  signac_object <- RunSVD(signac_object)
  signac_object <- RunUMAP(signac_object, reduction = "lsi", dims = 2:30)
  signac_object <- FindNeighbors(object = signac_object, reduction = 'lsi', dims = 2:30)
  signac_object <- FindClusters(object = signac_object, verbose = FALSE, algorithm = 3)
}
# Add the gene information to the object
gene_activity_matrix <- function(signac_object){
  gene.activities <- GeneActivity(signac_object)
  # add the gene activity matrix to the Seurat object as a new assay
  signac_object[['RNA']] <- CreateAssayObject(counts = gene.activities)
  # normalize RNA assay
  signac_object <- NormalizeData(
    object = signac_object,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(signac_object$nCount_RNA)
  )
}
# Read and process saved nebula results
read_nebula_results <- function(path){
  library(tidyverse)
  library(data.table)
  files <- list.files(path, pattern = "*.csv")
  results <- list()
  for (file in files){
    cell <- gsub('.csv','',basename(file))
    nebula_res <-
      data.table::fread(file.path(path, file), sep = 'auto')
    nebula_res <- nebula_res |> dplyr::select(-'V1')
    nebula_res$cell <- cell
    nebula_res$p_adj <- p.adjust(nebula_res$p_groupneonatal, method = "BH")
    nebula_res <- nebula_res |> dplyr::filter(p_adj < 0.05)
    results[[cell]] <- nebula_res
  }
  return(results)
}

# internal
run_ChIPseeker <- function(peaks){
  edb <- EnsDb.Hsapiens.v86
  seqlevelsStyle(edb) <- "UCSC"
  pk <- StringToGRanges(peaks)
  peakAnno.edb <- annotatePeak(
    pk,
    tssRegion = c(-3000, 3000),
    TxDb = edb,
    annoDb = "org.Hs.eg.db"
  )
  peak_annotation <- annoGR2DF(peakAnno.edb@anno)

  # create peak column in the dataframe
  peak_annotation$peak <- paste(
    peak_annotation$chr,
    peak_annotation$start,
    peak_annotation$end,
    sep = "-"
  )
  peak_annotation
}

select_closest_peak_to_gene <- function(annotation_dataframe) {
  # Select few columns and put peaks as row names
  peak_annotation_short <-
    annotation_dataframe |>
    dplyr::mutate(distance = abs(distanceToTSS)) |>
    dplyr::select(peak, SYMBOL, annotation, distanceToTSS, distance) |>
    remove_rownames() |>
    column_to_rownames("peak")

  # Arrange Genes by distance to TSS
  peak_annotation_short <-
    peak_annotation_short[order( peak_annotation_short[,1],
                                 peak_annotation_short[,4] ),]

  # Remove duplicates by chosing closest gene-peak
  peak_annotation_short <-
    peak_annotation_short[!duplicated(peak_annotation_short$SYMBOL), ]  # Extract unique rows

  # Select only gene column and peaks as row names
  peak_annotation_short |>
    dplyr::select(SYMBOL)
}

annotate_peaks_pseudobulk <- function(pseudobulk_df) {
  # get peaks (rownames of pseudobulk df)
  peak_annotation <- run_ChIPseeker(rownames(pseudobulk_df))
  peak_annotation_short <- select_closest_peak_to_gene(peak_annotation)

  # Merge scaled peak matrix with gene names
  # Put gene names as row names
  merge(pseudobulk_df, peak_annotation_short, by = 0) |>
    dplyr::filter(!is.na(SYMBOL)) |>
    remove_rownames() |>
    column_to_rownames("SYMBOL") |>
    dplyr::select(!Row.names)
}

closest_gene_peak_pair <- function(pseudobulk_df){
  library(fgsea)
  library(Repitools)
  library(ChIPseeker)
  library(EnsDb.Hsapiens.v86)
  rnames_at <- rownames(pseudobulk_df)
  edb <- EnsDb.Hsapiens.v86
  seqlevelsStyle(edb) <- "UCSC"
  pk <- Signac::StringToGRanges(gsub("\\.", "-", rnames_at))

  # Peak annotation with ChIPSeeker
  peakAnno.edb <-
    annotatePeak(pk, tssRegion=c(-3000, 3000),
                 TxDb=edb, annoDb="org.Hs.eg.db")
  peak_annotation <- annoGR2DF(peakAnno.edb@anno)

  # Create peak column in the dataframe
  peak_annotation$peak <-
    paste(peak_annotation$chr,
          peak_annotation$start,
          peak_annotation$end, sep = "-")

  # Select few columns and put peaks as row names
  peak_annotation_short <-
    peak_annotation |>
    dplyr::mutate(distance = abs(distanceToTSS)) |>
    dplyr::select(peak, SYMBOL, annotation, distanceToTSS, distance)

  peak_annotation_short <-
    peak_annotation_short[!duplicated(peak_annotation_short$peak),]

  peak_annotation_short <-
    peak_annotation_short |>
    remove_rownames() |>
    column_to_rownames(var="peak")

  # Arrange Genes by distance to TSS
  peak_annotation_short <-
    peak_annotation_short[order( peak_annotation_short[,1],
                                 peak_annotation_short[,4] ),]
  # Remove duplicates by chosing closest gene-peak
  peak_annotation_short <-
    peak_annotation_short[!duplicated(peak_annotation_short$SYMBOL), ]
  # Extract unique rows

  # Select only gene column and peaks as row names
  peak_annotation_short <-
    peak_annotation_short |>
    dplyr::select(SYMBOL)
  return(peak_annotation_short)
}
label_imputation <- function(tissue){
  require(tidyverse)
  require(data.table)
  require(VIM)
  path_label <- '/media/AGFORTELNY/people/rohit/projects/epigenetic_potential/scglue'
  combined_obs <- read.csv(file.path(path_label, paste0(tissue,'_combined_obs.csv')))
  combined.obsm <- read.csv(file.path(path_label, paste0(tissue,'_combined_obsm.csv')),
                            header = FALSE)
  names(combined.obsm) <- c('UMAP1', 'UMAP2')
  combined_obs['UMAP1'] <- combined.obsm['UMAP1']
  combined_obs['UMAP2'] <- combined.obsm['UMAP2']
  data <-
    combined_obs |>
    mutate(name = X) |>
    dplyr::select(name, rna_annotation, UMAP1, UMAP2)

  data[data$rna_annotation =="", 'rna_annotation'] <- NA
  data$tech <- ifelse(is.na(data$rna_annotation), 'ATAC', 'RNA')
  knn_imp <- kNN(data, impNA = TRUE)
  knn_atac <- knn_imp[knn_imp$tech =='ATAC', ]
  return(knn_atac)
}
save_peaks_by_celltype <- function(limma_res, peak_gene_pair_df, outpath){
  peak_gene_pair_df <-
    tibble::rownames_to_column(peak_gene_pair_df, "peak")
  res_filt <-
    limma_res |>
    dplyr::filter(adj.P.Val < 0.05)
  res_filt_split <- split(res_filt, res_filt[,"cell"])
  for (name in names(res_filt_split)){
    message(name)
    print(dim(res_filt_split[[name]]))
    if (dim(res_filt_split[[name]])[1] >200){
      goi <- res_filt_split[[name]]$rn
      #print(head(goi))
      peaks <- peak_gene_pair_df[peak_gene_pair_df$SYMBOL %in% goi, "peak"]
      peaks_gr <- Signac::StringToGRanges(peaks)
      GenomeInfoDb::sortSeqlevels(peaks_gr)
      peaks_gr <- sort(peaks_gr)
      print(peaks_gr)
      df <-
        data.frame(seqnames=seqnames(peaks_gr),
                   starts=start(peaks_gr)-1,
                   ends=end(peaks_gr))
      if (!dir.exists(outpath)) {dir.create(outpath, recursive = TRUE)}
      write.table(
        df, file=file.path(outpath, paste0(name,'.bed')),
        quote=F, sep="\t", row.names=F, col.names=F)
    }
  }
}
create_homer_data <- function(outpath, tissue, tissue_pseudo_bulk, tissue_res){
  # create outpath
  base_path <- paste0(
    outpath, '/',tissue,'/','motifs')
  # Get peak-gene pair
  psatac_peak_gene_pair <- closest_gene_peak_pair(tissue_pseudo_bulk)
  # save peaks by cell type in bed files
  save_peaks_by_celltype(tissue_res, psatac_peak_gene_pair, base_path)
}
homer_db <- function(tissue_path){
  path <- tissue_path
  directories <- list.dirs(path = path, recursive = FALSE, full.names = FALSE)
  motif_df_list <- list()
  for (dir in directories){
    #message(dir)
    cur_dir <- file.path(path, dir)
    if (file.exists(file.path(cur_dir, "knownResults.txt" ))){
      #message(file.path(cur_dir, "knownResults.txt" ))
      infile <- file.path(cur_dir, "knownResults.txt" )
      kr <- read.delim(infile, sep = '\t', header = TRUE)
      names(kr) <- gsub('[[:digit:]]+', '', names(kr))
      kr <-
        kr |>
        mutate(cell = dir) |>
        separate(Motif.Name, sep = '/',
                 c('motif_name', 'motif_source', 'database')) |>
        separate(motif_name, sep = ',',
                 c('name', 'subfamily')) |>
        mutate(name = str_replace_all(name, "\\(", "_")) |>
        mutate(name = str_replace_all(name, "\\)", "")) |>
        separate(name, sep = '_', c('motif_name','family')) |>
        unite(name, c('motif_name','family'), sep = "_", remove = FALSE) |>
        mutate(X..of.Target.Sequences.with.Motif = str_replace_all(
          X..of.Target.Sequences.with.Motif, "%", ""))  |>
        mutate(X..of.Background.Sequences.with.Motif = str_replace_all(
          X..of.Background.Sequences.with.Motif, "%", "")) |>
        mutate(odds_ratio = as.double(X..of.Target.Sequences.with.Motif)/as.double(X..of.Background.Sequences.with.Motif)) |>
        dplyr::filter(`P.value` < 0.05 )
      motif_df_list[[dir]] <- kr
    }
  }
  motif_db <- do.call(rbind, motif_df_list)
  motif_db <- motif_db |> arrange(family)
  motif_db <- motif_db[order(motif_db$family),]
  return(motif_db)
}
motif_dotplot <- function(homer_database, tissue){
  homer_plot <-
    homer_database |>
    ggplot(aes(x = cell, y = forcats::fct_reorder(name, family))) +
    geom_point(color = -log10(homer_database$P.value),
               size = homer_database$odds_ratio) +
    geom_point(data=homer_database %>% dplyr::filter(`q.value..Benjamini.` < 0.05),
               pch=21, size= 4, stroke = 1.5, colour="navy") +
    theme_minimal() +
    theme(axis.text.x = element_text(
      angle = 90, vjust = 0.5, hjust=1)) +
    #theme(axis.text.y = element_text(colour = c("yellow", "blue"))) +
    #scale_colour_gradientn(colours = terrain.colors(10)) +
    labs(x = 'cell type', y = 'motif', title = tissue,
         caption = 'Color by P-value')
  return(homer_plot)
}
