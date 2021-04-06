##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param run
##' @return
##' @author dylanmr
##' @export


gen_seur <- function(hto = NULL, endogen, run, aligner = "kallisto", meta_path = NULL,
                     min.features = 400, min.cells = 10, mito = 20, ribo = 20) {
  require(Seurat)
  require(magrittr)
  require(AnnotationDbi)
  require(org.Mm.eg.db)

  if (aligner == "kallisto") {

    # Read in kallisto counts
    require(tidyverse)
    ct_mat <- BUSpaRse::read_velocity_output(
      spliced_dir = file_in(endogen),
      spliced_name = "spliced",
      unspliced_dir = file_in(endogen),
      unspliced_name = "unspliced"
    )

    # Calculate barcode rank and classify cells
    tot_count <- Matrix::colSums(ct_mat$unspliced)
    bc_rank <- DropletUtils::barcodeRanks(ct_mat$spliced)
    bc_uns <- DropletUtils::barcodeRanks(ct_mat$unspliced)

    # Filter matrices on dropletUtils
    bcs_use <- names(tot_count)[tot_count > bc_uns$inflection]
    tot_genes <- Matrix::rowSums(ct_mat$unspliced)

    # Filter to keep genes found at least 1x
    genes_use <- rownames(ct_mat$unspliced)[tot_genes > 0]
    sf <- ct_mat$spliced[genes_use, bcs_use]
    uf <- ct_mat$unspliced[genes_use, bcs_use]

    # Generate matrix of summed counts
    summed <- sf + uf

    # Convert rownames from ens to gene names
    ann <- data.frame(mapIds(org.Mm.eg.db, sapply(strsplit(rownames(summed), "[.]"), "[", 1),
      column = "SYMBOL", keytype = "ENSEMBL"
    ))
    ann$ens <- as.character(rownames(ann))
    colnames(ann)[1] <- "gene"
    ann$gene <- as.character(ann$gene)
    converted <- ifelse(is.na(ann$gene), ann$ens, ann$gene)

    # rename matrices
    rownames(sf) <- converted
    rownames(uf) <- converted
    rownames(summed) <- converted

    if (!is.null(hto)) {

      # Read HTO mat
      hto <- BUSpaRse::read_count_output(dir = file_in(hto), name = "cells_x_features", tcc = FALSE)

      # filter hto and endo libraries to contain same cells
      joint.bcs <- dplyr::intersect(colnames(hto), colnames(summed))
      summed <- summed[, joint.bcs]
      hto <- as.matrix(hto[, joint.bcs])
      seur <- CreateSeuratObject(counts = summed, min.cells = min.cells, project = paste0("CP_", run))
      seur[["hto"]] <- CreateAssayObject(hto)

      # need to think about quantile here
      seur %>%
        NormalizeData(assay = "hto", normalization.method = "CLR") %>%
        HTODemux(assay = "hto", positive.quantile = 0.99) -> seur
      
    } else {
      
      seur <- CreateSeuratObject(counts = summed, min.cells = min.cells, min.features = min.features, project = paste0("CP_", run))
      seur[["percent.mt"]] <- PercentageFeatureSet(seur, pattern = "^mt-")
      seur[["percent.ribo"]] <- PercentageFeatureSet(seur, pattern = "^Rp[sl][[:digit:]]")
      seur <- subset(seur, subset = percent.mt < mito & percent.ribo < ribo)
      
    }

    # Build seurat object with assays for spliced and unspliced data

    sf %>%
      set_rownames(make.unique(rownames(sf))) %>%
      .[rownames(sf), colnames(summed)] -> spliced
    seur[["spliced"]] <- CreateAssayObject(counts = spliced)

    uf %>%
      set_rownames(make.unique(rownames(uf))) %>%
      .[rownames(seur), colnames(summed)] -> unspliced

    seur[["unspliced"]] <- CreateAssayObject(counts = unspliced)
    
  } else if (aligner == "cr") {
    
    mat <- Seurat::Read10X(data.dir = file_in(endogen))

    if (!is.null(hto)) {

      # Read HTO mat
      hto <- BUSpaRse::read_count_output(dir = file_in(hto), name = "cells_x_features", tcc = FALSE)

      # filter hto and endo libraries to contain same cells
      joint.bcs <- dplyr::intersect(colnames(hto), colnames(mat))
      mat <- mat[, joint.bcs]
      hto <- as.matrix(hto[, joint.bcs])

      # create seurat object (cannot filter out cells at this step with HTO - assays will no longer match in size)
      seur <- CreateSeuratObject(mat, min.cells = min.cells, project = run)
      seur[["hto"]] <- CreateAssayObject(hto)

      # need to think about quantile here
      seur %>%
        NormalizeData(assay = "hto", normalization.method = "CLR") %>%
        HTODemux(assay = "hto", positive.quantile = 0.99) -> seur
      
    } else {
      
      seur <- CreateSeuratObject(mat, min.cells = min.cells, min.features = min.features, project = run)
      seur[["percent.mt"]] <- PercentageFeatureSet(seur, pattern = "^mt-")
      seur[["percent.ribo"]] <- PercentageFeatureSet(seur, pattern = "^Rp[sl][[:digit:]]")
      seur <- subset(seur, subset = percent.mt < mito & percent.ribo < ribo)
      
    }
  }


  # get metadata from googlesheet
  if (!is.null(meta_path)) {
    
    meta <- get_metadata(meta_path)
    # join metadata with seurat object
    seur@meta.data <- cbind(seur@meta.data, meta[match(as.character(seur$hash.ID), meta$hash.ID), ])
    
  }

  # return final object
  return(seur)
}
