##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param data nested tibble
##' @param networkType signed, hybrid, unsigned; default is signed
##' @param clust_method clustering method for hclust; use complete for method rather than average (gives better results)
##' @param pamStage
##' @param merge_thresh

##' @return
##' @author dylanmr
##' @export
##'

# cant figure out how to use future map in this code

get_wgcna_modules <- function(model_prepped = model_prepped, networkType = "signed",
                              var.feats = 3000, clust_method = "complete",
                              pamStage = F, rsq_cut = 0.8, deepSplit = 4,
                              minClusterSize = 15, merge_thresh = 0.2) {

  print(output_folder)
  
  # create output folder for specific parameters
  outfolder <- sprintf("%swgcna/cm_%s_ps_%s_ds_%s_mms_%s/", output_folder, clust_method, pamStage, deepSplit, minClusterSize)
  dir.create(outfolder, showWarnings = FALSE)

  # generate rsq data frame and remove highly connected rsq vals
  model_prepped %>%
    transmute(
      # remove columns with metadata
      dat = map(data, ~ dplyr::select(.x, 1:all_of(var.feats))),

      # generate data frame with information on rsq, connectivity and pwr
      fitindex = map(dat, function(x) {
        pickSoftThreshold(
          data = x, dataIsExpr = TRUE, verbose = 1,
          powerVector = 1:40, corOptions = list(use = "p"),
          networkType = networkType
        ) %$%
          fitIndices %>%
          data.frame() %>%
          janitor::clean_names()
      }),

      # convert cell names to clean
      id = map(value, janitor::make_clean_names),

      # generate scale free topology plots and dont save in df
      sft_plot = walk2(fitindex, id, function(x, y) {
        pdf(paste0(outfolder, y, "_scalefreetop.pdf"))
        cex1 <- 0.9
        plot(x[["power"]], -sign(x[["slope"]]) * x[["sft_r_sq"]],
          xlab = "Soft Threshold (power)",
          ylab = "Scale Free Topology Model Fit, signed R^2",
          type = "n", main = paste("Scale independence")
        )
        text(x[, 1], -sign(x[, 3]) * x[, 2],
          cex = cex1, col = "red"
        )
        abline(h = 0.80, col = "red")
        dev.off()
      }),

      # select pwr based on rsq cutoff
      pwr = map(fitindex, function(x) {
        filt <- dplyr::filter(x, median_k <= quantile(median_k, 0.95, na.rm = T))
        if (sum(filt[["sft_r_sq"]] >= rsq_cut) > 1) {
          filt <- filt[filt[["sft_r_sq"]] >= rsq_cut, ]
          filt[which.min(filt$power), 1]
        } else {
          filt %>%
            dplyr::arrange(desc(sft_r_sq)) %>%
            dplyr::slice(1) %>%
            pull(power)
        }
      }),

      # generate tom
      tom = map2(pwr, dat, ~ generate_tom(sft = .x, dat = .y, networkType = networkType)),

      # generate gene tree
      gt = map(tom, ~ hclust(as.dist(.x), method = clust_method)),

      # generate modules
      mods = map2(gt, tom, ~ cutreeDynamic(
        dendro = .x, distM = as.matrix(.y),
        method = "hybrid", pamStage = pamStage,
        deepSplit = deepSplit,
        minClusterSize = minClusterSize
      )),

      # generate module dendrogram plot
      ## not sure why pwalk fails - for now removing column at end of function
      dend_plot = purrr::pmap(list(gt, mods, id), function(x, y, z) {
        pdf(paste0(outfolder, z, "_dend_mods.pdf"))
        WGCNA::plotDendroAndColors(x, y, "Dynamic Tree Cut",
          dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
          main = "Gene dendrogram and module colors"
        )
        dev.off()
      }),

      # calculate ME expression
      mes = map2(dat, mods, ~ WGCNA::moduleEigengenes(.x, .y)$eigengenes),

      # plot similarity of modules
      sim_mods = walk2(mes, id, function(x, y) {
        pdf(paste0(outfolder, y, "_mod_simil.pdf"))
        me_dis <- 1 - cor(x)
        me_tree <- hclust(as.dist(me_dis), method = "average")
        plot(me_tree, main = "Clustering of module eigengenes", xlab = "", sub = "")
        merge_thresh <- merge_thresh
        abline(h = merge_thresh, col = "red")
        dev.off()
      }),

      # merge modules on threshold of similarity
      mergedMEs = map2(dat, mods, ~ mergeCloseModules(.x, .y, cutHeight = merge_thresh, verbose = 1)),

      # calculate importance of each gene to each module
      kme = map2(dat, mergedMEs, ~ WGCNA::signedKME(.x, .y$newMEs)),

      # generate csv file of module genes and
      gene_mod_df = pmap(list(dat, mergedMEs, kme, id), function(x, y, z, a) {
        y[["colors"]] %>%
          # add ME to module name
          paste0("ME", .) %>%
          # assign each gene to module
          bind_cols(mod = ., gene = colnames(x)) %>%
          # create a gene module combined column
          unite(gene, mod, col = "match") %>%
          # join with kme values
          left_join(z %>%
            # save gene names
            rownames_to_column("gene") %>%
            # get kme value for each gene for each module
            pivot_longer(cols = contains("ME")) %>%
            # extract module names
            mutate(ME = substring(name, 2)) %>%
            # join gene and module name and only extract matches with gene module pairs from above
            unite(gene, ME, col = "match", remove = F)) %>%
          # extract a table of gene, module and kme values
          dplyr::select(gene, ME, value) %>%
          filter(ME != "ME0") %>%
          arrange(ME, -value) %>%
          write_csv(path = paste0(outfolder, a, "gene_mods.csv"))
      }),
      # save metadata column for downstream analysis
      meta = meta
    ) %>%
    dplyr::select(-dat, -dend_plot, -sft_plot)
}
