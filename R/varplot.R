##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param vpobj
##' @param id
##' @param mes
##' @param meta
##' @param outfolder
##' @return
##' @author dylanmr
##' @export
varplot <- function(wgcna, object_name) {
  
  library(BiocParallel)
  
  wgcna %>%
    mutate(
      # run variance partition
      varpart = map2(mergedMEs, meta, function(mod, met) {
        form <- ~ (1|trt) + (1|orig.ident)
        modules <- t(mod[["newMEs"]])
        variancePartition::fitExtractVarPartModel(modules, form, met, BPPARAM = SerialParam())
      }),
      plot = pmap(list(varpart, mergedMEs, meta, id, object_name), function(x, y, z, j, k) {
        if (dim(data.frame(x))[1] < 10) {
          a <- plotPercentBars(x[order(-x[["trt"]]), ][1:dim(data.frame(x))[1], ])
        } else {
          a <- plotPercentBars(x[order(-x[["trt"]]), ][1:10, ])
        }
        b <- plotVarPart(x)
        trt_max <- which.max(x[["trt"]])
        data.frame( Expression = y[["newMEs"]][,trt_max], trt = z[["trt"]], sample = z[["orig.ident"]]) %>%
          group_by(trt) %>%
          sample_frac(size = .25) -> GE
        bymod <- ggplot(GE, aes(x=fct_relevel(trt, "PF","R50","FGF"), y=Expression, fill= sample)) + 
          ggbeeswarm::geom_quasirandom(alpha=0.5, shape=21) + theme_classic() + 
          ggtitle(paste0("ME:", colnames(y[["newMEs"]])[trt_max]))
        
        samp_max <- which.max(x[["orig.ident"]])
        data.frame( Expression = y[["newMEs"]][,samp_max], trt = z[["trt"]], sample = z[["orig.ident"]]) %>%
          group_by(trt) %>%
          sample_frac(size = .25) -> GE
        bysamp <- ggplot(GE, aes(x=sample, y=Expression, color= trt)) + 
          ggbeeswarm::geom_quasirandom() + theme_classic() + 
          ggtitle(paste0("ME:", colnames(y[["newMEs"]])[samp_max]))
        
        top <- cowplot::plot_grid(a, b, axis = "tb", align = "hv", rel_widths = c(1.5, 1))
        bot <- cowplot::plot_grid(bymod, bysamp, axis = "tb", align = "hv")
        
        plot <- cowplot::plot_grid(top, bot, ncol = 1, align = "hv", axis = "tbr")
        title <- cowplot::ggdraw() +
          draw_label(paste0("Variance Partition of ", j),
                     fontface = "bold",
                     x = 0, hjust = 0
          ) +
          theme(plot.margin = margin(0, 0, 0, 7))
        cowplot::plot_grid(
          title, plot,
          ncol = 1,
          # rel_heights values control vertical title margins
          rel_heights = c(0.1, 1))
        k <- paste0(str_remove(k, pattern = "wgcna_"),"_")
        outname <- paste0("output/wgcna/varplots/", k,j, "_vp_plot.pdf")
        ggsave(filename = file_out(outname), w=10, h=8)
      }),
      id=id
    )} 