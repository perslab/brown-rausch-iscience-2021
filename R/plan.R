run_names =  c("34_R50E", "4_R50E",
               "12_R50E","37_R50E", 
               "30_R50E", "3_R50E")

data_folder = "/data/sc-10x/data-runs/170904-schwartz-fgf1/"
output_folder = "/raid5/home/cbmr/lhv464/rausch_brown_r50/output/"
epath = paste0(data_folder, run_names, "-5000_cells/outs/filtered_gene_bc_matrices/mm10")
glia_path = "/projects/dylan/bentsen-rausch-2019/data/glia/glia_seur_filtered.RDS"

      
the_plan <-
  
  drake_plan(
    
    # read in data and pre-process/cluster
    srt = target(gen_seur(endogen = epath, run = run_names, aligner = "cr") %>%
                   Seurat::AddMetaData(., run_scrublet(., run_names)$dbl_call, col.name = "dbl_call") %>%
                   Seurat::AddMetaData(., run_scrublet(., run_names)$dbl_score, col.name = "dbl_score") %>%
                   quickclus(),
                transform = map(epath = !!epath, run_names = !!run_names, .id = run_names)),
    
    # processed FGF1 and PF data
    ext_data = readRDS("/projects/dylan/bentsen-rausch-2019/data/filtglia.RDS"),
    
    # project labels from campbell for cluster assignment
    campref = gen_campbell_ref(),
  
    # merge r50 objects, remove neurons and RBC
    merged_r50 = target(merge(x = list(srt)[[1]], y = c(list(srt)[-1]), add.cell.ids = !!run_names, merge.data = T) %>%
                          subset(., subset = dbl_call == F) %>%
                          quickclus(seurobj = ., clusid = "orig.ident", idents = c("12_R50E", "37_R50E"), invert=T) %>%
                          subset(., subset=SCT_snn_res.0.2 %in% c(3,11), invert=T),
                         transform = combine(srt)),
    
    # adjusted
    merged_seur = merge(merged_r50, ext_data) %>% 
      quickclus() %>%
      set_ident(., "0.2") %>%
      RenameIdents(., .[[]] %>% 
                     dplyr::count(orig_labels, SCT_snn_res.0.2) %>% 
                     dplyr::group_by(SCT_snn_res.0.2) %>% 
                     top_n(1) %>% dplyr::select(2,1) %>% 
                     mutate(lab = as.character(orig_labels)) %>% 
                     ungroup() %>% dplyr::select(1,3) %>% deframe()),
    
    # need to pass labels from pre-labeled to new samples, then run sil to remove outliers
    # Compute silhouette coefficient
    filt_glia = runsil(merged_seur),
    
    # overview of data 
    
    # target_name = target(
    #   command = {
    #     rmarkdown::render(knitr_in("doc/data_overview.Rmd"))
    #     file_out("doc/data_overview.html")
    #   }
    # )
    
    # build renormalized tibbles
    ## goal: prep data for modelling
    ## steps: take each celltype, subset, run sctransform, store data, store metadata
    # input: glia_cluster
    model_prepped = gen_model_obj(filt_glia, "orig_labels"),
    
    # use wgcna to identify sets of coexpressed genes
    wgcna = target(get_wgcna_modules(model_prepped, clust_method = cm, pamStage = ps, deepSplit = ds, minClusterSize = mms),
                   transform = cross(
                     cm = c("complete", "average"),
                     ps = c(T,F),
                     ds = c(2,4),
                     mms = c(5,15,50)
                   )
    ),

    varplot = target(varplot(wgcna, .id_chr), transform = map(wgcna))
    
    # run pseudobulk analysis
    
  
)




# ,
# 
#     # run classifier to identify and split glia/neurons
#     # define input grid above
#     neur_classifier = classify_by_geneset(seurobj = merged_seur, 
#                                                  features = list(c("Snhg11","Snap25")),
#                                                  names = "neur", 
#                                                  k = 2, 
#                                                  group = "SCT_snn_res.0.8"),
#     
#     # process campbell dataset for projection
#     neur_ref = gen_camp_neuron_ref(),
#     
#     # process fgf1 nuclei for glia projection
#     glia_ref = readRDS(glia_path) %>% 
#       SCTransform() %>%
#       RunPCA(),
#     
#     # recluster
#     # define input grid above
#     reclustered = target(quickclus(seurobj = merged_seur, clusid = "SCT_snn_res.0.8", 
#                                    invert = invert, idents = neur_classifier %>%
#                                      dplyr::filter(assignment=="incr") %>%
#                                      pull(id)) %>% 
#                            project_labels(query = ., ref = celltype_ref) %>%
#                            set_ident(., "1") %>%
#                            RenameIdents(., rename_by_camp(.)),
#                          transform = map(invert = c(T,F),  celltype_ref = c(glia_ref, neur_ref),
#                                          .names = c("glia_cluster","neur_cluster"))),
#     
#     # use silhouette coefficient to remove ambiguous cells 
#     
#     
#     # build renormalized tibbles
#     ## goal: prep data for modelling
#     ## steps: take each celltype, subset, run sctransform, store data, store metadata
#     # input: glia_cluster
#     model_prepped = gen_model_obj(filt_glia),
#     
#     
