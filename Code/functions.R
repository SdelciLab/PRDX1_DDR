#### functions used in this project ####
Load_DIA_NN_Data <- function(report_DIA_tsv_file, Samples_df){
    read_tsv(here::here("Datasets","Raw",report_DIA_tsv_file)) %>% 
        as.data.frame() %>% 
        subset(Lib.Q.Value<= 0.01 & Lib.PG.Q.Value <= 0.01 ) %>% janitor::clean_names() %>% 
        dplyr::select(matches("pg_max_lfq|run|protein_group|^genes$")) %>% 
        distinct() %>% 
        subset(!str_detect(protein_group,paste0(contaminant_uniprots,collapse = "|"))) %>% 
        remove_rownames() %>%
        
        mutate(run = str_replace_all(run, c("P11690-1" = "P11591-1",
                                            "P11690-2" = "P11592-1",
                                            "P11690-3" = "P11593-1"))) %>% 
      mutate(run = run %>% str_match("-([:graph:]{6})-.$") %>% .[,2] %>% tolower()) %>%  
  left_join(Samples_df %>% 
              mutate(across(everything(),tolower))) %>% 
        pivot_wider(-run, names_from = Condition ,values_from = pg_max_lfq) %>% 
     
        # group_split(protein_group) %>%
        # map_dfr(.x = .,~.x %>% discard_uniprot_isoforms("protein_group","genes" )) %>%
        # group_by(protein_group) %>%
        # mutate(Is_duplicated = n() > 1,
        #        New_Uniprot = if_else(Is_duplicated == F,protein_group,old_pg)) %>% 
        column_to_rownames("protein_group")%>% 
        dplyr::select(any_of(Samples_df$Condition)) %>% 
        # discard_single_isoforms %>% 
        as.matrix() #%>% log2()
    
    
    
    
}
map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

DEP_DIA <- function(input_matrix,dataset_name){
    #this required a not normalised and not logged matrix with column names of condition in _rep format
   
         
  
  dataset_name = list_files_to_analyse[1][[1]]
            input_matrix <- Methods_DIA[[1]]

       input_matrix <- input_matrix %>% as.data.frame()
       set.seed(2023)
       experimental_design_DIA <-  data.frame(
        label = colnames(input_matrix),
        condition =  str_remove_all(colnames(input_matrix),"_[:graph:]*$"),
        replicate = str_remove_all(colnames(input_matrix),"^[:graph:]*_") %>% as.numeric()
    )
    data_unique_Etop <- input_matrix %>% rownames_to_column("name") %>% 
        left_join(HUMAN_9606 %>% 
                      subset(Type == "Gene_Name") %>% 
                      dplyr::select(-Type) %>% 
                    subset(!duplicated(Uniprot)), 
                  by = c("name" = "Uniprot")) %>%
        subset(!duplicated(name))
    # performs batch correction
    if((experimental_design_DIA %>% group_by(condition) %>% 
       dplyr::summarise(Total_Rep = length(unique(replicate))) %>% 
       pull(Total_Rep) %>% min())>3){
      n_accepted_NA_per_condition  = 3
     ComBAT = T
    }else{
      n_accepted_NA_per_condition  = 1
       ComBAT = F
    }
      output_folder = 'Output'
    Quant_columns <- which(colnames(data_unique_Etop) %in%colnames(input_matrix))# get LFQ column numbers
    data_se <- make_se(data_unique_Etop, Quant_columns, experimental_design_DIA)
    data_se_parsed <- make_se_parse(data_unique_Etop, Quant_columns)
    plot_frequency(data_se)+ggtitle(glue::glue("Protein_overlap ",dataset_name))
    ggsave(here::here('Output',glue::glue("Protein_overlap ",dataset_name,".pdf")))
    data_filt <- filter_missval(data_se, thr = 1)
    #data_filt2 <- filter_missval(data_se, thr = 1)
    plot_numbers(data_filt)+ggtitle(glue::glue("Protein_numbers ",dataset_name))
    ggsave(here::here(output_folder,glue::glue("Protein_numbers ",dataset_name,".pdf")))
    plot_coverage(data_filt)
    data_filt@assays@data@listData[[1]][is.nan(data_filt@assays@data@listData[[1]])] <- NA 
    pdf(here::here(output_folder,glue::glue("Protein_Missingness ",dataset_name,".pdf"))) 
    plot_missval(data_filt)
    dev.off()
    data_norm <- normalize_vsn(data_filt)
    DEP::meanSdPlot(data_norm)
    ggsave(here::here(output_folder,glue::glue("normalize_vsn ",dataset_name,".pdf")))
    
    plot_normalization(data_se, data_norm)+ggtitle(glue::glue("Protein_norm ",dataset_name))
    ggsave(here::here(output_folder,glue::glue("Protein_normalisation ",dataset_name,".pdf")))
    
    pca_res <- prcomp(data_norm@assays@data@listData[[1]]  %>% na.omit() %>% t(), scale=TRUE)
    var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)
    
    pca_res$x %>% 
        as.data.frame %>%
        rownames_to_column("Sample") %>% 
        mutate(Condition = str_remove(Sample,"_.$")) %>% 
        ggplot(aes(x=PC1,y=PC2, label = Sample, colour = Condition )) + geom_point(size=4) +
        ggrepel::geom_label_repel()+
        theme_bw(base_size=32) + 
        labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
             y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
        theme(legend.position="top") +
        ggtitle(dataset_name)+ 
        theme(plot.title = element_text(size = 20))
    ggsave(here::here(output_folder,glue::glue(dataset_name," PCA.pdf")))
    
    if(data_norm@assays@data@listData[[1]] %>% is.na() %>% any()){
        png(here::here(output_folder,glue::glue("Protein_Missingness_Abundance ",dataset_name,".pdf")), width = 2500, height = 3800,res  =300) 
        plot_detect(data_norm)
        dev.off()
        #ggsave(here::here(output_folder,glue::glue("Protein_missingness ",dataset_name,".pdf")))
          data_imp <- DEP::impute(data_norm, fun = "MinProb", q = 0.01)
         # to_impute <- data_norm@assays@data@listData[[1]] %>% subset(., rowSums(is.na(.))>n_accepted_NA_per_condition)
          # to_not_impute <- data_norm@assays@data@listData[[1]] %>% subset(., rowSums(is.na(.))<=(n_accepted_NA_per_condition))
          if(ComBAT == T){
            batch = if_else(str_match(colnames(input_matrix),"_(.)$")[,2] %>% as.numeric() <4,2,1)
            
            multiple_imputation <-  cbind(impute.mi(tab = data_norm@assays@data@listData[[1]][,batch == 2],#methodMNAR = "impute.pa",
                                                   conditions = experimental_design_DIA$condition[batch == 2] %>% as.factor(),
                                                   repbio = experimental_design_DIA$replicate[batch == 2] %>% as.factor()),
                  impute.mi(tab = data_norm@assays@data@listData[[1]][,batch == 1],#methodMNAR = "impute.pa",
                                                   conditions = experimental_design_DIA$condition[batch == 1] %>% as.factor(),
                                                   repbio = experimental_design_DIA$replicate[batch == 1] %>% as.factor())) 
            
          }else{multiple_imputation <-  impute.mi(tab = data_norm@assays@data@listData[[1]],#methodMNAR = "impute.pa",
                                                        conditions = experimental_design_DIA$condition %>% as.factor(),
                                                        repbio = experimental_design_DIA$replicate %>% as.factor())
          
            
          }
          
         rownames(multiple_imputation) <- rownames(data_norm@assays@data@listData[[1]])
         colnames(multiple_imputation) <- colnames(data_norm@assays@data@listData[[1]])
         data_imp@assays@data@listData[[1]] <- multiple_imputation#  rbind(to_not_impute,multiple_imputation)
         data_imp@assays@data@listData[[1]] <- data_imp@assays@data@listData[[1]][rownames(data_norm@assays@data@listData[[1]]),]
    }else{
        data_imp <- data_norm
    }
    plot_imputation(data_norm, data_imp)
    ggsave(here::here(output_folder,glue::glue("Protein_imputted ",dataset_name,".pdf")))
    
    if(ComBAT == T){
      edata = data_imp@assays@data@listData[[1]]
      data_norm@assays@data@listData[[1]] <- data_norm@assays@data@listData[[1]] %>%
        as.data.frame() %>% 
        dplyr::select(!matches("dmso_(4|5)")) %>% 
        as.matrix()

  
      
      # parametric adjustment
      combat_edata1 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
      
    
      input_matrix_batch <-2^combat_edata1 %>% as.data.frame() %>% 
        dplyr::select(!matches("dmso_(4|5)"))
      experimental_design_DIA <-  data.frame(
        label = colnames(input_matrix_batch),
        condition =  str_remove_all(colnames(input_matrix_batch),"_[:graph:]*$"),
        replicate = str_remove_all(colnames(input_matrix_batch),"^[:graph:]*_") %>% as.numeric()
      )
      data_unique_Etop <- input_matrix_batch %>% rownames_to_column("name") %>% 
        left_join(HUMAN_9606 %>% 
                    subset(Type == "Gene_Name") %>% 
                    dplyr::select(-Type) %>% 
                    subset(!duplicated(Uniprot)), 
                  by = c("name" = "Uniprot")) %>%
        subset(!duplicated(name))
      
      Quant_columns <- which(colnames(data_unique_Etop) %in%colnames(input_matrix_batch))# get LFQ column numbers
      data_se <- make_se(data_unique_Etop, Quant_columns, experimental_design_DIA)
      data_filt <- filter_missval(data_se, thr = 1)
      data_filt@assays@data@listData[[1]][is.nan(data_filt@assays@data@listData[[1]])] <- NA 
      data_norm_batch <- normalize_vsn(data_filt)
      plot_normalization(data_se,data_norm_batch)
         data_norm_batch@assays@data@listData[[1]][is.na(data_norm@assays@data@listData[[1]])] <- NA
     
      
      if(data_norm_batch@assays@data@listData[[1]] %>% is.na() %>% any()){
          data_imp <- DEP::impute(data_norm_batch, fun = "MinProb", q = 0.01)
        # to_impute <- data_norm@assays@data@listData[[1]] %>% subset(., rowSums(is.na(.))>n_accepted_NA_per_condition)
        # to_not_impute <- data_norm@assays@data@listData[[1]] %>% subset(., rowSums(is.na(.))<=(n_accepted_NA_per_condition)
          batch = if_else(str_match(colnames(input_matrix),"_(.)$")[,2] %>% as.numeric() <4,2,1)
          
          multiple_imputation <-  impute.mi(tab = data_norm_batch@assays@data@listData[[1]],#methodMNAR = "impute.pa",
                                                  conditions = experimental_design_DIA$condition %>% as.factor(),
                                                  repbio = experimental_design_DIA$replicate %>% as.factor())
          
        
        rownames(multiple_imputation) <- rownames(data_norm_batch@assays@data@listData[[1]])
        colnames(multiple_imputation) <- colnames(data_norm_batch@assays@data@listData[[1]])
        data_imp@assays@data@listData[[1]] <- multiple_imputation#  rbind(to_not_impute,multiple_imputation)
        data_imp@assays@data@listData[[1]] <- data_imp@assays@data@listData[[1]][rownames(data_norm_batch@assays@data@listData[[1]]),]
      }else{
        data_imp <- data_norm_batch
      }
         pca_res <- prcomp(data_imp@assays@data@listData[[1]]  %>% na.omit() %>% t(), scale=TRUE) 
         #plot_missval(data_filt)
         # Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
         
         var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)
         pca_res$x %>% 
           as.data.frame %>%
           rownames_to_column("Sample") %>% 
           mutate(Condition = str_remove(Sample,"_.$")) %>% 
           ggplot(aes(x=PC1,y=PC2, label = Sample, colour = Condition )) + geom_point(size=4) +
           ggrepel::geom_label_repel()+
           theme_bw(base_size=32) + 
           labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
                y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
           theme(legend.position="top") +
           ggtitle(glue::glue("Batch_Corrected",dataset_name))+ 
           theme(plot.title = element_text(size = 20))
         ggsave(here::here(output_folder,glue::glue(dataset_name,"batch PCA.pdf")))
         
      
      
    }
    data_diff_all_contrasts <- DEP::test_diff(data_imp, type = "all")
    dep <- add_rejections(data_diff_all_contrasts, alpha = 0.05, lfc = 0)
    
    Significant_protein_list <- data.frame(ProteinGroup = dep@elementMetadata$name,
                                           dmso_X0 = dep@elementMetadata$dmso_vs_x0h_significant,
                                           T0_X24 = dep@elementMetadata$x0h_vs_x24h_significant,
                                           dmso_X240 = dep@elementMetadata$dmso_vs_x24h_significant) %>% 
      mutate(Uniprot= ProteinGroup %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
      left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot)))
      
    
    Significant_proteins <- data_imp@assays@data@listData[[1]] %>% 
        as.data.frame() %>% 
        rownames_to_column("ProteinGroup") %>% 
        #mutate(Uniprot= ProteinGroup %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
        #left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot))) %>% 
        #left_join(Interesting_proteins) %>% 
        #subset(!is.na(Behaviour)) %>% 
        #pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
        #mutate(Condition= factor(Condition, levels= paste(rep(c("DMSO","T0","T24"), each= 3), rep(1:3,3),sep="_"))) %>% 
        #group_by(ProteinGroup) %>% pivot_wider(names_from = "Condition",values_from = Abundance) %>% 
        #mutate(DMSO = mean(c(DMSO_1,DMSO_2,DMSO_3))) %>% mutate(across(where(is.numeric), ~.x-DMSO)) %>%
        #dplyr::select(-DMSO) %>% pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
        left_join(dep@elementMetadata$significant %>% set_names(dep@elementMetadata$name) %>% enframe(name = "ProteinGroup", "Significant")) %>% 
        subset(Significant == T) %>% dplyr::select(-Significant) %>% 
        pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
      #   mutate(Condition = str_remove_all(Condition,"_.")) %>% 
      # ungroup() %>% 
      #   group_by(ProteinGroup,Condition) %>% 
        # dplyr::summarise(Mean_Abundance = mean(Abundance, na.rm = T)) %>% 
        pivot_wider(names_from = "Condition",values_from = "Abundance") %>% 
        mutate(Uniprot= ProteinGroup %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
        left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot))) %>% 
        #mutate(duplicated = BiocGenerics::duplicated(Uniprot))
        ungroup %>% 
        mutate(ID = if_else(duplicated(ID),paste0(ID,"_1"),ID)) %>% 
        mutate(ID = if_else(duplicated(ID),paste0(ID,"_1"),ID)) %>% 
        mutate(ID = if_else(duplicated(ID),paste0(ID,"_1"),ID)) %>% 
        subset(!is.na(ID)) %>% 
        column_to_rownames("ID") %>%
        dplyr::select(where(is.numeric)) %>%
         mutate(Rowmean = rowMeans(.),
                across(where(is.numeric),~.x- Rowmean)) %>% 
         dplyr::select(-Rowmean) %>% 
        as.matrix() 
    paletteLength <- 20
    myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
    # length(breaks) == length(paletteLength) + 1
    # use floor and ceiling to deal with even/odd length pallettelengths
    myBreaks <- c(seq(min(Significant_proteins, na.rm = T), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(Significant_proteins, na.rm = T)/paletteLength, max(Significant_proteins, na.rm = T), length.out=floor(paletteLength/2)))
    
    png(here::here(output_folder,glue::glue("Heatmap_Significant ",dataset_name,".pdf")), width = 2500, height = 3800,res  =300) 
    pheatmap::pheatmap(Significant_proteins[,experimental_design_DIA$label %>% sort()],cluster_cols = F,fontsize_row = 6, clustering_distance_rows = "euclidean", cluster_rows = T,
                       # scale = "row",
                       main = glue::glue(dataset_name, " \nSignificant Proteins Normalised to DMSO - imputted"),color=myColor, breaks=myBreaks)
    dev.off()

    Comparisons_list <- list()
    for(i in (dep@elementMetadata %>% names() %>% str_subset("diff") )){
                  # i = (dep@elementMetadata %>% names() %>% str_subset("diff"))[2]
        contrast <- str_remove_all(i,"_diff")
        print(contrast)
        conditions <- contrast %>% str_match("([:graph:]*)_vs_([:graph:]*)$") %>% .[2:3]
        non_missing_in_all_comparison <- c(data_norm@assays@data@listData[[1]] %>% as.data.frame() %>% dplyr::select(contains(conditions[1])) %>% 
            subset(., rowSums(is.na(.))<=n_accepted_NA_per_condition) %>% rownames(),
          data_norm@assays@data@listData[[1]] %>% as.data.frame() %>% dplyr::select(contains(conditions[2])) %>% 
              subset(., rowSums(is.na(.))<=n_accepted_NA_per_condition) %>% rownames()) %>% unique()
        Imputted <- c(data_norm@assays@data@listData[[1]] %>% as.data.frame() %>% dplyr::select(contains(conditions[1])) %>% 
                          subset(., rowSums(is.na(.))>1) %>% rownames(),
                      data_norm@assays@data@listData[[1]] %>% as.data.frame() %>% dplyr::select(contains(conditions[2])) %>% 
                          subset(., rowSums(is.na(.))>1) %>% rownames()) %>% unique()
        volcano_df <-  data.frame(log2_FC = dep@elementMetadata %>%  .[(glue::glue(contrast,"_diff"))] %>% unlist(),
                                  Uniprot = dep@elementMetadata$name,
                                  significant = dep@elementMetadata %>%  .[(glue::glue(contrast,"_significant"))] %>% unlist(),
                                  p.adj = dep@elementMetadata%>%  .[(glue::glue(contrast,"_p.adj"))] %>% unlist() ,
                                  p.val = dep@elementMetadata%>%  .[(glue::glue(contrast,"_p.val"))] %>% unlist()) %>% 
             subset(Uniprot %in% non_missing_in_all_comparison) %>% 
          
            mutate(Imputted_comparison = Uniprot %in% Imputted,
                   Single_Uniprot = Uniprot %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
          left_join(Metabolic_proteins, by  = c("Single_Uniprot" = "Uniprot") ) %>% 
          dplyr::rename(Metabolic_library = Behaviour) %>% 
            left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type), by  = c("Single_Uniprot" = "Uniprot")) %>% 
          mutate(Metabolic_library = if_else(is.na(Metabolic_library),"Non-Metabolic", Metabolic_library))
        Comparisons_list[[i]] <- volcano_df
    
        
    }
    
    
    return(list(Imputted = data_imp@assays@data@listData[[1]],
                Unimputted = data_norm@assays@data@listData[[1]], 
                DEPs = set_names(Comparisons_list, dep@elementMetadata %>% names() %>% str_subset("diff"))))
}
correp <- function(data_matrix,subset_selection){
  pacman::p_load(e1071,cluster,CORREP,NbClust,advclust) #
  #example from  https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/CORREP/inst/doc/CORREP.pdf
  #The format is = each column is a condition, each row is a Gene Replicate in that condiition and and all the replicates together
  # data_matrix <- Data_matrices$Imputted %>% dplyr::select(matches("_(1|2|3)$"))
   # subset_selection <- Significant_proteins
  # subset_selection = rownames(Data_matrices$Imputted)
   # data_matrix <-Data_matrices$Imputted  %>% 
   #   dplyr::select(matches("_(1|2|3)$"))
  data_matrix <- data_matrix[rownames(data_matrix) %in%  subset_selection, ]
  conditions = Data_matrices$Imputted  %>% colnames() %>%
    str_remove_all("_[:digit:]$") %>% unique()
  replicates <- Data_matrices$Imputted  %>% colnames() %>%
    str_match("_([:digit:]*)$") %>% .[,2] %>% as.numeric() %>% max()
  Pivotted_conditions <- map(.x= conditions, ~data_matrix %>% 
                               dplyr::select(contains(.x)) %>% 
                               rownames_to_column("Uniprot") %>% 
                               
                               pivot_longer(-Uniprot,names_to = "Replicate",values_to = .x) %>% 
                               mutate(Replicate = str_remove_all(Replicate,paste0(.x,"_"))) %>% 
                               unite(Uniprot,c(Uniprot,Replicate), sep = ", ")) %>% 
    purrr::reduce(left_join)
  
  d0.std <- apply(Pivotted_conditions %>% column_to_rownames("Uniprot"), 1, function(x) x/sd(x))
  input <- t(d0.std)
  M <- cor.balance(input, m=3, G=nrow(data_matrix)) 
  colnames(M) <- rownames(data_matrix)
  rownames(M) <- rownames(data_matrix)
  # rows_to_keep <- !(rownames(M) %in%Significant_proteins)
  # Correlated_protein <- M[rows_to_keep,Significant_proteins] 
  # Correlated_proteins <-  rownames(M)[matrixStats::rowSums2(Correlated_protein>0.65)>1]
  M_dist <- 1-M#[c(Correlated_proteins,Significant_proteins),c(Correlated_proteins,Significant_proteins)]
  # d <- as.dist(M_dist)
  # g<- diana(d)
  pheatmap(M)
  clusters <-  NbClust::NbClust(method = "complete", distance = "euclidean", data = M, index =  "silhouette")$Best.partition %>% 
    enframe(name = "gene", value = "cluster")
  fuzziness <- advclust::fuzzy.CM(M#[c(Correlated_proteins,Significant_proteins),c(Correlated_proteins,Significant_proteins)]
                        ,K=6 #clusters$cluster %>% max()
                        ,m=2,max.iteration=100,threshold=1e-5,RandomNumber=1234)
  Confident_items <-  fuzziness@member[(((fuzziness@member>0.4) %>% rowSums2()) <2) & (((fuzziness@member>0.8) %>% rowSums2()) ==1),] %>% 
    rownames()
  clusters <- fuzziness@hard.label[Confident_items]
  M_confi <-  M[Confident_items,Confident_items]
  row_ha = ComplexHeatmap::rowAnnotation(clusters = as.character(clusters))
  
  
  list(correp_matrix = M,
       correp_matrix_conf = M_confi ,
       clusters = clusters,
       heatmap = ComplexHeatmap::Heatmap(M_confi, left_annotation = row_ha,show_column_names = F,show_row_names = F, 
                                         name = "Protein Correlation", column_title = "Clustering of Significant Proteins"))
  
}


PeCoRa_function <- function(DIA_report_file,Samples, control_condition = "no Input"){
    #pacman::p_load(diann, PeCorA)
  # control_condition = NULL
      # DIA_report_file <- list_files_to_analyse$DIA_report_file_method2
     # Samples=Samples
    peptide_file <- read_tsv(here::here("Datasets","Raw",DIA_report_file)) 
    
    
    peptides.maxlfq <- diann_maxlfq(peptide_file[peptide_file$Lib.Q.Value <= 0.01 & peptide_file$Lib.PG.Q.Value <= 0.01,], 
                                    group.header="Stripped.Sequence",
                                    id.header = "Precursor.Id", 
                                    quantity.header = "Precursor.Normalised")
    PepCorA_setup <- peptides.maxlfq %>% as.data.frame()%>% rownames_to_column("Peptide") %>% 
        pivot_longer(-Peptide, names_to = "Condition_rep", values_to = "Normalised.Area") %>% 
         mutate(run = str_match(Condition_rep,"-([:graph:]{6})-[:graph:].raw$")[,2] %>% tolower()) %>% 
        
        left_join(Samples %>% mutate(run = tolower(run))) %>%
        separate(Condition, into = c("Condition","BioReplicate"), sep = "_") %>% 
        left_join(peptide_file %>% as.data.frame() %>% dplyr::select(Modified.Sequence,Stripped.Sequence,Protein.Group) %>% distinct(),
                  by = c("Peptide" = "Stripped.Sequence")) %>% as.data.frame() %>% mutate(Normalised.Area = replace_na(Normalised.Area,0)) %>% 
        mutate(BioReplicate = as.numeric(BioReplicate),
               Condition = as.factor(Condition)) %>% 
        dplyr::rename(Protein=Protein.Group,
                      Peptide.Modified.Sequence = Modified.Sequence)
    area_column_name <- which(str_detect(PepCorA_setup %>% colnames(),"ormali"))
    if(control_condition == "no Input"){
        control_condition <-PepCorA_setup$Condition[1] %>% as.character()
    }
    scaled_peptides <- PeCorA_preprocessing(PepCorA_setup,
                                            area_column_name = area_column_name,
                                            threshold_to_filter = 100,
                                            control_name = control_condition)
    disagree_peptides <- PeCorA::PeCorA(scaled_peptides)
    list(disagree_peptides = disagree_peptides,
         scaled_peptides = scaled_peptides,
         disagree_proteins = disagree_peptides %>% left_join(peptide_file %>% as.data.frame() %>% dplyr::select(Genes,Protein.Group) %>% distinct(),
                                                             by = c("protein" = "Protein.Group")))
    
    
    
}

Output_proteomic_ruler <- function(Samples_df,dataset_name,report_file, Uniprot_length_Mass){
 # report_file <- list_files_to_analyse$DIA_report_file_method2_w_old
 # dataset_name <- "testing"
  # Samples_df <- Samples
 peptide_file <- read_tsv(here::here("Datasets","Raw",report_file)) %>% 
   subset(Lib.Q.Value <= 0.01 & Lib.PG.Q.Value <= 0.01) %>%group_by(Protein.Group) %>% 
   janitor::clean_names()
ruler <- left_join(peptide_file%>% 
   dplyr::select(matches("pg_max_lfq|run|protein_group")) %>% 
   distinct() %>% 
   subset(!str_detect(protein_group,paste0(contaminant_uniprots,collapse = "|"))) %>% 
   mutate(run = str_replace_all(run, c("P11690-1" = "P11591-1",
                                       "P11690-2" = "P11592-1",
                                       "P11690-3" = "P11593-1")),
          run = str_match(run,"-([:graph:]{6})-")[,2] %>% tolower()) %>%  
   left_join(Samples_df) %>%  
   pivot_wider(-run, names_from = Condition ,values_from = pg_max_lfq),
   peptide_file %>% group_by(protein_group) %>% dplyr::summarise(Unique_razor_pept = sum(proteotypic),
             Total_pept = n()) ) %>% 
    # rename(Uniprot = protein_group) %>%  
  mutate(Uniprot = str_remove_all(protein_group,";[:graph:]*$")) %>% 
  left_join(Uniprot_length_Mass, by = "Uniprot") 
write_tsv(ruler,here::here("Datasets","Processed",glue::glue("For_Proteomic_ruler_",dataset_name,".txt")))

pept_condition <-  peptide_file %>% mutate(run = str_replace_all(run, c("P11690-1" = "P11591-1",
                                                                        "P11690-2" = "P11592-1",
                                                                        "P11690-3" = "P11593-1")),
                                           run = str_match(run,"-([:graph:]{6})-")[,2] %>% tolower()) %>%  
                     left_join(Samples_df) %>% separate(Condition, c("condition", "replicate")) %>%  group_by(protein_group,condition) %>% dplyr::summarise(Unique_razor_pept = sum(proteotypic),
                                                                          Total_pept = n())  %>% 
  # rename(Uniprot = protein_group) %>%  
  mutate(Uniprot = str_remove_all(protein_group,";[:graph:]*$")) %>% 
  left_join(Uniprot_length_Mass, by = "Uniprot") 

 
  write_tsv(pept_condition,here::here("Datasets","Processed",glue::glue("peptide_condition",dataset_name,".txt")))
  
}

Import_proteomic_ruler <- function(prot_ruler_file){
  # prot_ruler_file <- "From_proteomic_ruler.txt"
  input_matrix <- here::here("Datasets","Processed",prot_ruler_file) %>% 
    read.delim() %>% .[-1,] %>% 
    janitor::clean_names() %>% 
    subset(absolute_quantification_accuracy != "low") %>% 
    remove_rownames() %>% 
    dplyr::select(protein_group, contains("copy")) %>% as.data.frame() %>% 
    column_to_rownames("protein_group") %>%
    mutate(across(everything(), ~as.numeric(.x))) %>% 
    set_names(.,str_remove_all(colnames(.),"copy_number_")) %>% 
  mutate(across(where(is.numeric),~log2(as.numeric(.x))),
         across(where(is.numeric), ~if_else(is.infinite(.x), NaN,.x))) %>% 
    as.matrix() %>% 
    proDA::median_normalization()
  # experimental_design_DIA <-  data.frame(
  #   label = colnames(input_matrix),
  #   condition =  str_remove_all(colnames(input_matrix),"_[:graph:]*$"),
  #   replicate = str_remove_all(colnames(input_matrix),"^[:graph:]*_") %>% as.numeric()
  # )
  # data_unique_Etop <- input_matrix %>% rownames_to_column("name") %>% 
  #   left_join(HUMAN_9606 %>% 
  #               subset(Type == "Gene_Name") %>% 
  #               dplyr::select(-Type), 
  #             by = c("name" = "Uniprot")) %>%
  #   subset(!duplicated(name))
  # 
  # Quant_columns <- which(colnames(data_unique_Etop) %in%colnames(input_matrix))# get LFQ column numbers
  # data_se <- make_se(data_unique_Etop, Quant_columns, experimental_design_DIA)
  # data_se_parsed <- make_se_parse(data_unique_Etop, Quant_columns)
  # data_filt <- filter_missval(data_se, thr = 1)
  # data_filt@assays@data@listData[[1]][is.nan(data_filt@assays@data@listData[[1]])] <- NA 
  # plot_missval(data_filt)
  # data_norm <- normalize_vsn(data_filt)
  # data_norm@assays@data@listData[[1]]
  
}

Jaccard_index_list<-function(list_of_vectors, max_jacc = 0.5,steps = 1){

  #removes top step most similar with other pathways and most similar until max_jacc is reached
  #samples as cols
  # list_of_vectors = terms_to_reduce
  row_max = 1
  # top_similar = 1
  # max_jacc = 0.3
  # steps = 1
  while(row_max>max_jacc){
    # list_of_vectors <- enrichment_df_BP %>% subset(`p.adjust`<0.05 ) %>% 
    #   dplyr::select(core_enrichment,Description) %>% subset(!duplicated(Description)) %>% 
    #   pull(core_enrichment,Description)  %>% purrr::map(str_split,"/") %>% flatten()
    index<-map(.x = list_of_vectors, ~.x %>% list(.) %>% rep(.,length(list_of_vectors)) %>% 
                 map2_dbl(.x = .,.y = list_of_vectors,~bayesbio::jaccardSets(.x,.y))) %>% 
      imap_dfr(.x = ., ~set_names(.x,names(list_of_vectors)) %>% enframe(name = "Pathway2",value = "JaccIndex") %>% mutate(Pathway1 = .y)) %>% 
      subset(JaccIndex != 1) %>% as.data.table() 
    setkey(index,Pathway1,Pathway2)
    index <- index[Pathway1>Pathway2]
    row_max <- index$JaccIndex %>% max()
    index <- index[JaccIndex == row_max][,top_similar:= fifelse(str_detect(Pathway1,"nucleo|eactive"),Pathway2,Pathway1)]
    
    top_similar <- index$top_similar
      # pivot_wider(names_from = "Pathway2", values_from = "JaccIndex") %>% 
      # column_to_rownames("Pathway1") %>% as.matrix()
    # diag(index) <- 0
    # row_max <- index %>% matrixStats::rowMaxs() %>% set_names(row.names(index)) %>% sort(decreasing = T) %>% .[1]
    #top_similar <- index %>% matrixStats::rowSums2() %>% set_names(row.names(index)) %>% sort(decreasing = T) %>% names() %>% .[1:steps]
    # top_similar <- c(top_similar)
    list_of_vectors[which(names(list_of_vectors)%in%top_similar)]<-NULL
    print(row_max)
  }
  list_of_vectors
}
