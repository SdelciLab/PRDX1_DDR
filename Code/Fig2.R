#### Etop DIA ####
pacman::p_load(here, tidyverse, DEP,pheatmap,diann,sva,imp4p,proDA,
org.Hs.eg.db,clusterProfiler,ggridges,SubCellBarCode,eulerr,scales,data.table,
visNetwork,matrixStats,magick,openxlsx, janitor,seqinr,UpSetR)

#Change the inference to Protein level
#Double run through NN
#### Loading####

source(here::here("Code","functions.R"))

#genes involved in human metabolism downloaded from
# https://github.com/SysBioChalmers/Human-GEM/blob/main/model/genes.tsv
Human_GEM_genes <- here::here("Datasets","Raw","Human_GEM_model_genes.tsv") %>% read_tsv() %>%  
pull(geneSymbols) %>% unique()

# Proteomic Ruler of WCE from CCLE
# add citation 
Proteomic_Ruler <- here::here("Datasets","Processed","CCLE_prot_Ruler.txt") %>% read.delim() %>% .[-1,] %>% 
dplyr::select(matches("Copy|Uniprot_Acc|accuracy"))%>% 
remove_rownames() %>% 
column_to_rownames("Uniprot_Acc") %>% 
#mutate(across(where(is.numeric), as.numeric)) %>% 
set_names(.,str_remove_all(names(.),"Copy\\.number\\.")) %>% 
mutate(across(contains("_"),~log2(as.numeric(.x))),
across(where(is.numeric), ~if_else(is.infinite(.x), NaN,.x))) %>% 
subset(!is.nan(U2OS_BONE)) %>% 
# subset(.,rowSums(is.na(.))<(ncol(.)/3)) %>%
subset(Absolute.quantification.accuracy != "low") %>%
dplyr::select(-Absolute.quantification.accuracy) %>%
janitor::clean_names()



Uniprot_length_Mass <- here::here("Datasets","Raw", "Uniprot_Molecular_sizes.tab") %>% 
read_tsv() %>% 
dplyr::select(-`Entry name`) %>% 
set_names(c("Uniprot","Length","Mass")) %>% 
mutate(Mass = Mass/1000)

HUMAN_9606 <- read_tsv(here::here("Datasets","Raw", "HUMAN_9606_idmapping.dat"),
col_names = FALSE) %>% set_names(c("Uniprot","Type","ID"))

Human_hsa <- HUMAN_9606 %>% 
subset(Type == "KEGG") %>% 
mutate(ID = str_remove_all(ID,"hsa:")) %>% 
dplyr::select(-Type) %>% 
subset(!duplicated(Uniprot))

# Sabatini library from AddGene
Sabatini <- openxlsx::read.xlsx(here::here("Datasets","Raw","metabolism_gene_list_Sabatini.xlsx")) %>% 
pull(Gene.Symbol)
Sabatini_Uniprot <- HUMAN_9606 %>% subset((ID %in% Sabatini) & Type == "Gene_Name") %>%
subset(!duplicated(Uniprot)) %>% pull(Uniprot, ID)

# Hits_Joanna <- openxlsx::read.xlsx(here::here("Datasets","Raw",
#               "Candidates_proteins.xlsx")) %>% pull(`71.genes`) %>% c("FDX2")
# Candidates_Joanna<-openxlsx::read.xlsx(here::here("Datasets","Raw",
#                   "Candidates_proteins.xlsx"), sheet =2 , colNames = F)[,1] %>% unlist() %>%c("FDX2")
# Candidates_Joanna <- c(All_interesting_genes, Screens_new$Etop_High_low_10_cnv_negsg.mle.gene_summary_fixed_names.txt %>%
# subset(`high|fdr`<0.05 & abs(`high|beta`)>0.2) %>% pull(Gene)) %>% unique()
# early_fast <- c("Ku80","PRKDC","XRCC4","PARP1","Ku70","LIG4",
# "PCNA","RAD50","ATM")
# early_slow <- c("MDC1",'53BP1',"MRE11","BRCA1","BRCA2","RPA","BARD1","PAXIP1")
# late_slow <- c("RAD54",'RAD52',"RAD51")
# Proteins_known_behaviour<- data.frame(ID = c(early_fast,early_slow,late_slow),
#       Behaviour = rep(c("early_fast","early_slow","late_slow"), c(length(early_fast),
#                                                                   length(early_slow),
#                                                                   length(late_slow))))
# Uniprot_known_behaviour <- left_join(Proteins_known_behaviour, HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type)) %>% 
# group_split(Behaviour) %>% 
# set_names(.,map_chr(.x = .,~.x %>% pull(Behaviour) %>% unique())) %>% 
# map(~pull(.x,Uniprot) %>% unique() %>% na.omit())
# KORE<- openxlsx::read.xlsx(here::here("Datasets","Raw","KORE_DDR_genes_uniprot_IDS.xlsx"))[,2]
# 
# DDR_GO<- subset(HUMAN_9606, ID %in% (read_csv(here::here("Datasets","Raw","QuickGO-annotations-DSB repair_without duplicates.csv")) %>%
#          pull(SYMBOL))) %>% pull(Uniprot) %>% unique()
list_of_metabolism = list("Metabolic_library" = Sabatini_Uniprot,
                          "Human_GEM" = HUMAN_9606 %>% subset(ID %in% Human_GEM_genes) %>% 
                              pull(Uniprot) %>% unique())

# lists_of_interest = list("Early_Fast" = Uniprot_known_behaviour$early_fast,
# "Early_Slow"=  Uniprot_known_behaviour$early_slow,
# "Late_slow" =Uniprot_known_behaviour$late_slow,
# "KORE" = KORE,
# "DDR_GO" = DDR_GO,
# "Hits_Joanna" = subset(HUMAN_9606, ID %in% (Hits_Joanna)) %>% pull(Uniprot) %>% unique(),
# "Candidates_Joanna" = subset(HUMAN_9606, ID %in% (Candidates_Joanna)) %>% pull(Uniprot) %>% unique() )

add_multiple_behaviour <- function(list_Uniprot,behaviour){
# list_Uniprot <- list_of_metabolism[[1]]
# behaviour <- names(list_of_metabolism  )[1]
Interesting_proteins <- data.frame(Uniprot= c(list_Uniprot %>% unlist(use.names = F))) %>% 
mutate(Behaviour = " ") %>%
mutate(Behaviour = if_else(Uniprot %in% list_Uniprot, paste(Behaviour, behaviour, sep = ";"),Behaviour))
Interesting_proteins
}
Metabolic_proteins <- purrr::imap_dfr(list_of_metabolism,add_multiple_behaviour) %>%
    distinct() %>% 
    group_by(Uniprot) %>% 
    mutate(Behaviour = paste0(Behaviour, collapse = ""),
           Behaviour = str_remove_all(Behaviour,"^ \\;")) %>% 
    distinct()

# Interesting_proteins <- purrr::imap_dfr(lists_of_interest,add_multiple_behaviour) %>%
# distinct() %>% 
# group_by(Uniprot) %>% 
# mutate(Behaviour = paste0(Behaviour, collapse = ""),
# Behaviour = str_remove_all(Behaviour,"^ \\;")) %>% 
# distinct()
Compartment_colours = c(Chromatin = "darkred",
                        Cytoplasm = "darkblue",
                        Mitochondria = "grey50",
                        Nuclear = "grey50",
                        Secretory = "grey50")
# colours_compartmens <- c("#ADE5F0","#2FFEE0","#29FEFD","#34D5F7","#00D8DF","#C9B3A4","#C58F6C",
#   "#ECECEC","#D1D1D1","#A8A8A8","#6B6B6B",
#   "#FFDD17","#FFBF3D","#FFA095","#FF846D") %>% set_names(markerProteins$Compartments %>% levels()) 
  

  # Complexes <- read.delim(here::here("Datasets","Raw","coreComplexes.txt")) %>% subset(Organism == "Human") %>% 
  #     janitor::clean_names() %>% 
  #     pull(subunits_uni_prot_i_ds,complex_name)
  # Complexes <- Complexes %>% str_split(";") %>% set_names(names(Complexes))
  #####
  list_files_to_analyse <- list(
  # DIA_report_file_1 = "
# ",
                                # DIA_report_file_method2 = "P11833_P11841_Method_2_report.tsv",
                                # DIA_report_file_method2_w_old = "Method2_original_report.tsv",
                                DIA_report_file_method2_old_better = "Method2_original_report_new.tsv")
                                # DIA_report_file_P11685 = "P11685_Discovery_report.tsv")
  
  Samples <- data.frame(run=  c("P11833","P11834","P11835","P11836","P11837" ,"P11838","P11839","P11840" ,"P11841",
                                "p11591","p11592","p11593","p11685","p11688" ,"p11691","p11686","p11689" ,"p11692"),
                        Condition = c("DMSO_1","DMSO_2","DMSO_3","0H_1", "0H_2","0H_3","24H_1", "24H_2" ,"24H_3",
                                      "DMSO_4","DMSO_5","DMSO_6","0H_4", "0H_5","0H_6","24H_4", "24H_5" ,"24H_6"))%>% 
      dplyr::mutate(dplyr::across(tidyselect::everything(),janitor::make_clean_names))
  
  # Combined_batch <- map(.x = list_files_to_analyse[[3]], 
  #                  ~Load_DIA_NN_Data(.x,Samples))
  
  # MaxQuant contaminat files
  # used in the loading function to remove them
  contaminant_uniprots <- seqinr::read.fasta(here::here(
      "Datasets","Raw","contaminants.fasta"),seqtype ="AA") %>% names() 
  
  Methods_DIA <- map(.x = list_files_to_analyse[1],
                     ~Load_DIA_NN_Data(.x,Samples))
  # Methods_1_report <- read_tsv(here::here("Datasets", "Raw",list_files_to_analyse$DIA_report_file_1))
  # Methods_1_report %>% 
  #   janitor::clean_names() %>% 
  #   subset(str_detect(run,"P11833")) %>% 
  #   group_by(protein_names) %>% 
  #   glimpse() %>% 
  #   summarise(Amount_of_evidence = n()) %>% 
  #   mutate(Shared = if_else(str_detect(protein_names,";"),"Shared","Unique")) %>% 
  #  # View()
  #   ggplot(aes(x = Shared, y = log2(Amount_of_evidence), colour = Shared, label =protein_names )) +
  #   geom_boxplot()+
  #   geom_point()+
  #   ggrepel::geom_label_repel(data = . %>% subset(log2(Amount_of_evidence)>7.5))
  # Methods_DIA$DIA_report_file_1 %>% row.names() %>% str_detect(";") %>% table()
  # plot(eulerr::euler(list("Unique_groups" =Methods_DIA$DIA_report_file_1 %>% row.names() %>% str_subset(";", negate = T),
  #                          "Shared_groups" = Methods_DIA$DIA_report_file_1 %>% row.names() %>% str_subset(";") %>% 
  #                            str_split(";",simplify = T) %>% .[.%in% (Methods_DIA$DIA_report_file_1 %>% row.names() %>% str_subset(";", negate = T))]  %>% unique()),
  #                     shape = "ellipse"), quantities = TRUE, main = "Uniprots in Unique and Shared Protein entries from DIANN")
  # 
  
  Methods_DIA_DEP$DIA_report_file_method2_old_better[-3] %>% write.xlsx(overwrite = T, file = here::here("Datasets","Processed","Batch_corrected_matrices.xlsx"), rowNames = T)
  Methods_DIA_DEP$DIA_report_file_method2_old_better$DEPs %>% write.xlsx(overwrite = T, file = here::here("Datasets","Processed","Batch_corrected_DEPs.xlsx"),rowNames = T)
  
  #testing <- Load_DIA_NN_Data(list_files_to_analyse$DIA_report_file_method2_w_old,Samples)
  pept_condition <- map(.x = list_files_to_analyse, ~
                            read_tsv(here::here("Datasets","Processed",glue::glue("peptide_condition",.x,".txt")))) 
  Proteomic_Ruler_input <- Import_proteomic_ruler("From_proteomic_ruler.txt")
  Temp_proteomic_ruler_sample <-Proteomic_Ruler_input %>%
      as.data.frame() %>% 
      rownames_to_column("PG") %>%
      mutate(Uniprot = str_remove_all(PG, ";[:graph:]*$"))
  Temp_proteomic_ruler_sample_2 <- Temp_proteomic_ruler_sample %>% 
      left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% 
                    subset(!duplicated(Uniprot)), by =  "Uniprot")%>%
      left_join(markerProteins[,1:2], by = c("ID" = "Proteins")) %>% 
      pivot_longer(cols = contains("_"), names_to = "Condition", values_to = "Abundance")
  ggplot(Temp_proteomic_ruler_sample_2 %>% subset(Condition == "dmso_1"),
         aes(x = reorder(PG, Abundance), y= Abundance, colour = Compartments))+
      geom_col()+
      facet_wrap("Compartments")+
      ggtitle("Copy Numbers on chrom DMSO_1")
  
  Temp_proteomic_ruler <-     left_join(Temp_proteomic_ruler_sample,Proteomic_Ruler %>% rownames_to_column("Uniprot") %>% 
                                            dplyr::select(Uniprot,u2os_bone)) %>% 
      mutate(across(where(is.numeric), ~.x-u2os_bone)) %>% 
      subset(!is.na(u2os_bone))
  DMSO_ruler <- Temp_proteomic_ruler %>% 
      left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% 
                    subset(!duplicated(Uniprot)), by =  "Uniprot") %>% 
      
      dplyr::select(matches("dmso|ID"))
  DMSO_ruler$Mean_DMSO = DMSO_ruler %>% dplyr::select(-ID) %>%as.matrix() %>%   rowMeans2 (.,na.rm = T) %>% scale()
  DMSO_ruler <- dplyr::select(DMSO_ruler, ID, Mean_DMSO)
  # Missing_U2OS_test <- list(Missing_U2OS = Temp_proteomic_ruler %>% subset(is.na(u2os_bone)) %>% pull(Uniprot) %>% unique(),
  #                           # Missing_Etop_DIA =  Temp_proteomic_ruler %>% subset(!is.na(u2os_bone)) %>% pull(Uniprot) %>% unique()
  #                      Universe = Temp_proteomic_ruler %>% pull(Uniprot) %>% unique())
  # 
  # ego <- enrichGO(gene          = Missing_U2OS_test$Missing_U2OS,
  #                 universe      = Missing_U2OS_test$Universe,
  #                 OrgDb         = org.Hs.eg.db,
  #                 ont           = "CC",
  #                 keyType = "UNIPROT",
  #                 pAdjustMethod = "BH",
  #                 pvalueCutoff  = 0.05,
  #                 qvalueCutoff  = 0.05,
  #                 readable      = TRUE)
  # #enrichplot::upsetplot(ego %>% simplify(), n  = 20)
  # dotplot(ego, showCategory=30) + ggtitle("Proteins Found by Chromatin DIA but not U2OS CCLE_TMT")
  Proteomic_Ruler %>% 
      dplyr::select(u2os_bone ) %>% 
      rownames_to_column("Uniprot") %>% 
      na.omit() %>% 
      left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% 
                    subset(!duplicated(Uniprot)), by =  "Uniprot")%>%
      left_join(markerProteins[,1:2], by = c("ID" = "Proteins")) %>% 
      # pivot_longer(cols = where(is.numeric), names_to = "Condition", values_to = "Abundance") %>% 
      # subset(absolute_quantification_accuracy != "low") %>% 
      ggplot( aes(x = reorder(Uniprot, u2os_bone), y= u2os_bone, colour = Compartments))+
      geom_col()+
      facet_wrap("Compartments")+
      ggtitle("WCE copy numbers u2os")
  
  
  Temp_proteomic_ruler <- Temp_proteomic_ruler %>% 
      subset(!is.na(u2os_bone)) %>% 
      left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% 
                    subset(!duplicated(Uniprot)), by =  "Uniprot")%>%
      left_join(markerProteins[,1:2], by = c("ID" = "Proteins")) %>% 
      pivot_longer(cols = where(is.numeric), names_to = "Condition", values_to = "Abundance") %>% 
      subset(Condition != "u2os_bone")
  for(i in 1:3){
      # i = 1
      ggplot(Temp_proteomic_ruler %>% subset(Condition == glue::glue("dmso_{i}")) %>% 
                 mutate(Compartment = case_when(
                     Compartments %in% c("N1","N2","N4") ~"Nuclear",
                     Compartments %in% c("N3") ~"Chromatin",
                     str_detect(Compartments, "^C") ~"Cytoplasm",
                     str_detect(Compartments, "^S") ~"Secretory",
                     str_detect(Compartments, "^M") ~"Mitochondria",
                 ),
                 Abundance = scale(Abundance)[,1]) %>% 
                 subset(!is.na(Compartment)),
             aes(x = Compartment, y= Abundance, colour = Compartment, label = ID))+
          geom_boxplot()+
          geom_point(data = . %>% subset(ID == "COX4I1"))+
          ggrepel::geom_label_repel(data = . %>% subset(ID == "COX4I1"))+
          # geom_jitter(alpha= 0.1)+
          
          theme(axis.ticks.x = element_blank(),
                axis.text.x = element_blank(), legend.position = "none")+
          # facet_wrap("Compartment")+ 
          labs(x = "Proteins with known localisation", y = "Relative Enrichemnt" )+ theme_bw()+
          scale_colour_manual(values = Compartment_colours)+
          ggtitle("Protocol Results Enrichement of Chromatin Protein and Depletion of Cytoplasmic Proteins")
      ggsave(here::here("Output","Figures",glue::glue("Chromatin_enrichment_DMSO_{i}.pdf")))
  }
  ggplot(Temp_proteomic_ruler %>% subset(str_detect(Condition,"dmso")) %>% 
             group_by(Uniprot, ID, Compartments) %>% summarise(Abundance = mean(Abundance,na.rm = T))%>% ungroup %>% 
             mutate(Compartment = case_when(
                 Compartments %in% c("N1","N2","N4") ~"Nuclear",
                 Compartments %in% c("N3") ~"Chromatin",
                 str_detect(Compartments, "^C") ~"Cytoplasm",
                 str_detect(Compartments, "^S") ~"Secretory",
                 str_detect(Compartments, "^M") ~"Mitochondria",
             ),
             Abundance = scale(Abundance)[,1]) %>% 
             subset(!is.na(Compartment)),
         aes(x = Compartment, y= Abundance, colour = Compartment, label = ID))+
      geom_boxplot()+
      geom_point(data = . %>% subset(ID %in% c("BRD4","EEF2","COX4I1","H2AFY","DDX5","H1FX","TOP2B","KDM2A","NDUFS6","NDUFA8","UQCR10","SLC7A5")))+
      ggrepel::geom_label_repel(data = . %>% subset(ID %in% c("BRD4","EEF2","COX4I1","H2AFY","DDX5","H1FX","TOP2B","KDM2A","NDUFS6","NDUFA8","UQCR10","SLC7A5")))+
      # geom_jitter(alpha= 0.1)+
      
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank(), legend.position = "none")+
      # facet_wrap("Compartment")+ 
      labs(x = "Proteins with known localisation", y = "Relative Enrichemnt" )+ theme_bw()+
      scale_colour_manual(values = Compartment_colours)+
      ggtitle("Protocol Results Enrichement of Chromatin Protein and Depletion of Cytoplasmic Proteins")
  ggsave(here::here("Output","Figures",glue::glue("Chromatin_enrichment_DMSO_avg.pdf")))
  
  dmso_gsea <- gseGO(geneList     = Temp_proteomic_ruler %>% subset(Condition == "dmso_1" & !is.na(Abundance) & !duplicated(Uniprot)) %>% arrange(-Abundance) %>%
                         pull(Abundance,Uniprot) %>%  scale() %>% .[,1] ,
                     OrgDb        = org.Hs.eg.db,
                     ont          = "CC",
                     minGSSize    = 100,
                     keyType = "UNIPROT",
                     maxGSSize    = 500,
                     pvalueCutoff = 0.05,
                     verbose      = FALSE)
  enrichplot::ridgeplot(dmso_gsea , showCategory = 68)+ggtitle("enrichment of GO_CC in dmso_1 normalised to U2OS cells ")+
      scale_fill_gradient2(high = scales::muted("red"),mid = scales::muted("blue"))
  ggsave(here::here("Output","Figures","Chromatin_enrichment_GSEA_DMSO.pdf"))
  #enrichplot::upsetplot(ego %>% simplify(), n  = 20)
  dmso_gsea_sample <- gseGO(geneList     = Temp_proteomic_ruler_sample_2 %>% subset( Condition == "dmso_1" & !is.na(Abundance) & !duplicated(Uniprot)) %>% arrange(-Abundance) %>%
                                pull(Abundance,Uniprot) %>%  scale() %>% .[,1] ,
                            OrgDb        = org.Hs.eg.db,
                            ont          = "CC",
                            minGSSize    = 100,
                            keyType = "UNIPROT",
                            maxGSSize    = 500,
                            pvalueCutoff = 0.05,
                            verbose      = FALSE)
  ridgeplot(dmso_gsea_sample %>% simplify, showCategory = 68)+ggtitle("enrichment of GO_CC in dmso_1 unnormalised to U2OS cells ")
  
  
  
  Temp_proteomic_ruler %>%  
      group_by(Uniprot,ID,Compartments ) %>%
      dplyr::summarise(Mean_Abundance_u2os = mean(Abundance, na.rm = T),.groups = "keep") %>% 
      ungroup() %>% 
      left_join(Temp_proteomic_ruler_sample_2 %>% 
                    group_by(Uniprot,ID,Compartments ) %>%
                    dplyr::summarise(Mean_Abundance_sample = mean(Abundance, na.rm = T),.groups = "keep") %>% 
                    ungroup()) %>% subset(!is.na(Compartments)) %>% 
      ggplot(aes (x = reorder(Uniprot,Mean_Abundance_sample),y = Mean_Abundance_sample, fill = Compartments))+
      geom_col(width = 1)+  facet_grid(~Compartments, scales = "free_x")+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      scale_fill_manual(values = colours_compartmens)+
      labs(x= "Proteins Identified on Chromatin",y ="Copy Numbers on Chromatin" )+
      ggtitle("Protein Copy Numbers Detected on Chromatin")
  
  
  
  
  
  Temp_proteomic_ruler %>%  
      group_by(Uniprot,ID,Compartments ) %>%
      dplyr::summarise(Mean_Abundance_u2os = mean(Abundance, na.rm = T),.groups = "keep") %>% 
      ungroup() %>% 
      left_join(Temp_proteomic_ruler_sample_2 %>% 
                    group_by(Uniprot,ID,Compartments ) %>%
                    dplyr::summarise(Mean_Abundance_sample = mean(Abundance, na.rm = T),.groups = "keep") %>% 
                    ungroup()) %>% 
      mutate(across(contains("Abundance"), scale)) %>% 
      # subset(is.na(Compartments)) %>% 
      ggplot(aes(x = Mean_Abundance_u2os, y = Mean_Abundance_sample, label = ID, colour = Compartments))+
      geom_point()+
      ggrepel::geom_text_repel(data = . %>% subset(str_detect(ID, "^RPL")),max.overlaps = 50 )+
      geom_abline(slope=1, intercept = 0)+
      facet_wrap("Compartments")+
      ggtitle("% of Protein on chromatin vs Copy numbers on chromatin",
              "Using proteomic Ruler normalised to U2OS CCLE_TMT")
  Temp_proteomic_ruler %>%  
      group_by(Uniprot,ID,Compartments ) %>%
      dplyr::summarise(Mean_Abundance_u2os = mean(Abundance, na.rm = T),.groups = "keep") %>% 
      ungroup() %>% 
      left_join(Temp_proteomic_ruler_sample_2 %>% 
                    group_by(Uniprot,ID,Compartments ) %>%
                    dplyr::summarise(Mean_Abundance_sample = mean(Abundance, na.rm = T),.groups = "keep") %>% 
                    ungroup()) %>% 
      mutate(across(contains("Abundance"), scale)) %>% 
      subset(Compartments == "N1") %>% 
      ggplot(aes(x = Mean_Abundance_u2os, y = Mean_Abundance_sample, label = ID))+
      geom_point()+
      ggrepel::geom_text_repel(data = . %>% subset(str_detect(ID, "^RPL")),max.overlaps = 50 )+
      geom_abline(slope=1, intercept = 0)+
      facet_wrap("Compartments")+
      ggtitle("% of Protein on chromatin vs Copy numbers on chromatin",
              "Using proteomic Ruler normalised to U2OS CCLE_TMT")
  
  
  Output_proteomic_ruler(data_norm@assays@data@listData[[1]] %>% as.data.frame(), "testing_temp", Uniprot_length_Mass )
  DIA_matrices <- map(.x = list_files_to_analyse, 
                      ~Load_DIA_NN_Data(.x,Samples%>% 
                                            mutate(across(everything(),janitor::make_clean_names))))
  DIA_DEP <- map2(DIA_matrices,list_files_to_analyse,DEP_DIA)
  
  
  Methods_DIA_DEP$DIA_report_file_P11685$DEPs$x0h_vs_x24h_diff[,c("log2_FC","Uniprot", 'significant')] %>%
      dplyr::rename(P11685_Log2_FC_24_vs_0 = log2_FC,
                    significant_P11685 = significant) %>% 
      inner_join( Methods_DIA_DEP$DIA_report_file_method2$DEPs$x0h_vs_x24h_diff[,c("log2_FC","Uniprot", 'significant')] %>% 
                      dplyr::rename(Method2_Log2_FC_24_vs_0 = log2_FC,
                                    significant_method2 = significant)) %>% 
      mutate(significant = case_when(
          significant_method2 == T  & significant_P11685==T ~ "Both",
          significant_method2 == T  & significant_P11685==F ~ "method2",
          significant_method2 == F  & significant_P11685==T ~ "P11685",
          significant_method2 == F  & significant_P11685==F ~ "neither"
          
      )) %>% 
      ggplot(aes(x = P11685_Log2_FC_24_vs_0, y =  Method2_Log2_FC_24_vs_0,colour =  significant ))+
      geom_point(alpha = 0.1, data = . %>% subset(significant == "neither"))+
      geom_point(alpha = 1, data = . %>% subset(significant != "neither"))+
      geom_point(alpha = 1, size = 30, data = . %>% subset(significant == "both"))+
      
      theme_bw()+
      ggtitle("Correlation of Fold change between P11685 Samples and method 2 P11833 samples",
              subtitle = "24hrs vs 0 hrs chromatin")
  
  Method1_P11833_41_DEP$DEPs$x0h_vs_x24h_diff[,c("log2_FC","Uniprot")] %>%
      dplyr::rename(Method1_Log2_FC_24_vs_0 = log2_FC) %>% 
      inner_join(Method2_P11833_41_DEP$DEPs$x0h_vs_x24h_diff[,c("log2_FC","Uniprot")] %>% 
                     dplyr::rename(Method2_Log2_FC_24_vs_0 = log2_FC)) %>% 
      ggplot(aes(x = P11685_Log2_FC_24_vs_0, y =  Method2_Log2_FC_24_vs_0))+
      geom_point()
  
  
  Pecora_result <- PeCoRa_function(list_files_to_analyse$DIA_report_file_method2, Samples)
  
  ego <- enrichGO(gene          = Pecora_result$disagree_peptides %>% subset(adj_pval<0.05) %>% pull(protein) %>% 
                      str_remove_all(";[:graph:]*$"),
                  universe      = Pecora_result$disagree_peptides  %>% pull(protein) %>% 
                      str_remove_all(";[:graph:]*$"),
                  OrgDb         = org.Hs.eg.db,
                  keyType = "UNIPROT",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = FALSE)
  barplot(ego, showCategory = 100)+ggtitle("Protein Groups with Significant PeCoRA Peptides")
  
  PeCorA::PeCorA_plotting(Pecora_result$disagree_peptides,Pecora_result$disagree_peptides[1,],Pecora_result$scaled_peptides)
  row.names(testing) <- testing$peptide
  ego <- enrichGO(gene          = testing %>% subset(adj_pval<0.01) %>% pull(protein) %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-.") %>% unique(),
                  universe      = testing %>% pull(protein) %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-.")%>% unique(),
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  keyType = "UNIPROT",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  #enrichplot::upsetplot(ego %>% simplify(), n  = 20)
  dotplot(ego%>% simplify(), showCategory=30) + ggtitle("dotplot for ORA")
  hist(log2(df_input[, "Precursor.Normalised"]), 100, main = "Histogram of log2 intensities", 
       col = "steelblue", border = "steelblue", freq = FALSE)
  
  
  
  plot(eulerr::euler(list("Human_GEM" = Human_GEM_genes,
                          "Metabolic_Screen" = Sabatini),
                     shape = "ellipse"), quantities = TRUE)
  ego <- enrichGO(gene          = setdiff(Human_GEM_genes,Sabatini),
                  universe      = unique(c(Human_GEM_genes,Sabatini)),
                  OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                  keyType = "SYMBOL",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = FALSE)
  barplot( clusterProfiler::simplify(ego),showCategory=20) 
  
  
  Data_matrices <- map(.x = 1:2,
                       ~openxlsx::read.xlsx(here::here("Datasets","Processed", "Batch_corrected_matrices.xlsx"), sheet = .x, rowNames = T)) %>% 
      set_names(getSheetNames(here::here("Datasets","Processed", "Batch_corrected_matrices.xlsx")))  
  natalia_enzymes <- Data_matrices$Unimputted[,1:3]
  natalia_enzymes$av = mean(c(natalia_enzymes$dmso_1,natalia_enzymes$dmso_2,natalia_enzymes$dmso_3),na.rm = T)
  natalia_enzymes_prot =  rownames(natalia_enzymes)[!is.na(natalia_enzymes$av)] %>% str_split(";") %>% unlist()
  Interesting_proteins_nat = c("CBSL","CBS","DNMT",'MTHFR','MAT1A','BHMT','PEMT',
                               "DHFR","SHMT1"  ,"SHMT2" , "MAT2B" , "MAT2A"  ,"MTRR" ,'CHDH',
                               "PEMT",   "AHCYL1" ,"AHCYL2", "AHCY"  , "TYMS")
  HUMAN_9606 %>% 
      subset(Type == "Gene_Name" & Uniprot %in% natalia_enzymes_prot) %>% subset(str_detect(ID, paste(Interesting_proteins_nat,collapse = "|" ))) %>% pull(ID)
  
  list_of_genes <- c("RAD18","BARD1","P4HA2","TPI1","PTPN14","TGM2","BLVRB","IMPDH2","ENO1","PLPP3","CTDSPL2")
  list_of_uniprots_w_genes <- HUMAN_9606 %>%
      subset(ID %in% list_of_genes & Type == "Gene_Name") %>% 
      dplyr::select(-Type) %>% 
      distinct 
  list_of_uniprots <-list_of_uniprots_w_genes %>% pull(Uniprot)
  Data_matrices$Imputted %>% 
      subset(rownames(.)=="Q06830") %>% as.data.frame() %>% 
      rownames_to_column("Uniprot") %>% 
      pivot_longer(-Uniprot,names_to = "Sample",values_to = "Abundance") %>% 
      mutate(Condition = str_remove_all(Sample,"_.$"),
             Replicate = str_remove_all(Sample,"[:graph:]*_"))%>% 
      ggplot(aes(x = Condition, y = Abundance, group = Uniprot))+
      geom_point(alpha = 0.2)+theme_bw()+
      annotate("text",label = "PRDX1", x="x24h", y=20.25 ,size = 10, colour = muted("red"))+
      geom_smooth(fill = muted("red"),
                  colour = muted("red"), alpha = 0.1)+
      scale_x_discrete(labels= c("DMSO","No Release","Release"))+
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=18,face="bold"))+
      ggtitle("PRDX1 Abundance across the conditions")
  ggsave(here::here("Output","Figures", 'PRDX1_abundace_proteomics_across_conditions.pdf') ) 
  Mito_complex_genes <- openxlsx::read.xlsx(here::here("Datasets","Processed","Mito_complex_genes_enrichment.xlsx"))
  Mito_complex_genes <- Mito_complex_genes %>% subset(!is.na(Complex))
  openxlsx::write.xlsx(Mito_complex_genes, here::here("Datasets","Processed","Enriched_ETC_proteomics_screen.xlsx"))
  eulerr::euler(list(Survival_screen = Mito_complex_genes %>% subset(comparison == 'D14_BP') %>% pull(Genes)  , 
                     Chromatin_24h_release = Mito_complex_genes %>% subset(comparison != 'D14_BP') %>% pull(Genes))) %>% plot()
  
  Mito_complex_genes_Significant_complex <- Mito_complex_genes %>% subset(comparison == "dmso_vs_x24h_diff" & !is.na(Complex))
  lines <- c("#832424","#c73e3e","#a30303","#730505") %>% set_names(Mito_complex_genes_Significant_complex$Complex %>% unique())
  Data_matrices$Imputted %>% 
      mutate(Rowmean = rowMeans(.),
             across(where(is.numeric),~.x- Rowmean)) %>% 
      dplyr::select(-Rowmean) %>% 
      as.data.frame() %>% 
      rownames_to_column("Uniprot") %>% 
      
      pivot_longer(-Uniprot,names_to = "Sample",values_to = "Abundance") %>% 
      left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% 
                    subset(!duplicated(Uniprot)), by =  "Uniprot")%>%
      mutate(Condition = str_remove_all(Sample,"_.$"),
             Replicate = str_remove_all(Sample,"[:graph:]*_"))%>%  
      
      group_by(Uniprot,ID,Condition) %>% 
      summarise(Abundance = mean(Abundance)) %>% 
      inner_join(Mito_complex_genes_Significant_complex, by = c("ID" = "Genes")) %>% 
      ggplot(aes(x = Condition, y = Abundance, group =Uniprot, colour= Complex))+
      geom_point(alpha = 0.1)+
      geom_line(alpha = 0.1)+
      theme_bw()+
      # annotate("text", x="x24h", y=20.25 ,size = 10, colour = muted("red"))+
      geom_smooth(alpha = 0.05, aes(group = Complex, fill = Complex,colour = Complex))+
      scale_x_discrete(labels= c("DMSO","No Release","Release"))+
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=18,face="bold"))+
      ggtitle("Mitocomplexes Abundance across the conditions")+
      scale_fill_manual(values = lines)+
      scale_colour_manual(values = lines)
  ggsave(here::here("Output","Figures", 'mitocomplex_abundace_proteomics_across_conditions.pdf'),width = 10,height = 10 ) 
  
  Data_matrices$Imputted %>% 
      subset(rownames(.)%in%list_of_uniprots) %>% as.data.frame() %>% 
      rownames_to_column("Uniprot") %>% 
      pivot_longer(-Uniprot,names_to = "Sample",values_to = "Abundance") %>% 
      mutate(Condition = str_remove_all(Sample,"_.$"),
             Replicate = str_remove_all(Sample,"[:graph:]*_"))%>% 
      left_join(list_of_uniprots_w_genes) %>% 
      mutate(ID = factor(ID, levels = list_of_genes)) %>% 
      ggplot(aes(x = Condition, y = Abundance, group = Uniprot))+
      geom_point(alpha = 0.2)+theme_bw()+
      # annotate("text",label = "PRDX1", x="x24h", y=20.25 ,size = 10, colour = muted("red"))+
      geom_smooth(fill = muted("red"),
                  colour = muted("red"), alpha = 0.1)+
      scale_x_discrete(labels= c("DMSO","No Release","Release"))+
      theme(axis.text=element_text(size=12),strip.text = element_text(size=18),
            axis.title=element_text(size=12,face="bold"))+
      ggtitle("Abundance across the conditions")+facet_wrap("ID")
  ggsave(here::here("Output","Figures", 'Enzymes_abundace_proteomics_across_conditions.pdf') ) 
  
  Volcano_DFs <- map(.x = 1:3,
                     ~openxlsx::read.xlsx(here::here("Datasets","Processed", "Batch_corrected_DEPs.xlsx"), sheet = .x, rowNames = T)) %>% 
      set_names(openxlsx::getSheetNames(here::here("Datasets","Processed", "Batch_corrected_DEPs.xlsx")))
  Significant_proteins <- map(.x = Volcano_DFs, ~.x %>% 
                                  subset(p.adj<0.05) %>% 
                                  pull(Uniprot)) %>% unlist() %>% unique()
  
  Significant_genes<- map(.x = Volcano_DFs, ~.x %>% 
                              subset(p.adj<0.05) %>% 
                              pull(ID)) %>% unlist() %>% unique()
  
  
  Data_matrices$Imputted %>% 
      .[rownames(.) %chin% Significant_proteins,] %>% 
      pheatmap(scale = "row", cluster_cols = F)
  ##### CORREP#####
  Significant_proteins_matrix <- Data_matrices$Imputted %>% 
      as.data.frame() %>% 
      rownames_to_column("ProteinGroup") %>% 
      subset(ProteinGroup %in% Significant_proteins) %>% 
      mutate(Uniprot= ProteinGroup %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
      left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot))) %>% 
      #mutate(duplicated = BiocGenerics::duplicated(Uniprot))
      ungroup 
  
  Significant_proteins_names <- Significant_proteins_matrix %>% pull(ID, ProteinGroup)
  Significant_proteins_matrix <-  Significant_proteins_matrix %>%  column_to_rownames("ProteinGroup") %>%
      # set_names(.,str_remove_all(colnames(.),"_[:digit:]$")) %>%
      dplyr::select(where(is.numeric)) %>%
      mutate(Rowmean = rowMeans(.),
             across(where(is.numeric),~.x- Rowmean)) %>% 
      dplyr::select(-Rowmean) %>% 
      as.matrix() 
  paletteLength <- 20
  myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  # length(breaks) == length(paletteLength) + 1
  # use floor and ceiling to deal with even/odd length pallettelengths
  myBreaks <- c(seq(min(Significant_proteins_matrix, na.rm = T), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(Significant_proteins_matrix, na.rm = T)/paletteLength, max(Significant_proteins_matrix, na.rm = T),
                    length.out=floor(paletteLength/2)))
  
  
  Correp_correlation <- correp(Data_matrices$Imputted %>% dplyr::select(matches("_(1|2|3)$")),Significant_proteins)
  pdf(here::here("Output","Figures","Cluster_heatmap.pdf"))
  Correp_correlation$heatmap
  dev.off()
  rowMeans <- Data_matrices$Imputted %>% rowMeans()
  Interesting_Proteins <- c(c("RAD18","BARD1", "TOP2A","BRCA1","PCNA","PSMD2","NUCKS1","FXR1","DNAJC2", "TGM2"),HUMAN_9606$ID %>% 
                                str_subset("^NCAP|^SMC")) %>% unlist()
  clusters_df <- Data_matrices$Imputted %>% 
      as.data.frame() %>% 
      mutate(across(everything(), ~.x-rowMeans)) %>% 
      rownames_to_column("gene") %>% 
      pivot_longer(-gene, names_to = "condition", values_to = "abundance") %>% 
      inner_join(Correp_correlation$clusters %>% enframe("gene","cluster")) %>% 
      mutate(cluster = paste0("Cluster ", cluster),
             cluster_title = case_when(
                 cluster == "Cluster 1" ~ "Depleted Upon Etop Release",
                 cluster == "Cluster 2" ~ "Accummulates Upon Treatment and drops on Release",
                 cluster == "Cluster 3" ~ "Accumulates Upon Etop Treatment",
                 cluster == "Cluster 4" ~ "Depleted Upon Treatment and Recovers on Release",
                 cluster == "Cluster 5" ~ "Accumulates Upon Etop Release",
                 cluster == "Cluster 6" ~ "Depleted Upon Etop Treatment"
             ),
             condition = factor(condition, levels = .$condition %>% unique()%>%  sort)) %>% 
      left_join(Significant_proteins_names %>% enframe("gene","ID"))
  openxlsx::write.xlsx(clusters_df,overwrite = T, here::here("Datasets","Processed","correp_proteomics_clusters.xlsx"))
  clusters_df_avg <- clusters_df %>% 
      group_by(condition,cluster_title) %>% 
      summarise(abundance = median(abundance)) %>% 
      mutate(ID = "none")
  level_order <- clusters_df$condition %>% unique() %>% sort #this vector might be useful for other plots/analyses
  ggplot(data=clusters_df,aes(x = factor(condition,level = level_order),  y= abundance, colour = cluster_title, group = gene , label = ID))+
      geom_line(alpha = 0.3)+
      geom_point(alpha = 0.3)+
      facet_wrap("cluster_title")+
      theme(strip.text.x = element_text(size = 14, colour = "grey20"))+
      guides(colour = "none")+
      ggrepel::geom_label_repel(data = . %>% subset((condition == "x24h_5") & (ID %in% Interesting_Proteins)) %>% 
                                    group_by(ID) %>%  slice(which.min(abundance)))+
      geom_line(data=clusters_df_avg, aes(x=condition, y=abundance, group = cluster_title ), colour="grey10", alpha = 0.8 )+
      scale_x_discrete(labels= case_when( str_detect(level_order,"dmso")~ "DMSO",
                                          str_detect(level_order,"x0")~ "ETOP",
                                          str_detect(level_order,"x24")~ "Release"),guide = guide_axis(n.dodge = 2))+
      scale_color_manual(values = c("Depleted Upon Etop Treatment" = "darkblue",
                                    "Depleted Upon Etop Release" = "darkblue",
                                    "Accumulates Upon Etop Release"= "darkred",
                                    "Accumulates Upon Etop Treatment" ="darkred",
                                    "Depleted Upon Treatment and Recovers on Release" = "darkblue",
                                    "Accummulates Upon Treatment and drops on Release" = "darkred"))+
      ggtitle("Etoposide-responsive Proteins Cluster into different Behaviours")+theme_bw()
  ggsave(here::here("Output","Figures","Significant Protein Clusters.pdf"))
  ggplot(data=clusters_df %>% mutate(RAD_18 = if_else(ID == "RAD18",T,F)),
         aes(x = factor(condition,level = level_order),  y= abundance, colour = RAD_18, group = gene , label = ID))+
      geom_line(alpha = 0.3)+
      geom_point(alpha = 0.3)+
      facet_wrap("cluster")+
      theme(strip.text.x = element_text(size = 14, colour = "grey20"))+
      guides(colour = "none")+
      ggrepel::geom_label_repel(data = . %>% subset((condition == "x24h_5") & (ID %in% Interesting_Proteins)) %>% 
                                    group_by(ID) %>%  slice(which.min(abundance)))+
      geom_line(data=clusters_df_avg, aes(x=condition, y=abundance, group = cluster ), colour="grey10", alpha = 0.8 )+
      scale_x_discrete(labels= case_when( str_detect(level_order,"dmso")~ "DMSO",
                                          str_detect(level_order,"x0")~ "ETOP",
                                          str_detect(level_order,"x24")~ "Release"),guide = guide_axis(n.dodge = 2))+
      scale_color_manual(values = c(`TRUE`= "darkred",
                                    `FALSE` ="grey90" ))+
      ggtitle("Etoposide-responsive Proteins Cluster into different Behaviours")+theme_bw()
  
  Volcano_DFs$x0h_vs_x24h_diff %>% 
      mutate(ID = if_else(Uniprot == "P11388;Q02880","TOP2A; TOP2B",ID)) %>% 
      # dplyr::mutate(Metabolic_library = dplyr::if_else(Metabolic_library == "Non-Metabolic", "Non-Metabolic","Metabolic")) %>% 
      ggplot(aes(x = -log2_FC, y = -log10(p.val), label = ID, colour = Metabolic_library))+
      geom_point(data = . %>% subset(Metabolic_library =="Non-Metabolic" ), alpha = 0.05)+
      geom_point(data = . %>% subset(Metabolic_library !="Non-Metabolic" ), alpha = 0.1)+
      geom_point(data = . %>% subset(ID %in% c("BRCA1","TOP2A","SMC4","PSMD2", "ATF3")))+
      geom_point(data = . %>% subset(significant == T))+
      # geom_point(data = . %>% subset(Imputted_comparison == T), colour = "grey50")+
      #ggrepel::geom_label_repel(data = . %>% subset(ID %in% c("BRCA1","TOP2A","SMC4","PSMD2", "ATF3","TOP2A; TOP2B","PCNA")))+
      ggrepel::geom_label_repel(data = . %>% subset(significant == T&  ((Uniprot %in% Interesting_proteins$Uniprot)|
                                                                            (ID %in% c("BRCA1","TOP2A","SMC4","PSMD2", "ATF3","TOP2A; TOP2B","PCNA")))))+ 
      ggrepel::geom_label_repel(data = . %>% subset(significant == T&  Metabolic_library!="Non-Metabolic"))+
      scale_color_manual(values = c(`Non-Metabolic` = "darkblue",
                                    `Metabolic_library ;Human_GEM` = "#610101",
                                    `Human_GEM` = "#B80202",
                                    `Metabolic_library` = "#F7BD00"))+
      labs(x = "Fold Change on Chromatin", y = "-log10(pvalue)")+
      # annotate("text", x = c(-0.5,0.5), y=0, label = rev(conditions))+
      annotate("text", x = c(-1.5,1.5), y=0, label = rev(c("Release 24Hrs","Etop 3hrs")))+
      lims(x = c(-3.19,3.19))+
      theme(#legend.position = "none",
          panel.background = element_blank())+
      ggtitle(glue::glue("Chromatome Changes upon Etoposide Release"),
              subtitle = "Many Metabolism-Related Protein Detected on Chromatin")
  ggsave(here::here("Output","Proteomic_Volcano_x0h_vs_x24h_diff.pdf"))
  Human_hsa <- inner_join( HUMAN_9606 %>% 
                               subset(Type == "KEGG") %>% 
                               mutate(ID = str_remove_all(ID,"hsa:")) %>% 
                               dplyr::select(-Type) %>% 
                               dplyr::rename(KEGG = ID) %>% 
                               subset(!duplicated(Uniprot)),
                           HUMAN_9606 %>% 
                               subset(Type == "Gene_Name") %>% 
                               dplyr::select(-Type) %>% 
                               subset(!duplicated(Uniprot))) %>% 
      dplyr::select(-Uniprot) %>% 
      subset(!duplicated(KEGG)) %>% 
      subset(!duplicated(ID))
  convert_KEGG_into_Symbol <- function(x){
      # x <- "163/476/483/481/476/1213/1175/160/161/5566/1211"
      subset(Human_hsa, KEGG %in% ( x %>% str_split("/",simplify = T))) %>% pull(ID) %>% unique()
  }
  core_enrichment_KEGG <- openxlsx::read.xlsx(here::here("Datasets","Processed","KEGG-Multi-Omic_reduced.xlsx")) %>% 
      mutate(Genes = purrr::map(core_enrichment,convert_KEGG_into_Symbol))
  Volcano_DFs$dmso_vs_x24h_diff %>% 
      mutate(ID = if_else(Uniprot == "P11388;Q02880","TOP2A; TOP2B",ID),
             KEGG = case_when(
                 ID %in% unlist(subset(core_enrichment_KEGG,
                                       Description == "Cell cycle" & comparison == "dmso_vs_x24h_diff") %>%
                                    pull(Genes))~"Cell cycle",
                 TRUE~NA_character_
             ),
             Significant_alpha = if_else(significant ==T|!is.na(KEGG),T,F)) %>% 
      dplyr::mutate(Metabolic_library = dplyr::if_else(Metabolic_library == "Non-Metabolic", "Non-Metabolic","Metabolic")) %>% 
      ggplot(aes(x = -log2_FC, y = -log10(p.val), label = ID, colour = KEGG, alpha = Significant_alpha))+
      # geom_point(data = . %>% subset(Metabolic_library =="Non-Metabolic" ), alpha = 0.05)+
      # geom_point(data = . %>% subset(Metabolic_library !="Non-Metabolic" ), alpha = 0.1)+
      # geom_point(data = . %>% subset(ID %in% c("BRCA1","TOP2A","SMC4","PSMD2", "ATF3")))+
      geom_point(data = . %>% subset(is.na(KEGG)), size = 3.5)+
      geom_point(data = . %>% subset(is.na(KEGG) & Metabolic_library == "Metabolic"), size = 3.5,shape = 21,colour = "black",stroke  =1.5)+
      
      geom_point(data = . %>% subset(!is.na(KEGG)), size = 3.5)+
      geom_point(data = . %>% subset(!is.na(KEGG) & Metabolic_library == "Metabolic"), size = 3.5,shape = 21,colour = "black",stroke  =1.5)+
      
      # geom_point(data = . %>% subset(Imputted_comparison == T), colour = "grey50")+
      #ggrepel::geom_label_repel(data = . %>% subset(ID %in% c("BRCA1","TOP2A","SMC4","PSMD2", "ATF3","TOP2A; TOP2B","PCNA")))+
      ggrepel::geom_text_repel(data = . %>% subset(significant == T),max.overlaps = 2, size = 8)+
      ggrepel::geom_text_repel(data = . %>% subset(!is.na(KEGG)), size = 8)+ 
      # scale_colour_manual(values = 
      #                       c("DNA maintenance - DNA damage response" = "#6EA376",
      #                         "Chemical carcinogenesis - reactive oxygen species" = "#A95049" ,
      #                         "Citrate cycle (TCA cycle)" ="#FFA100"),
      # )+
      labs(x = "Fold Change on Chromatin", y = "-log10(pvalue)")+
      # annotate("text", x = c(-0.5,0.5), y=0, label = rev(conditions))+
      annotate("text", x = c(-1.5,1.5), y=0, label = rev(c("24h recovery","DMSO")),size = 7, colour = "#8F9191")+
      lims(x = c(-3.19,3.19))+
      theme(legend.position = "none", panel.border = element_rect(colour = "black", fill=NA, size=1),
            text = element_text(size = 25),
            panel.background = element_blank()
      )+
      ggtitle(glue::glue("Genes causing enrichemnt of Cell cycle KEGG in DMSO vs 24hrs"),
              subtitle = "They are all anaphase related")
  cell_cycle <- openxlsx::read.xlsx("C:/Users/skourtis/OneDrive - CRG - Centre de Regulacio Genomica/Bioinformatics Projects/Cell_cycle_oni/Datasets/Processed/periodicity_proteins_data.xlsx") %>% 
      dplyr::select(ID,cluster) %>% unique() %>% 
      mutate(cluster = as.character(cluster)) %>% subset(!(ID == 'TOP2A' & cluster ==1))
  Volcano_DFs$dmso_vs_x24h_diff %>% 
      left_join(cell_cycle) %>% 
      mutate(ID = if_else(Uniprot == "P11388;Q02880","TOP2A; TOP2B",ID),
             # KEGG = case_when(
             #   ID %in%cell_cycle~"Cell cycle",
             #   TRUE~NA_character_
             # ),
             Significant_alpha = if_else(significant ==T|!is.na(cluster),T,F)) %>% 
      dplyr::mutate(Metabolic_library = dplyr::if_else(Metabolic_library == "Non-Metabolic", "Non-Metabolic","Metabolic")) %>% 
      ggplot(aes(x = -log2_FC, y = -log10(p.val), label = ID, colour = cluster, alpha = Significant_alpha))+
      geom_point()+facet_wrap("cluster")+labs(x = "Fold Change on Chromatin", y = "-log10(pvalue)")+
      # annotate("text", x = c(-0.5,0.5), y=0, label = rev(conditions))+
      labs( x ="DMSOvs24Hrs") +
      ggrepel::geom_text_repel(data = . %>% subset(!is.na(cluster) ), size = 5)+ 
      ggtitle("Periodicity genes in DMSO vs 24hrs")
  
  Volcano_DFs$dmso_vs_x24h_diff %>% 
      left_join(cell_cycle) %>% 
      mutate(ID = if_else(Uniprot == "P11388;Q02880","TOP2A; TOP2B",ID),
             # KEGG = case_when(
             #   ID %in%cell_cycle~"Cell cycle",
             #   TRUE~NA_character_
             # ),
             Significant_alpha = if_else(significant ==T|!is.na(cluster),T,F)) %>% 
      dplyr::mutate(Metabolic_library = dplyr::if_else(Metabolic_library == "Non-Metabolic", "Non-Metabolic","Metabolic")) %>% 
      ggplot(aes(x = -log2_FC, y = -log10(p.val), label = ID, colour = cluster, alpha = Significant_alpha))+
      geom_point()+facet_wrap("cluster")+labs(x = "Fold Change on Chromatin", y = "-log10(pvalue)")+
      # annotate("text", x = c(-0.5,0.5), y=0, label = rev(conditions))+
      labs( x ="DMSOvs24Hrs") +
      # ggrepel::geom_text_repel(data = . %>% subset(!is.na(cluster) ), size = 5)+ 
      ggtitle("Periodicity genes in DMSO vs 24hrs")
  Volcano_DFs$x0h_vs_x24h_diff %>% 
      mutate(ID = if_else(Uniprot == "P11388;Q02880","TOP2A; TOP2B",ID),
             KEGG = case_when(
                 ID %in% unlist(subset(core_enrichment_KEGG, 
                                       Description == "Chemical carcinogenesis - reactive oxygen species" & comparison == "x0h_vs_x24h_diff") %>%
                                    pull(Genes))~"Chemical carcinogenesis - reactive oxygen species",
                 ID %in% unlist(subset(core_enrichment_KEGG,
                                       Description == "Citrate cycle (TCA cycle)" & comparison == "D14_BP") %>%
                                    pull(Genes))~"Citrate cycle (TCA cycle)",
                 ID %in% c(unlist(subset(core_enrichment_KEGG,
                                         Description %in% c("Homologous recombination",
                                                            "Fanconi anemia pathway") & comparison == "x0h_vs_x24h_diff") %>%
                                      pull(Genes)),
                           unlist(subset(core_enrichment_BP,
                                         Description %in% c("cellular response to DNA damage stimulus","DNA replication") & comparison == "x0h_vs_x24h_diff") %>%
                                      pull(Genes)))~"DNA maintenance - DNA damage response",
                 # Gene %in% unlist(subset(core_enrichment_KEGG,
                 #                         Description == "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis" & comparison == "dmso_vs_x24h_diff") %>%
                 #                      pull(Genes))~"Glycosylphosphatidylinositol (GPI)-anchor biosynthesis",
                 TRUE~NA_character_
             ),
             Significant_alpha = if_else(significant ==T|!is.na(KEGG),T,F)) %>% 
      dplyr::mutate(Metabolic_library = dplyr::if_else(Metabolic_library == "Non-Metabolic", "Non-Metabolic","Metabolic")) %>% 
      ggplot(aes(x = -log2_FC, y = -log10(p.val), label = ID, colour = KEGG, alpha = Significant_alpha))+
      # geom_point(data = . %>% subset(Metabolic_library =="Non-Metabolic" ), alpha = 0.05)+
      # geom_point(data = . %>% subset(Metabolic_library !="Non-Metabolic" ), alpha = 0.1)+
      # geom_point(data = . %>% subset(ID %in% c("BRCA1","TOP2A","SMC4","PSMD2", "ATF3")))+
      geom_point(data = . %>% subset(is.na(KEGG)), size = 3.5)+
      geom_point(data = . %>% subset(is.na(KEGG) & Metabolic_library == "Metabolic"), size = 3.5,shape = 21,colour = "black",stroke  =1.5)+
      
      geom_point(data = . %>% subset(!is.na(KEGG)), size = 3.5)+
      geom_point(data = . %>% subset(!is.na(KEGG) & Metabolic_library == "Metabolic"), size = 3.5,shape = 21,colour = "black",stroke  =1.5)+
      
      # geom_point(data = . %>% subset(Imputted_comparison == T), colour = "grey50")+
      #ggrepel::geom_label_repel(data = . %>% subset(ID %in% c("BRCA1","TOP2A","SMC4","PSMD2", "ATF3","TOP2A; TOP2B","PCNA")))+
      ggrepel::geom_text_repel(data = . %>% subset(significant == T&  is.na(KEGG)),max.overlaps = 2, size = 8)+
      ggrepel::geom_text_repel(data = . %>% subset(significant == T&  !is.na(KEGG)), size = 8)+ 
      scale_colour_manual(values = 
                              c("DNA maintenance - DNA damage response" = "#6EA376",
                                "Chemical carcinogenesis - reactive oxygen species" = "#A95049" ,
                                "Citrate cycle (TCA cycle)" ="#FFA100"),
      )+
      labs(x = "Fold Change on Chromatin", y = "-log10(pvalue)")+
      # annotate("text", x = c(-0.5,0.5), y=0, label = rev(conditions))+
      annotate("text", x = c(-1.5,1.5), y=0, label = rev(c("24h recovery","No release")),size = 7, colour = "#8F9191")+
      lims(x = c(-3.19,3.19))+
      theme(legend.position = "none", panel.border = element_rect(colour = "black", fill=NA, size=1),
            text = element_text(size = 25),
            panel.background = element_blank()
      )+
      ggtitle(glue::glue("Chromatome Changes upon Etoposide Release"),
              subtitle = "Many Metabolism-Related Protein Detected on Chromatin")
  ggsave(here::here("Output","Figures","Annotated_KEGG_Proteomic_Volcano_x0h_vs_x24h_diff_illustrator.pdf"), height = 10, width = 10)
  Volcano_DFs$dmso_vs_x0h_diff %>% 
      mutate(ID = if_else(Uniprot == "P11388;Q02880","TOP2A; TOP2B",ID),
             KEGG = case_when(
                 ID %in% unlist(subset(core_enrichment_KEGG, 
                                       Description == "Chemical carcinogenesis - reactive oxygen species" & comparison == "x0h_vs_x24h_diff") %>%
                                    pull(Genes))~"Chemical carcinogenesis - reactive oxygen species",
                 ID %in% unlist(subset(core_enrichment_KEGG,
                                       Description == "Citrate cycle (TCA cycle)" & comparison == "D14_BP") %>%
                                    pull(Genes))~"Citrate cycle (TCA cycle)",
                 ID %in% c(unlist(subset(core_enrichment_KEGG,
                                         Description %in% c("Homologous recombination",
                                                            "Fanconi anemia pathway") & comparison == "x0h_vs_x24h_diff") %>%
                                      pull(Genes)),
                           unlist(subset(core_enrichment_BP,
                                         Description %in% c("cellular response to DNA damage stimulus","DNA replication") & comparison == "x0h_vs_x24h_diff") %>%
                                      pull(Genes)))~"DNA maintenance - DNA damage response",
                 # Gene %in% unlist(subset(core_enrichment_KEGG,
                 #                         Description == "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis" & comparison == "dmso_vs_x24h_diff") %>%
                 #                      pull(Genes))~"Glycosylphosphatidylinositol (GPI)-anchor biosynthesis",
                 TRUE~NA_character_
             ),
             Significant_alpha = if_else(significant ==T|!is.na(KEGG),T,F)) %>% 
      dplyr::mutate(Metabolic_library = dplyr::if_else(Metabolic_library == "Non-Metabolic", "Non-Metabolic","Metabolic")) %>% 
      ggplot(aes(x = -log2_FC, y = -log10(p.val), label = ID, colour = KEGG, alpha = Significant_alpha))+
      # geom_point(data = . %>% subset(Metabolic_library =="Non-Metabolic" ), alpha = 0.05)+
      # geom_point(data = . %>% subset(Metabolic_library !="Non-Metabolic" ), alpha = 0.1)+
      # geom_point(data = . %>% subset(ID %in% c("BRCA1","TOP2A","SMC4","PSMD2", "ATF3")))+
      geom_point(data = . %>% subset(is.na(KEGG)), size = 3.5)+
      geom_point(data = . %>% subset(is.na(KEGG) & Metabolic_library == "Metabolic"), size = 3.5,shape = 21,colour = "black",stroke  =1.5)+
      
      geom_point(data = . %>% subset(!is.na(KEGG)), size = 3.5)+
      geom_point(data = . %>% subset(!is.na(KEGG) & Metabolic_library == "Metabolic"), size = 3.5,shape = 21,colour = "black",stroke  =1.5)+
      
      # geom_point(data = . %>% subset(Imputted_comparison == T), colour = "grey50")+
      #ggrepel::geom_label_repel(data = . %>% subset(ID %in% c("BRCA1","TOP2A","SMC4","PSMD2", "ATF3","TOP2A; TOP2B","PCNA")))+
      ggrepel::geom_text_repel(data = . %>% subset(significant == T&  is.na(KEGG)), size = 8)+ 
      ggrepel::geom_text_repel(data = . %>% subset(significant == T&  !is.na(KEGG)), size = 8)+ 
      # scale_color_manual(values = c(`Non-Metabolic` = "darkblue",
      #                               `Metabolic_library ;Human_GEM` = "#610101",
      #                               `Human_GEM` = "#B80202",
      #                               `Metabolic_library` = "#F7BD00"))+
      labs(x = "Fold Change on Chromatin", y = "-log10(pvalue)")+
      # annotate("text", x = c(-0.5,0.5), y=0, label = rev(conditions))+
      annotate("text", x = c(-1.5,1.5), y=0, label =rev(c("No release","DMSO")),size = 7, colour = "#8F9191")+
      lims(x = c(-3.8,3.8))+
      theme(legend.position = "none", panel.border = element_rect(colour = "black", fill=NA, size=1),
            text = element_text(size = 25),
            panel.background = element_blank()
      )+
      scale_colour_manual(values = 
                              c("DNA maintenance - DNA damage response" = "#6EA376",
                                "Chemical carcinogenesis - reactive oxygen species" = "#A95049" ,
                                "Citrate cycle (TCA cycle)" ="#FFA100"),
      )+
      ggtitle(glue::glue("Chromatome Changes upon Etoposide Release"),
              subtitle = "Many Metabolism-Related Protein Detected on Chromatin")
  ggsave(here::here("Output","Figures","Annotated_KEGG_Proteomic_Volcano_x0h_vs_DMSO_diff_Illustrator.pdf"), height = 10, width = 10)
  
  
  # Volcano_DFs$x0h_vs_x24h_diff %>% 
  #   mutate(ID = if_else(Uniprot == "P11388;Q02880","TOP2A; TOP2B",ID)) %>% 
  #   # dplyr::mutate(Metabolic_library = dplyr::if_else(Metabolic_library == "Non-Metabolic", "Non-Metabolic","Metabolic")) %>% 
  #   ggplot(aes(x = -log2_FC, y = -log10(p.val), label = ID, size = -log10(p.val)))+
  # 
  #   geom_point(data = . %>% subset(significant == T), colour = "black")+
  #   geom_point(data = . %>% subset(significant == F), colour = "black", alpha = 0.1)+
  #   # ggrepel::geom_label_repel(data = . %>% subset(ID %in% c("IMPDH2")), size = 4)+
  #   # scale_color_manual(values = c(`Non-Metabolic` = "darkblue",
  #   #                               `Metabolic_library ;Human_GEM` = "#610101",
  #   #                               `Human_GEM` = "#B80202",
  #   #                               `Metabolic_library` = "#F7BD00"))+
  #   labs(x = "Fold Change on Chromatin", y = "-log10(pvalue)")+
  #   # annotate("text", x = c(-0.5,0.5), y=0, label = rev(conditions))+
  #   # annotate("text", x = c(-1.5,1.5), y=0, label = rev(c("Release 24Hrs","Etop 3hrs")))+
  #   lims(x = c(-3.19,3.19))+
  #   theme(#legend.position = "none",
  #     panel.background = element_blank())+
  #   ggtitle(glue::glue("Chromatome Changes upon Etoposide Release"),
  #           subtitle = "Many Metabolism-Related Protein Detected on Chromatin")
  # ggsave(here::here("Output","Proteomic_Volcano_x0h_vs_x24h_diff_IMPDH2.pdf"))
  
  Volcano_DFs$dmso_vs_x24h_diff %>% 
      mutate(ID = if_else(Uniprot == "P11388;Q02880","TOP2A; TOP2B",ID)) %>% 
      dplyr::mutate(Protein_Type = dplyr::if_else(Metabolic_library == "Non-Metabolic", "Non-Metabolic","Metabolic")) %>% 
      ggplot(aes(x = -log2_FC, y = -log10(p.val), label = ID, colour = Metabolic_library))+
      geom_point(data = . %>% subset(Protein_Type =="Non-Metabolic" ), alpha = 0.05)+
      geom_point(data = . %>% subset(Protein_Type !="Non-Metabolic" ), alpha = 0.1)+
      # geom_point(data = . %>% subset(ID %in% c("BRCA1","TOP2A","SMC4","PSMD2", "ATF3")))+
      geom_point(data = . %>% subset(significant == T))+
      # geom_point(data = . %>% subset(Imputted_comparison == T), colour = "grey50")+
      #ggrepel::geom_label_repel(data = . %>% subset(ID %in% c("BRCA1","TOP2A","SMC4","PSMD2", "ATF3","TOP2A; TOP2B","PCNA")))+
      ggrepel::geom_label_repel(data = . %>% subset(significant == T&  ((Uniprot %in% Interesting_proteins$Uniprot)|
                                                                            (ID %in% c("BRCA1","TOP2A","SMC4","PSMD2", "ATF3","TOP2A; TOP2B","PCNA")))))+ 
      ggrepel::geom_label_repel(data = . %>% subset(significant == T&  Protein_Type!="Non-Metabolic"))+
      # scale_color_manual(values = c(`Non-Metabolic` = "darkblue",
      #                               `Metabolic` = "darkred"))+
      scale_color_manual(values = c(`Non-Metabolic` = "darkblue",
                                    `Metabolic_library ;Human_GEM` = "#610101",
                                    `Human_GEM` = "#B80202",
                                    `Metabolic_library` = "#F7BD00"))+
      labs(x = "Fold Change on Chromatin", y = "-log10(pvalue)")+
      # annotate("text", x = c(-0.5,0.5), y=0, label = rev(conditions))+
      annotate("text", x = c(-1.5,1.5), y=0, label = rev(c("Release 24Hrs","DMSO")))+
      lims(x = c(-3.19,3.19))+
      theme(#legend.position = "none",
          panel.background = element_blank())+
      ggtitle(glue::glue("Chromatome Changes upon Etoposide Release"),
              subtitle = "Many Metabolism-Related Protein Detected on Chromatin")
  ggsave(here::here("Output","Figures","Proteomic_Volcano_dmso_vs_x24h_diff.pdf"))
  
  Volcano_DFs$dmso_vs_x0h_diff %>% 
      mutate(ID = if_else(Uniprot == "P11388;Q02880","TOP2A; TOP2B",ID)) %>% 
      dplyr::mutate(Protein_Type = dplyr::if_else(Metabolic_library == "Non-Metabolic", "Non-Metabolic","Metabolic")) %>% 
      ggplot(aes(x = -log2_FC, y = -log10(p.val), label = ID, colour = Metabolic_library))+
      geom_point(data = . %>% subset(Protein_Type =="Non-Metabolic" ), alpha = 0.05)+
      geom_point(data = . %>% subset(Protein_Type !="Non-Metabolic" ), alpha = 0.1)+
      geom_point(data = . %>% subset(ID %in% interesing_proteins))+
      geom_point(data = . %>% subset(significant == T))+
      # geom_point(data = . %>% subset(Imputted_comparison == T), colour = "grey50")+
      ggrepel::geom_label_repel(data = . %>% subset(significant == T&  ((Uniprot %in% Interesting_proteins$Uniprot)|
                                                                            (ID %in% c("BRCA1","TOP2A","SMC4","PSMD2", "ATF3","TOP2A; TOP2B","PCNA",
                                                                                       "PSMC3","DNAJC2","NCAPH", "NCAPG","MDK")))))+ 
      ggrepel::geom_label_repel(data = . %>% subset(significant == T &Metabolic_library != "Non-Metabolic"))+
      # scale_color_manual(values = c(`Non-Metabolic` = "darkblue",
      #                               `Metabolic` = "darkred"))+
      scale_color_manual(values = c(`Non-Metabolic` = "darkblue",
                                    `Metabolic_library ;Human_GEM` = "#610101",
                                    `Human_GEM` = "#B80202",
                                    `Metabolic_library` = "#F7BD00"))+
      labs(x = "Fold Change on Chromatin", y = "-log10(pvalue)")+
      # annotate("text", x = c(-0.5,0.5), y=0, label = rev(conditions))+
      annotate("text", x = c(-1.5,1.5), y=0, label = rev(c("Etop 3 Hrs","DMSO")))+
      ggtitle("Chromatome Changes upon Etoposide Treatment") +
      lims(x = c(-3.19,3.19))+
      theme(#legend.position = "none",
          panel.background = element_blank())
  
  ggsave(here::here("Output","Figures", "Proteomic_Volcano_dmso_vs_x0h_diff.pdf"))
  
  
  #Multiple GSEA
  list_of_diff_genes <- clusters_df %>% 
      group_split(cluster) %>% 
      set_names(.,map_chr(.x = .,~paste0("cluster_",.x %>% pull(cluster) %>% unique()))) %>% 
      map(.x = .,~.x %>% pull(gene) %>% str_remove_all(";[:graph:]*$")%>% unique() )
  enrich_clusters <- function(gene_list){
      enrichGO(gene         = gene_list,
               OrgDb         = org.Hs.eg.db,
               universe = rownames(Data_matrices$Imputted)%>% str_remove_all(";[:graph:]*$")%>% unique(),
               keyType       = 'UNIPROT',
               ont           = "BP",
               # minGSSize = 5,
               # maxGSSize = 50,
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.01,
               qvalueCutoff  = 0.05)
  }
  ck <- compareCluster(geneCluster = list_of_diff_genes, fun = enrich_clusters )
  terms_to_reduce <- ck@compareClusterResult %>% subset(`p.adjust`<0.05 ) %>% 
      dplyr::select(geneID,Description) %>% subset(!duplicated(Description)) %>% 
      pull(geneID,Description)  %>% purrr::map(str_split,"/") %>% flatten()
  reduced_terms_CC <- Jaccard_index_list(terms_to_reduce,max_jacc = 0.8)
  ck@compareClusterResult <- ck@compareClusterResult %>% subset(Description %in% (reduced_terms_CC %>% names()) )
  
  ck@compareClusterResult %>% #View()
      mutate(GeneRatio = sapply(strsplit(GeneRatio, split = "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))) %>%
      group_by(Cluster) %>% 
      arrange(p.adjust) %>% 
      slice(1:20) %>%
      arrange(Cluster) %>% mutate(
          Description = factor(Description)
      ) %>% 
      ggplot(aes(x = Cluster, y= reorder(Description,Cluster), colour = p.adjust, size = GeneRatio))+
      geom_point()+
      ggtitle("CORREP Clusters Enrichment BP simplified")
  ggsave(here::here("Output", "Etop_DIA_EC_two_methods","latest","Significant Protein Clusters ORA.png"))
  
  
  
  Clusters <- clusters_df %>% 
      dplyr::select(gene,cluster) %>% 
      distinct() %>% 
      arrange(cluster) %>% 
      column_to_rownames("gene") %>% 
      dplyr::select(cluster)
  
  png(here::here(output_folder,glue::glue("Heatmap_Significant ",dataset_name,".png")), width = 2500, height = 3800,res  =300) 
  pheatmap::pheatmap(Significant_proteins_matrix[Clusters %>% rownames(),Significant_proteins_matrix %>% colnames() %>% sort()],cluster_cols = F,fontsize_row = 6, 
                     clustering_distance_rows = "euclidean", cluster_rows = F,labels_col= c(rep(c("DMSO","Etop","Release"),times = c(4,6,6))),
                     annotation_row =  Clusters, labels_row = Significant_proteins_names,legend_labels = "Abundance",legend = TRUE,
                     gaps_row = cumsum(Clusters$cluster %>% table()),
                     # scale = "row",
                     main = glue::glue("Etoposide-responsive Proteins on Chromatin"),color=myColor, breaks=myBreaks)
  dev.off()
  # BiocManager::install("OmnipathR")
  library(OmnipathR)
  ## We query and store the interactions into a dataframe
  resources_to_use <- get_interaction_resources() %>% str_subset("BioG|HPRD|IntAct|CORUM")#|MPPI|IME
  interactions <-
      import_all_interactions(resources=resources_to_use,
                              organism = 9606)
  ## Warning in omnipath_check_param(param): The following resources are not available: STRING. Check the resource
  ## names for spelling mistakes.
  
  selected_genes <- All_interesting_genes
  selected_genes <- c(sort_unsort_hits,Depleted_d14_genes) %>% unique()
  interactions_Screen_hits <- interactions %>% subset(source_genesymbol %in% selected_genes|
                                                          target_genesymbol %in% selected_genes) %>% 
      dplyr::select(source_genesymbol, target_genesymbol) %>% distinct()
  interactions_Screen_hits <- rbind(interactions_Screen_hits,
                                    data.frame(source_genesymbol = c("PRDX1","PRDX1"),
                                               target_genesymbol = c("TERF1","PRNP")))
  #Cyto_network 
  # intersect(Volcano_DFs$dmso_vs_x0h_diff$ID, Candidates_Joanna)
  library(RCy3)
  library(igraph)
  cytoscapePing()
  Hits <- read.csv(here::here("Datasets","Processed","Hits.csv")) %>% deframe()
  Hits_names <- names(Hits)
  mypal <- colorRampPalette( c( scales::muted("blue"),"white", scales::muted( "red" )) )(50 )
  map2color<-function(x,pal,limits=NULL){
      if(is.null(limits)) limits=range(x)
      pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
  }
  map2color_center<-function(x,vector_low_mid_high){
      # if(is.null(limits)) limits=range(x)
      # vector_low_mid_high <- c(low = "blue",
      #                          mid = "white",
      #                          high = "red")
      # pal = mypal
      # x = Volcano_DFs$x0h_vs_x24h_diff$log2_FC
      # vector_low_mid_high = c("blue","white","red")
      # x = Volcano_DFs$dmso_vs_x0h_diff$log2_FC 
      names(x) <- 1:length(x)
      pal_high <- colorRampPalette( vector_low_mid_high[-1] )(100 )
      pal_low <- colorRampPalette( vector_low_mid_high[-3] )(100 )
      x_high <- x[x>0]
      x_low <- x[x<=0]
      limits_high=range(x_high)
      limits_low=range(x_low)
      x_high <-  pal_high[findInterval(x_high,seq(limits_high[1],limits_high[2],length.out=length(pal_high)+1), all.inside=TRUE)] %>% 
          set_names(names(x_high))
      x_low <-  pal_low[findInterval(x_low,seq(limits_low[1],limits_low[2],length.out=length(pal_low)+1), all.inside=TRUE)]%>% 
          set_names(names(x_low))
      c(x_high,x_low) %>% .[names(x)] %>% as.vector() 
      
  }
  #Here I had to remove all columns from edges and almost all from nodes to allow import
  Chromatin_genes <- Volcano_DFs %>% purrr::map(.x =.,~pull(.x,ID)) %>% unlist
  visnet <- list(nodes = read.csv(here::here("Hits_screen default  node.csv")) %>% 
                     dplyr::select(name) %>% distinct() %>% dplyr::rename(id = name),
                 edges = read.csv(here::here("Hits_screen default  edge.csv"))  %>% 
                     dplyr::select(name) %>% separate(name, into = c("from","extra"), sep = " \\(") %>%
                     separate(extra, into = c("extra","to"), sep = "\\) ") %>% dplyr::select(from,to))
  
  
  # suid <- getNetworkSuid()
  #remove extra columns from cytoscape and duplicate nodes
  # visnet <- createIgraphFromNetwork(suid) %>% toVisNetworkData
  # visnet <- list(nodes = data.frame(
  #   id = c(interactions_Screen_hits$source_genesymbol,interactions_Screen_hits$target_genesymbol) %>% unique()),
  #                edges = interactions_Screen_hits %>% set_names(c("from","to")))
  
  visnet$edges <- visnet$edges %>% 
      #subset(from %in% genes & to %in% genes) %>% 
      mutate(
          Hit = case_when(from %in% Hits_names ~ from,
                          to %in% Hits_names ~ to,
                          TRUE~"no_hit_involved"),
          Partner = case_when(from %in% Hits_names ~ to,
                              to %in% Hits_names ~ from,
                              TRUE~"no_hit_involved"),
      ) %>% subset(Hit != "no_hit_involved") %>% dplyr::select(from, to, Hit,Partner) %>% distinct() 
  Connecting_nodes <- c(visnet$edges$from,visnet$edges$to) %>% table() %>% .[.>1] %>% names 
  Connecting_nodes <- Connecting_nodes[!(Connecting_nodes %in% Hits_names )]
  visnet$edges <-   mutate(visnet$edges ,
                           Connecting = if_else(Partner %in% Connecting_nodes | Hit%in% Connecting_nodes ,1,0),
                           on_chrom = if_else(Partner %in% Chromatin_genes,1,0),
                           Hits_direct = if_else(Hit %in% Hits_names  & Partner %in% Hits_names,3,0)) %>% 
      mutate(         importance = as.numeric(Connecting) +as.numeric(on_chrom)+as.numeric(Hits_direct)) %>% 
      group_by(Hit) %>% arrange(Hit,-importance) %>% distinct(Partner,Hit,.keep_all = T) %>% 
      slice_head(n=10)
  # visnet$edges <-  visnet$edges %>% subset(Hit %in% names(Hits[Hits>0]))
  visnet$nodes <-  visnet$nodes %>% 
      left_join(Volcano_DFs$dmso_vs_x24h_diff %>% dplyr::select(ID,log2_FC) %>% 
                    mutate(color = map2color_center(-log2_FC,vector_low_mid_high = c(scales::muted("blue"),"#f2f5f7","#9e3a3a")))
                
                %>% distinct(ID, .keep_all= T) , by = c("id" = "ID")) %>% 
      #dplyr::select(id,label, matches("Int")) %>% 
      subset(id %in% c(visnet$edges$from,visnet$edges$to)) 
  visnet$nodes <-  mutate(visnet$nodes,shape = if_else(id %in% Hits_names,"rect","dot"),
                          color.border = if_else(!is.na(color),"grey70", "gray90"),
                          color.border = if_else(id %in% Hits_names &  !is.na(color),"gray40",color.border),
                          color.background = if_else(!is.na(color),color, "white"),
                          
                          value = if_else(id %in% Hits_names,15,7),
                          `font.size` = (9*value),
                          label = id) %>% 
      mutate(label = if_else(
          # (label %in% Chromatin_genes) & (label %in% Connecting_nodes)|
          # (abs(log2_FC)>0.3) & (label %in% Chromatin_genes) |
          label %in%Hits_names,label,NA_character_ ))
  # left_join(median_enrichement, by = c("id"= "ID")) %>% 
  # mutate(color = map2color(Abundance %>% replace_na(0),mypal),
  #        color = if_else(is.na(Abundance), "grey80", color )
  library(visNetwork)
  visNetwork(visnet$nodes, visnet$edges, height = "1500px", width = "1500px", main = "Positively selected Genes PPI on Chrom, 24vsDMSO, red = more on chrom at 24hrs, ")  
  # install.packages('igraph')
  # library(igraph)
  library(ggraph)
  visnet$nodes %>% mutate()
  net <- graph_from_data_frame(d=visnet$edges, vertices=visnet$nodes, directed=F) 
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  
  ggraph(net, layout="stress") +
      geom_edge_fan(color="gray90", width=0.3) + 
      geom_node_point(shape = 21, fill=V(net)$color.border, size=V(net)$value,
                      # labellabel = V(net)$label,
                      colour = V(net)$color.background ,stroke = 3) +
      geom_node_text(label = V(net)$label, size=10, color="gray30", repel=T) +
      theme_void()
  ggsave(here::here("Output","Figures","PPI_network_high.pdf"), width = 15,height = 10)
  visnet$nodes %>% ggplot(aes(x = id, colour = -log2_FC, y = log2_FC))+
      geom_point()+
      scale_colour_gradient2(mid =c(scales::muted("blue"),"#f2f5f7","#9e3a3a")[2] , 
                             high = c(scales::muted("blue"),"#f2f5f7","#9e3a3a")[3],low = c(scales::muted("blue"),"#f2f5f7","#9e3a3a")[1] )
  ggsave(here::here("Output","Figures","PPI_network_high_scale.pdf"), width = 15,height = 10)
  
  KEGG_enrichment <- function(df_input){
      # KEGG_ <- s_prot
      gene_list <- df_input %>%
          arrange(-log2_FC) %>%   pull(log2_FC,KEGG) %>%na.omit()
      
      kk2 <- gseKEGG(geneList     = gene_list,
                     organism     = 'hsa',
                     minGSSize    = 10,
                     pvalueCutoff = 1,
                     verbose      = FALSE)
      
  }
  function_enchrichment_all <- function(df_input){
      gseGO(geneList     =df_input %>% pull(log2_FC,ID) %>% sort(decreasing = T) ,
            OrgDb        = org.Hs.eg.db,
            ont          = "BP",
            keyType = "SYMBOL",
            #nPerm        = 1000,
            minGSSize    = 10,
            maxGSSize    = 500,
            pvalueCutoff = 1,
            verbose      = FALSE) }
  Volcano_DFs[["D14_BP"]] <- norm_genes %>% dplyr::rename(ID = Gene,
                                                          log2_FC =Etop_vs_DMS0 )
  Volcano_DFs[["D10_high_BP"]] <- Screens_new$Etop_High_low_10_cnv_negsg.mle.gene_summary_fixed_names.txt %>% 
      dplyr::rename(ID = Gene,
                    log2_FC =`high|beta` )
  
  enrichment <- map(Volcano_DFs[4],function_enchrichment_all)
  # enrichment[["D10_high_BP"]] <- D10_high_BP
  # enrichment[["D14_BP"]] <- D14_BP
  
  enrichment_df_BP <- enrichment %>% imap_dfr(.x = ., ~.x@result%>% as.data.frame %>% 
                                                  dplyr::select(Description,NES,p.adjust,core_enrichment)  %>% mutate(comparison = .y)) %>% 
      mutate(Significant = p.adjust<0.05) %>% 
      # subset(!(comparison %in%  c("dmso_vs_x24h_diff") )) %>% 
      group_by(Description) %>% 
      add_tally(Significant) %>% subset(n>0)%>% 
      # subset(Description %in%Simplified_terms) %>% 
      mutate(Type = case_when(
          p.adjust<0.05 & NES>0 ~ "Signi_Up",
          p.adjust<0.05 & NES<0 ~ "Signi_Down",
          p.adjust>0.05~ "Not_Signi",
          TRUE ~"Error"
      ))
  write.xlsx(enrichment_df_BP,here::here("Datasets","Processed","BP-Multi-Omic.xlsx"), overwrite =T)
  enrichment_df_BP <- openxlsx::read.xlsx(here::here("Datasets","Processed","BP-Multi-Omic.xlsx")) %>% mutate(NES = if_else(comparison == "D14_BP",-NES,NES),
                                                                                                              comparison =if_else(comparison == "D14_BP", "D14_BP_rev",comparison))
  terms_to_reduce <- enrichment_df_BP %>% subset(`p.adjust`<0.05 ) %>% 
      dplyr::select(core_enrichment,Description) %>% subset(!duplicated(Description)) %>% 
      pull(core_enrichment,Description)  %>% purrr::map(str_split,"/") %>% flatten()
  reduced_terms <- Jaccard_index_list(terms_to_reduce,max_jacc = 0.3)
  enrichment_df_BP <- enrichment_df_BP %>% subset(Description %in% (reduced_terms %>% names()) )
  write.xlsx(enrichment_df_BP,here::here("Datasets","Processed","BP-Multi-Omic_reduced.xlsx"), overwrite =T)
  Interesting_terms <- rbind(subset(enrichment_df_BP, n== 1 & `p.adjust`<0.05 ) %>% arrange(`p.adjust`) %>% 
                                 group_by(comparison) %>%
                                 slice(1:5),
                             subset(enrichment_df_BP, n> 1 & `p.adjust`<0.05 ) %>% arrange(`p.adjust`) %>% 
                                 group_by(comparison) %>%
                                 slice(1:3)) %>% pull(Description)
  enrichment_df_BP <- subset(enrichment_df_BP,Description %in% Interesting_terms)%>% 
      # subset(comparison != "dmso_vs_x24h_diff") %>% 
      mutate(NES = if_else(`p.adjust`<0.05,NES,0)) %>% 
      # mutate(NES = if_else(str_detect(comparison,"vs"),NES,-NES)) %>% 
      arrange(n,comparison) %>% 
      mutate(Description  = factor(Description, Description %>% unique))
  ggplot(enrichment_df_BP,aes(x = comparison,y  = reorder(Description,n), size = -log10(p.adjust), colour = NES))+
      geom_point()+
      theme_bw()+
      scale_colour_gradient2(high= muted("red"), low = muted("blue"))+
      scale_x_discrete(
          # labels=c("D10_high_BP"  = "High_vs_unsort",
          #                         "D14_BP"= 'D14Etop_vs_D14DMSO',
          #                         "x0h_vs_x24h_diff" = "Release_vs_Etop",
          #                         "dmso_vs_x0h_diff" = "Etop_vs_DMSO",
          #                         "dmso_vs_x24h_diff" = "Release_vs_DMSO"),
          guide = guide_axis(n.dodge = 2))+
      ggtitle("BP Biological Response to Etoposide Release")
  ggsave(here::here("Output","Figures","BP-Multi-Omic_reduced.pdf"), width  = 10, height = 15)
  
  
  Human_hsa_screen <- inner_join( HUMAN_9606 %>% 
                                      subset(Type == "KEGG") %>% 
                                      mutate(ID = str_remove_all(ID,"hsa:")) %>% 
                                      dplyr::select(-Type) %>% 
                                      dplyr::rename(KEGG = ID) %>% 
                                      subset(!duplicated(Uniprot)),
                                  HUMAN_9606 %>% 
                                      subset(Type == "Gene_Name") %>% 
                                      dplyr::select(-Type) %>% 
                                      subset(!duplicated(Uniprot))) %>% 
      dplyr::select(-Uniprot) %>% 
      subset(!duplicated(KEGG)) %>% 
      subset(!duplicated(ID))
  KEGG_volcanos_screen <- Volcano_DFs[4:5] %>% purrr::map(.x = .,~left_join(.x,Human_hsa_screen, by  = "ID") )
  KEGG_volcanos_prot <- Volcano_DFs[1:3]%>% purrr::map(.x = .,~left_join(.x %>% dplyr::select(-ID),Human_hsa, by  = c("Single_Uniprot"= "Uniprot") ) %>% dplyr::rename(KEGG= ID))
  KEGG_volcanos <- c(KEGG_volcanos_prot,KEGG_volcanos_screen)
  enrichment_KEGG <- map(KEGG_volcanos,KEGG_enrichment)
  enrichment_df <- enrichment_KEGG %>% imap_dfr(.x = ., ~.x@result%>% as.data.frame %>% 
                                                    dplyr::select(Description,NES,p.adjust,core_enrichment)  %>% mutate(comparison = .y)) %>% 
      mutate(Significant = p.adjust<0.05) %>% 
      group_by(Description) %>% 
      add_tally(Significant) %>% subset(n>0)%>% 
      # subset(Description %in%Simplified_terms) %>% 
      mutate(Type = case_when(
          p.adjust<0.05 & NES>0 ~ "Signi_Up",
          p.adjust<0.05 & NES<0 ~ "Signi_Down",
          p.adjust>0.05~ "Not_Signi",
          TRUE ~"Error"
      )) 
  enrichment_df <- enrichment_df %>% mutate(NES = if_else(comparison == "D14_BP",NES*(-1),NES))
  # ,
  # comparison =   case_when(
  #   comparison  == "D10_high_BP"  ~ "D10_high_BP",
  #   comparison  == "D14_BP"~ 'D14Etop_vs_D14DMSO_Rev',
  #   comparison  == "x0h_vs_x24h_diff" ~ "Release_vs_Etop",
  #   comparison  == "dmso_vs_x0h_diff" ~ "Etop_vs_DMSO",
  #   comparison  == "dmso_vs_x24h_diff" ~ "Release_vs_DMSO")
  # )
  
  write.xlsx(enrichment_df,here::here("Datasets","Processed","KEGG-Multi-Omic.xlsx"), overwrite = T)
  enrichment_df <- read.xlsx(here::here("Datasets","Processed","KEGG-Multi-Omic.xlsx")) 
  terms_to_reduce <- enrichment_df %>% subset(`p.adjust`<0.05 ) %>% 
      dplyr::select(core_enrichment,Description) %>% subset(!duplicated(Description)) %>% 
      pull(core_enrichment,Description)  %>% purrr::map(str_split,"/") %>% flatten()
  reduced_terms <- Jaccard_index_list(terms_to_reduce,max_jacc = 0.3)
  enrichment_df <- enrichment_df %>% subset(Description %in% (reduced_terms %>% names()) )
  write.xlsx(enrichment_df,here::here("Datasets","Processed","KEGG-Multi-Omic_reduced.xlsx"), overwrite =T)
  enrichment_df <- openxlsx::read.xlsx(here::here("Datasets","Processed","KEGG-Multi-Omic_reduced.xlsx"))
  
  Interesting_terms <- rbind(subset(enrichment_df, n== 1 & `p.adjust`<0.05 ) %>% arrange(`p.adjust`) %>% 
                                 group_by(comparison) %>%
                                 slice(1:5),subset(enrichment_df, n> 1 & `p.adjust`<0.05 ) %>% arrange(`p.adjust`) %>% 
                                 group_by(comparison) %>%
                                 slice(1:5)) %>% pull(Description)
  enrichment_df <- subset(enrichment_df,Description %in% Interesting_terms)%>% 
      # subset(comparison != "dmso_vs_x24h_diff") %>% 
      mutate(NES = if_else(`p.adjust`<0.05,NES,0),
             comparison= if_else(comparison == "D14_BP","D14_BP_rev",comparison)) %>% 
      # mutate(NES = if_else(str_detect(comparison,"vs"),NES,-NES)) %>% 
      arrange(n,comparison) %>% 
      mutate(Description  = factor(Description, Description %>% unique))
  ggplot(enrichment_df %>% #subset(comparison != "D10_high_BP")
             mutate(NES = if_else(str_detect(comparison,"_vs"),-NES,NES))
         ,aes(x = comparison,y  = reorder(Description,n), size = -log10(p.adjust), colour = NES))+
      geom_point()+
      theme_bw()+
      scale_colour_gradient2(high= muted("red"), low = muted("blue"),mid = "grey70")+
      scale_x_discrete(
          labels=c("x0h_vs_x24h_diff" = "Release vs No release",
                   "dmso_vs_x0h_diff" = "No release vs Untreated",
                   "dmso_vs_x24h_diff" = "Release vs Untreated",
                   "D10_high_BP" = "gH2X-High vs Unsorted",
                   "D14_BP_rev" = "D5 Etop vs DMSO"))+
      # guide = guide_axis(n.dodge = 2))+
      ggtitle("Biological Response to Etoposide Release")+
      theme(axis.text.y = element_text(size = 15),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 15))
  ggsave(here::here("Output","Figures","Multi-Omic-KEGGenrich_reduced.pdf"), height = 15, width = 10)
  # enrichment_0_24 <- enrichment$x0h_vs_x24h_diff %>%clusterProfiler::simplify()
  # Simplified_terms <- enrichment_0_24@result$Description
  # enrichment_dmso_etop <- enrichment$dmso_vs_x0h_diff#  %>% simplify()
  # enrichment_df <- enrichment %>% imap_dfr(.x = ., ~.x@result %>% mutate(comparison = .y)) %>% 
  #   subset(Description %in%Simplified_terms)
  # enrichment_df <- enrichment %>% imap_dfr(.x = ., ~.x@result%>% as.data.frame %>% 
  #                                            dplyr::select(Description,NES,p.adjust)  %>% mutate(comparison = .y)) %>% 
  #   mutate(Significant = p.adjust<0.05) %>% 
  #   group_by(Description) %>% 
  #   add_tally(Significant) %>% subset(n>0)%>% 
  #   subset(Description %in%Simplified_terms) %>% 
  #   mutate(Type = case_when(
  #     p.adjust<0.05 & NES>0 ~ "Signi_Up",
  #     p.adjust<0.05 & NES<0 ~ "Signi_Down",
  #     p.adjust>0.05~ "Not_Signi",
  #     TRUE ~"Error"
  #   )) %>%dplyr::select(comparison,Type,Description)  %>% 
  #   pivot_wider(names_from = "comparison", values_from = "Type" ) %>% 
  #   ungroup %>% 
  #   
  #   mutate(concat = paste0(dmso_vs_x0h_diff,dmso_vs_x24h_diff,x0h_vs_x24h_diff),
  #          Behaviour = if_else(str_detect(concat,"Up") & (str_detect(concat,"Down")|str_detect(concat,"Not")),"Unconsistent","Consistent" )) %>% 
  #   dplyr::select(-c(concat,Description)) %>%
  #   group_by(dmso_vs_x0h_diff,dmso_vs_x24h_diff,x0h_vs_x24h_diff,Behaviour) %>% 
  #   count()
  # 
  # pacman::p_load(ggalluvial)
  # ggplot(as.data.frame(enrichment_df),
  #        aes(y = freq, axis1 = dmso_vs_x0h_diff, axis2 = x0h_vs_x24h_diff, axis3 = dmso_vs_x24h_diff)) +
  #   geom_alluvium(aes(fill = Behaviour), width = 1/12) +
  #   geom_stratum(width = 1/12, fill = "grey60", color = "grey") +
  #   geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  #   scale_x_discrete(limits = c("dmso_vs_x0h_diff", "x0h_vs_x24h_diff","dmso_vs_x24h_diff"), expand = c(.05, .05)) +
  #   scale_fill_manual(values = c(Consistent  ="red",
  #                                Signi_Down = "grey20"))+
  #   theme_bw()+
  #   ggtitle("GSEA BP terms Proteomic Comparison")
  # terms <- c("nucleobase-containing small molecule metabolic process",
  #            "ATP metabolic process","response to drug","response to oxidative stress","lipid metabolic process","regulation of programmed cell death",
  #            "mitochondrial gene expression","proteolysis","nuclear division","cell cycle phase transition","DNA metabolic process","protein acylation",
  #            "DNA integrity checkpoint","regulation of response to DNA damage stimulus","chromatin organization involved in regulation of transcription")
  # new_terms <- c("purine-containing compound metabolic process","monocarboxylic acid metabolic process","autophagy","mitotic sister chromatid segregation",
  #                "proteasomal protein catabolic process","regulation of signal transduction",
  #                "carbohydrate metabolic process","neutrophil degranulation","cell division","RNA splicing",
  #                "microtubule cytoskeleton organization","histone modification","cellular response to DNA damage stimulus")
  # enrichment_0_24@result <-  enrichment_0_24@result %>% subset(Description %in% terms & p.adjust<0.05)
  # enrichment_dmso_etop@result <- enrichment_dmso_etop@result %>% subset(Description %in% new_terms & p.adjust<0.05)
  # enrichment_0_24@result %>% ggplot(aes(x = NES, y = reorder(Description,NES), colour = NES, size = -log10(p.adjust)))+geom_point()+
  #   guides(size = 'none')+
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #         # axis.ticks.x = element_blank(),
  #         # axis.text.x = element_blank(),
  #         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #   scale_color_gradient2(low = "#125CC7", high = "#FF2400")+
  #   labs(x= "Enrichement Score", y = "Categories Enriched")+
  #   ggtitle("Processes Perturbed by Etopside Release")
  # enrichment_dmso_etop@result %>% ggplot(aes(x = NES, y = reorder(Description,NES), colour = NES, size = -log10(p.adjust)))+geom_point()+
  #   guides(size = 'none')+
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #         # axis.ticks.x = element_blank(),
  #         # axis.text.x = element_blank(),
  #         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #   scale_color_gradient2(low = "#125CC7", high = "#FF2400")+
  #   labs(x= "Enrichement Score", y = "Categories Enriched")+
  #   ggtitle("Processes Perturbed by Etopside Treatment")
  # 
  # enrichplot::ridgeplot(testing,showCategory = 100,core_enrichment = F)
  # Simplified_terms <- enrichment$dmso_vs_x0h_diff %>% simplify()
  # significant_terms <-enrichment_df %>% subset(p.adjust <0.05) %>% pull(ID) 
  # terms_to_use = c("cell division","negative regulation of nucleobase-containing compound metabolic process",
  #                  "secretion","cation transport","biological adhesion","regulation of translation","amide transport",
  #                  "autophagy","histone modification","ATP metabolic process","oxoacid metabolic process","cellular response to DNA damage stimulus",
  #                  "response to oxidative stress","response to cytokine","DNA integrity checkpoint" ,"response to drug",
  #                  "regulation of carbohydrate metabolic process","lipid metabolic process","small molecule metabolic process" ,
  #                  "ATP metabolic process"  ,"phosphorus metabolic process"  
  #                  
  #                  )
  # plotting_enrichment <- enrichment_df %>% 
  #   subset(ID %in% significant_terms) %>% 
  #   subset(p.adjust <0.05 & Description %in% terms_to_use) %>% 
  #   dplyr::select(ID,Description,NES,comparison,p.adjust) %>% 
  #   group_by(Description) %>% 
  #   mutate(N = n(),
  #          NES_med= median(NES)) %>% 
  #   ungroup() %>% arrange(N,NES_med) %>% 
  #   mutate(Description = factor(Description, levels = unique(Description)))
  #   #,pvalue,p.adjust
  #   # pivot_wider(names_from = comparison, values_from = NES) %>% 
  #   # mutate(label = if_else(str_detect(Description, "damage|metabo|RNA|chrom|nucleo"),T,F))
  # 
  # plotting_enrichment %>% 
  #   ggplot(aes(x = comparison,y  = Description, size = -log10(p.adjust), colour = NES))+
  #   geom_point()+
  #   theme_bw()+
  #   scale_colour_gradient2(high= muted("red"), low = muted("blue"))+
  #   scale_x_discrete(labels=c("x0h_vs_x24h_diff" = "Etop_vs_Release",
  #                             "dmso_vs_x0h_diff" = "DMSO_vs_Etop",
  #                               "dmso_vs_x24h_diff" = "DMSO_vs_Release"))+
  #   ggtitle("Biological Response to Etoposide Release")
  # 
  # 
  # dmso_x0 <- enrichment$dmso_vs_x24h_diff %>% clusterProfiler::simplify()
  # dmso_x0@result <- dmso_x0@result %>% 
  #   subset(str_detect(Description,"repair|damage|metabol|chrom") & p.adjust<0.05)
  # enrichplot::ridgeplot(dmso_x0, showCategory = 68)+
  #   ggtitle("Gene Enrichment between DMSO and 24Hr Release")
  # 
  # #Diff Corr#
  # design_mat = data.frame( X0Hrs  = c(rep(1,6), rep(0,6)),
  #                          X24Hrs = c(rep(0,6),rep(1,6))) %>% as.matrix()
  # ddcor_res = DGCA::ddcorAll(inputMat = Data_matrices$Imputted %>% dplyr::select(!contains("dmso")) %>% 
  #                              .[,colnames(.) %>% sort()], design = design_mat,
  #                      compare = c("X0Hrs", "X24Hrs"),corrType = "pearson",dCorAvgMethod = "mean",
  #                      adjust = "none", nPerm = 0, nPairs = 100)
  # head(ddcor_res)
  # 
  # ddcor_res_genes <- Data_matrices$Imputted %>% dplyr::select(!contains("dmso")) %>%
  #   mutate(Rowmean = rowMeans(.),
  #          across(everything(),~.x- Rowmean)) %>% 
  #   dplyr::select(-Rowmean) %>% 
  #   rownames_to_column("Uniprot") %>% 
  #   
  #   pivot_longer(-Uniprot, names_to = "Condition", values_to = "Abundance")%>% 
  #   separate(Condition, c("Condition","Replicate")) %>%
  #   mutate(Single_Uniprot = str_remove_all(Uniprot,";[:graph:]*$")) %>% 
  #   left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% subset(!duplicated(Uniprot)), by = c("Single_Uniprot" = "Uniprot")) %>% 
  #   subset(Uniprot %in% (c(ddcor_res$Gene1,ddcor_res$Gene2) %>% unique()))
  # ddcor_res_genes %>% 
  #   subset(Uniprot %in% c("O15511;Q9BPX5","P45880")) %>% 
  #   ggplot(aes(x = Replicate, y = Abundance, color = ID, group = ID))   +
  #   geom_point(alpha = 0.5)+
  #   geom_line(alpha = 0.5)+
  #   facet_wrap(~Condition)+
  #   ggtitle("Protein Behaviour Etop Release")
  # 
  # 
  # Temp_proteomic_ruler %>%  
  #   group_by(Uniprot,ID,Compartments ) %>%
  #   dplyr::summarise(Mean_Abundance_u2os = mean(Abundance, na.rm = T),.groups = "keep") %>% 
  #   ungroup() %>% 
  #   left_join(Temp_proteomic_ruler_sample_2 %>% 
  #               group_by(Uniprot,ID,Compartments ) %>%
  #               dplyr::summarise(Mean_Abundance_sample = mean(Abundance, na.rm = T),.groups = "keep") %>% 
  #               ungroup()) %>% subset(!is.na(Compartments)) %>% 
  #   ggplot(aes (x = reorder(Uniprot,Mean_Abundance_sample),y = Mean_Abundance_sample, fill = Compartments))+
  #   geom_col(width = 1)+  facet_grid(~Compartments, scales = "free_x",space = "free_x")+
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x = element_blank(),
  #         axis.text.x = element_blank(),
  #         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #   scale_fill_manual(values = colours_compartmens)+
  #   labs(x= "Proteins Identified on Chromatin",y ="Copy Numbers on Chromatin" )+
  #   ggtitle("Protein Copy Numbers Detected on Chromatin")
  # 
  # Data_matrices$Unimputted %>% 
  #   rownames_to_column("Uniprot") %>% 
  #   mutate(Uniprot = str_remove_all(Uniprot,";[:graph:]*$")) %>% 
  #   left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% 
  #               subset(!duplicated(Uniprot)), by =  "Uniprot")%>%
  #   pivot_longer(cols = contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
  #   group_by(ID) %>% summarise(Median_Abundance = median(Abundance, na.rm = T)) %>% 
  #    left_join(markerProteins[,1:2], by = c("ID" = "Proteins"))  %>%  subset(!is.na(Compartments)) %>% 
  #   ggplot(aes (x = reorder(ID,Median_Abundance),y = Median_Abundance, fill = Compartments))+
  #   geom_col(width = 1)+  facet_grid(~Compartments, scales = "free_x",space = "free_x")+
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x = element_blank(),
  #         axis.text.x = element_blank(),
  #         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #   scale_fill_manual(values = colours_compartmens)+
  #   labs(x= "Proteins Identified on Chromatin",y ="LFQ intensity on Chromatin" )+
  #   ggtitle("LFQ intensity on Chromatin")
  # 
  # 
  # ggplot(Temp_proteomic_ruler %>% subset(Condition == "dmso_1") %>% 
  #          mutate(
  #            # Compartment = case_when(
  #          #   Compartments %in% c("N1","N2","N4") ~"Nuclear",
  #          #   Compartments %in% c("N3") ~"Chromatin",
  #          #   str_detect(Compartments, "^C") ~"Cytoplasm",
  #          #   str_detect(Compartments, "^S") ~"Secretory",
  #          #   str_detect(Compartments, "^M") ~"Mitochondria",
  #          # ),
  #          Abundance = scale(Abundance)[,1]) %>% 
  #          subset(!is.na(Compartments)),
  #        aes(x = reorder(PG, Abundance), y= Abundance, fill = Compartments))+
  #   geom_col(width= 1)+
  #   theme(axis.ticks.x = element_blank(),
  #         axis.text.x = element_blank(), legend.position = "none")+
  #   facet_grid(~Compartments, scales = "free_x",space = "free_x")+labs(x = "Proteins with known localisation", y = "Relative Enrichemnt" )+
  #   scale_fill_manual(values = colours_compartmens)+
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x = element_blank(),
  #         axis.text.x = element_blank(),
  #         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #   ggtitle("Protein Pool Percentage on Chromatin")
  # ###Combining Enrichment
  # simplified_terms <- enrichment$dmso_vs_x24h_diff %>%clusterProfiler::simplify()
  # simplified_terms_screen <- ego3 %>%clusterProfiler::simplify()
  # FELLA_result = data.frame(Description = c("nucleotide-excision repair",
  #                                           "DNA replication"),
  #                           p.adjust = 0.004,
  #                           NES = NA,
  #                           Omic = "Metabolomics")
  # multi_omic_BP <- rbind(enrichment$dmso_vs_x24h_diff@result %>% mutate(
  #   Omic = "Chromatin Proteomics",
  #   NES = NES*(-1)),
  #   ego3@result%>% mutate(
  #     Omic = "Screening Day 14")) %>% 
  #   dplyr::select(Description,NES,p.adjust,Omic) %>% 
  #   # rbind(FELLA_result
  #   # )  %>% 
  #   subset(p.adjust <0.05) %>% 
  #   # dplyr::select(ID,Description,NES,comparison,p.adjust) %>% 
  #   group_by(Description) %>% 
  #   mutate(N = n(),
  #          omic = case_when(Omic== "Chromatin Proteomics" ~ 1,
  #                           Omic == "Screening Day 14"~2)) %>% 
  #   ungroup() %>% arrange(N,omic) %>% 
  #   mutate(Description = factor(Description, levels = unique(Description)))
  # terms_to_plot <- c(intersect(multi_omic_BP %>% subset(N>1) %>% pull(Description),simplified_terms@result$Description),
  #                    intersect(multi_omic_BP %>% subset(N==1) %>% pull(Description),simplified_terms@result %>% 
  #                                arrange(p.adjust) %>% head(8) %>% pull(Description)),
  #                    intersect(multi_omic_BP %>% subset(N==1) %>% pull(Description),simplified_terms_screen@result%>% 
  #                                arrange(p.adjust) %>% head(8) %>% pull(Description))) %>% unique()
  # 
  # terms_to_plot <- c(FELLA_result$Description)
  # multi_omic_BP %>%    
  #   subset(Description %in% terms_to_plot) %>% 
  #   ggplot(aes(x = Omic,y  = Description, size = -log10(p.adjust), colour = NES))+
  #   geom_point()+
  #   theme_bw()+
  #   scale_colour_gradient2(high= muted("red"), low = muted("blue"))+
  #   # scale_x_discrete(labels=c("x0h_vs_x24h_diff" = "Etop_vs_Release",
  #   #                           "dmso_vs_x0h_diff" = "DMSO_vs_Etop",
  #   #                           "dmso_vs_x24h_diff" = "DMSO_vs_Release"))+
  #   ggtitle("Biological Pathways Across Omics")
  # 
  # ##### Metabolic terms #
  # FELLA_result_KEGG = rbind(data.frame(Description = c("Carbon metabolism",
  #                                           "Insulin signalling pathway"),
  #                           p.adjust = 0.004,
  #                           NES = NA,
  #                           Omic = "Metabolomics"),
  #                           data.frame(Description = c("Propanoate Metabolism",
  #                                                      "Pentose Phosphate Pathway",
  #                                                      "Purine catabolism",
  #                                                      "Glycogenesis",
  #                                                      "Transport peroxisomal",
  #                                                      "Urea Cycle",
  #                                                      "Transport extracellular"),
  #                                      p.adjust = 0.004,
  #                                      NES = c(-3.5,-2.7,-2.6,-2.1,-1.2,1.2,1.6),
  #                                      Omic = "Metabolomics"))
  # 
  # 
  # 
  # gene_list <- Volcano_DFs$dmso_vs_x24h_diff %>%
  #   dplyr::select(Uniprot,log2_FC) %>%
  #   # mutate(log2_FC = log2_FC*(-1)) %>% 
  #   left_join(Human_hsa) %>%
  #   na.omit()%>%arrange(-log2_FC) %>%   pull(log2_FC,ID)
  # 
  # dmdo_24 <- gseKEGG(geneList     = gene_list,
  #                organism     = 'hsa',
  #                minGSSize    = 10,
  #                pvalueCutoff = 0.05,
  #                verbose      = FALSE)
  # multi_omic_KEGG <- rbind(dmdo_24@result %>% mutate(
  #   Omic = "Chromatin Proteomics"),
  #   kk2@result%>% mutate(
  #     Omic = "Screening Day 14"))  %>% 
  #   subset(p.adjust <0.05) %>% 
  #   dplyr::select(Description,NES,p.adjust,Omic) %>% 
  # 
  #   group_by(Description) %>% 
  #   mutate(N = n(),
  #          omic = case_when(Omic== "Chromatin Proteomics" ~ 1,
  #                           Omic == "Screening Day 14"~2,
  #                           Omic == "Metabolomics" ~ 3)) %>% 
  #   ungroup() %>% arrange(N,omic) %>% 
  #   mutate(Description = factor(Description, levels = unique(Description)))
  # terms_to_plot <- list(multi_omic_KEGG %>% subset(N>1) %>% pull(Description) %>% as.character() %>% 
  #                         str_subset("Carbon|Oxida|Chemical"),
  #                    multi_omic_KEGG %>% mutate(Description= as.character(Description)) %>%  subset(omic == 1) %>% slice_max(order_by = p.adjust,n= 5) %>% pull(Description),
  #                    multi_omic_KEGG %>%  mutate(Description= as.character(Description)) %>%subset(omic == 2) %>% slice_max(order_by = p.adjust,n= 5) %>% pull(Description),
  #                    multi_omic_KEGG %>%  mutate(Description= as.character(Description)) %>%subset(omic == 3) %>% slice_max(order_by = p.adjust,n= 5) %>% pull(Description)) %>% unlist %>% unique()
  # 
  # 
  # multi_omic_KEGG %>%    
  #    subset(Description %in% terms_to_plot) %>% 
  #   ggplot(aes(x = Omic,y  = Description, size = -log10(p.adjust), colour = NES))+
  #   geom_point()+
  #   theme_bw()+
  #   scale_colour_gradient2(high= muted("red"), low = muted("blue"))+
  #   # scale_x_discrete(labels=c("x0h_vs_x24h_diff" = "Etop_vs_Release",
  #   #                           "dmso_vs_x0h_diff" = "DMSO_vs_Etop",
  #   #                           "dmso_vs_x24h_diff" = "DMSO_vs_Release"))+
  #   ggtitle("KEGG Across Omics")
  # 
  # kk2@result
  
  Metabolic_genes <- subset(HUMAN_9606,Uniprot %in%  Metabolic_proteins$Uniprot & Type =="Gene_Name") %>% pull(ID) %>% unique()
  pacman::p_load(eulerr)
  s4 <- list(Chromatin_Detected = Chromatin_genes %>% unique(),
             Metabolic_Library = Sabatini,
             Human_GEM =Human_GEM_genes ,
             CRISPR_hits = All_interesting_genes %>% unique(),
             Etop_Treatment_Chrom=Significant_genes%>% unique() )
  write.xlsx(s4,here::here("Datasets","Processed","Venn_Overlap.xlsx"))
  pdf(here::here("Output","Figures","Venn_GEM_Chrom_Metablolism.pdf"))
  plot(venn(s4))
  dev.off()
  pdf(here::here("Output","Figures","Euler_GEM_Chrom_Metablolism.pdf"))
  plot(euler(s4, shape = "ellipse"), quantities = TRUE)
  dev.off()
  pdf(here::here("Output","Figures","Venn_upset_plot.pdf"),width = 15, height =10)
  UpSetR::upset(fromList(s4), order.by = "freq")
  dev.off()
  
  ##### high and low hits ####
  hits_comparison <- list(D10_enriched = Screen_input_volc %>% subset(Significant  =="High") %>% pull(Gene),
                          D10_depleted= Screen_input_volc %>% subset(Significant  =="Low") %>% pull(Gene),
                          D14_depleted = subset(d14_genes,group =="bottom") %>% pull(Gene),
                          D14_enriched = subset(d14_genes,group =="top") %>% pull(Gene)) %>% map(as.character)
  enrich_clusters <- function(gene_list){
      enrichGO(gene         = gene_list,
               OrgDb         = org.Hs.eg.db,
               universe = d14_genes$Gene%>% unique(),
               keyType       = 'SYMBOL',
               ont           = "BP",
               # minGSSize = 5,
               # maxGSSize = 50,
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.01,
               qvalueCutoff  = 0.05)
  }
  ck <- compareCluster(geneCluster = hits_comparison, fun = enrich_clusters )
  
  ck@compareClusterResult %>% 
      mutate(GeneRatio = sapply(strsplit(GeneRatio, split = "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))) %>%
      group_by(Cluster) %>% 
      arrange(p.adjust) %>% 
      slice(1:20) %>% ungroup %>% 
      arrange(Cluster) %>% mutate(Description = factor(Description, levels = rev(unique(Description)))) %>% 
      ggplot(aes(x = Cluster, y= Description, colour = p.adjust, size = GeneRatio))+
      geom_point()+
      ggtitle("Screens Enriched-Depleted Enrichment BP",
              subtitle = "All terms in each timepoint share basically all genes \n D10 = DTYMK/DUT/TXNRD1/TYMS/RRM1 \n D14 = mitochondrial/phosphorylation \n No enriched terms for Depleted D10 and Depleted D14")
  ggsave(here::here("Output","Figures", "Screens Enriched-Depleted Enrichment BP.pdf"),height = 10, width = 15)
  
  ####PRDX1 PPI Cytoscape ####
  PRDX1_int <- fread(here::here("Datasets","Processed", "IntAct Network - PRDX1 default  node.csv"))
  PRDX1_int_enrich <- enrichGO(gene         = PRDX1_int$`shared name`,
                               OrgDb         = org.Hs.eg.db,
                               universe = c(interactions$source_genesymbol,interactions$target_genesymbol) %>% unique(),
                               keyType       = 'SYMBOL',
                               ont           = "BP",
                               # minGSSize = 5,
                               # maxGSSize = 50,
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 0.05)
  terms_to_reduce <- PRDX1_int_enrich@result %>% subset(`p.adjust`<0.05) %>% 
      # dplyr::select(core_enrichment,Description) %>%
      pull(geneID,Description)  %>% purrr::map(str_split,"/") %>% flatten()
  reduced_terms <- Jaccard_index_list(terms_to_reduce,max_jacc = 0.3)
  PRDX1_int_enrich_plot <- PRDX1_int_enrich@result %>% subset(Description %in% (reduced_terms %>% names()) ) %>% 
      separate(GeneRatio,into = c("num","denom")) %>% mutate(
          GeneRatio = as.numeric(num)/as.numeric(denom)
      )
  write.xlsx(PRDX1_int_enrich@result,here::here("Datasets","Processed","PPI_PRDX1_ORA.xlsx"))  
  ggplot(PRDX1_int_enrich_plot,aes(x = GeneRatio,y  = reorder(Description,GeneRatio), size = Count, colour = `p.adjust`))+
      geom_point()+
      theme_bw()+
      # scale_colour_gradient2(high= muted("red"), low = muted("blue"))+
      ggtitle("PRDX1 PPI")
  ggsave(here::here("Output","Figures","PRDX1_PPIs_enrichemt.pdf"),width = 10,height = 15)
  
  #Proteins correlated with PRDX1 behaviour
  data_matrix = Data_matrices$Imputted %>% dplyr::select(matches("_(1|2|3)$"))
  subset_selection = rownames(data_matrix)
  pacman::p_load(e1071,cluster,CORREP,NbClust,advclust) #
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
  
  
  
  Correlated_proteins <- M["Q06830",] %>% sort() %>% tail(15) %>% enframe(name = "Uniprot")
  Correlated_proteins_IMPDH2 <- M["P12268",] %>% sort() %>% tail(163) %>% enframe(name = "Uniprot")
  
  Correlated_proteins_IMPDH2 %>% 
      mutate(Single_Uniprot = str_remove_all(Uniprot,";[:graph:]*")) %>% 
      left_join(Metabolic_proteins, by = c("Single_Uniprot"= "Uniprot" )) %>% 
      left_join(subset(HUMAN_9606, Type == "Gene_Name"),by = c("Single_Uniprot"= "Uniprot" ) ) %>% 
      mutate(ID_behaviour = if_else(is.na(Behaviour), NA_character_,ID)) %>% 
      ggplot(aes( x= rank(value), y = value, colour = Behaviour, label = ID_behaviour ))+
      geom_point()+
      ggrepel::geom_text_repel( max.overlaps = 100) +
      lims(y = c(0,1))+
      ggtitle("Proteins most correlated with IMPDH2 behaviour")
  
  Correlated_proteins_IMPDH2 %>% 
      mutate(Single_Uniprot = str_remove_all(Uniprot,";[:graph:]*")) %>% 
      left_join(Metabolic_proteins, by = c("Single_Uniprot"= "Uniprot" )) %>% 
      left_join(subset(HUMAN_9606, Type == "Gene_Name"),by = c("Single_Uniprot"= "Uniprot" ) ) %>% 
      mutate(ID_behaviour = if_else(!is.na(Behaviour), NA_character_,ID)) %>% 
      ggplot(aes( x= rank(value), y = value, colour = Behaviour, label = ID_behaviour ))+
      geom_point()+
      ggrepel::geom_text_repel( max.overlaps = 100) +
      lims(y = c(0,1))+
      ggtitle("Proteins most correlated with IMPDH2 behaviour")
  
  IMPDH2_cor=M["P12268",] %>% sort(decreasing = T) 
  names(IMPDH2_cor) <-  names(IMPDH2_cor) %>% str_remove_all(";[:graph:]*")
  IMPDH2_cor <- IMPDH2_cor[!(IMPDH2_cor %>% names() %>% duplicated())]
  IMPDH2_correlated <- clusterProfiler::gseGO(IMPDH2_cor,OrgDb = org.Hs.eg.db,
                                              
                                              ont          = "BP",
                                              minGSSize    = 100,
                                              keyType = "UNIPROT",
                                              maxGSSize    = 500,
                                              pvalueCutoff = 0.05,
                                              eps = 0,
                                              verbose      = FALSE)
  IMPDH2_correlated_simpl <- simplify(IMPDH2_correlated)
  ridgeplot(IMPDH2_correlated_simpl)+
      ggtitle("GSEA on the IMPDH2 correlation scores of proteins, +ve = correlated, -ve = lowly/anti-correlated")
  
  ggsave(here::here("Output","Figures", 'Proteins_correlated_IMPDH2_gsea.pdf') ) 
  Data_matrices$Imputted %>% 
      subset(rownames(.)%in%Correlated_proteins_IMPDH2$Uniprot) %>% as.data.frame() %>% 
      rownames_to_column("Uniprot") %>% 
      pivot_longer(-Uniprot,names_to = "Sample",values_to = "Abundance") %>% 
      mutate(Condition = str_remove_all(Sample,"_.$"),
             Replicate = str_remove_all(Sample,"[:graph:]*_"))%>%
      inner_join(Correlated_proteins_IMPDH2) %>% 
      mutate(Uniprot = str_remove_all(Uniprot,";[:graph:]*")) %>%
      left_join(subset(HUMAN_9606, Type == "Gene_Name") ) %>% 
      arrange(-value) %>% 
      mutate(ID = factor(ID, levels = unique(ID))) %>% 
      ggplot(aes(x = Condition, y = Abundance, group = Uniprot))+
      geom_point(alpha = 0.2)+theme_bw()+
      # annotate("text",label = "PRDX1", x="x24h", y=20.25 ,size = 10, colour = muted("red"))+
      geom_smooth(fill = muted("red"),
                  colour = muted("red"), alpha = 0.1)+
      scale_x_discrete(labels= c("DMSO","No Release","Release"))+
      theme(axis.text=element_text(size=12),strip.text = element_text(size=18),
            axis.title=element_text(size=12,face="bold"))+
      ggtitle("Proteins most correlated with IMPDH2 behaviour")+facet_wrap("ID")
  ggsave(here::here("Output","Figures", 'Proteins_correlated_PRDX1.pdf') ) 
  
  Complexes_mito <- read.delim(here::here("Datasets","Raw","allComplexes.txt")) %>% dplyr::select(`Approved.symbol`,`Group.name`) %>% 
      dplyr::rename(Genes =`Approved.symbol`,
                    Complex =`Group.name` )
  
  Complexes_mito_proteins <- subset(HUMAN_9606, ID %in% Complexes_mito$Genes & Type == "Gene_Name") %>%
      pull(Uniprot) %>% unique()
  Data_matrices$Imputted %>% 
      subset(str_detect(rownames(.), paste0(Complexes_mito_proteins, collapse = "|"))) %>% as.data.frame() %>% 
      rownames_to_column("Uniprot") %>% 
      pivot_longer(-Uniprot,names_to = "Sample",values_to = "Abundance") %>% 
      mutate(Condition = str_remove_all(Sample,"_.$"),
             Replicate = str_remove_all(Sample,"[:graph:]*_"))%>%
      mutate(Uniprot = str_remove_all(Uniprot,";[:graph:]*")) %>%
      left_join(subset(HUMAN_9606, Type == "Gene_Name") ) %>% 
      # arrange(-value) %>% 
      mutate(ID = factor(ID, levels = unique(ID))) %>% 
      ggplot(aes(x = Condition, y = Abundance, group = Uniprot))+
      geom_point(alpha = 0.2)+theme_bw()+
      # annotate("text",label = "PRDX1", x="x24h", y=20.25 ,size = 10, colour = muted("red"))+
      geom_smooth(fill = muted("red"),
                  colour = muted("red"), alpha = 0.1)+
      scale_x_discrete(labels= c("DMSO","No Release","Release"))+
      theme(axis.text=element_text(size=12),strip.text = element_text(size=18),
            axis.title=element_text(size=12,face="bold"))+
      ggtitle("Proteins most correlated with PRDX1 behaviour")+facet_wrap("ID")
  ggsave(here::here("Output","Figures", 'ETC_Proteins.pdf'), width = 20, height = 20 ) 
  Data_matrices$Unimputted %>% 
      subset(str_detect(rownames(.), paste0(Complexes_mito_proteins, collapse = "|"))) %>% as.data.frame() %>% 
      rownames_to_column("Uniprot") %>% 
      pivot_longer(-Uniprot,names_to = "Sample",values_to = "Abundance") %>% 
      mutate(Condition = str_remove_all(Sample,"_.$"),
             Replicate = str_remove_all(Sample,"[:graph:]*_"))%>%
      mutate(Uniprot = str_remove_all(Uniprot,";[:graph:]*")) %>%
      left_join(subset(HUMAN_9606, Type == "Gene_Name") ) %>% 
      # arrange(-value) %>% 
      mutate(ID = factor(ID, levels = unique(ID))) %>% 
      ggplot(aes(x = Condition, y = Abundance, group = Uniprot))+
      geom_point(alpha = 0.2)+theme_bw()+
      # annotate("text",label = "PRDX1", x="x24h", y=20.25 ,size = 10, colour = muted("red"))+
      geom_smooth(fill = muted("red"),
                  colour = muted("red"), alpha = 0.1)+
      scale_x_discrete(labels= c("DMSO","No Release","Release"))+
      theme(axis.text=element_text(size=12),strip.text = element_text(size=18),
            axis.title=element_text(size=12,face="bold"))+
      ggtitle("Proteins most correlated with PRDX1 behaviour")+facet_wrap("ID")
  ggsave(here::here("Output","Figures", 'ETC_Proteins_Unimputted.pdf'), width = 20, height = 20 ) 
  Data_matrices$Unimputted %>% 
      subset(str_detect(rownames(.), paste0(Complexes_mito_proteins, collapse = "|"))) %>% as.data.frame() %>% 
      mutate(across(everything(),~if_else(is.na(.x),0,1))) %>% pheatmap(cluster_cols = F, main= "Presence Absence Mito complexes")
  
  
  Volcano_DFs$x0h_vs_x24h_diff %>% 
      mutate(ETC = if_else(Single_Uniprot %in% Complexes_mito_proteins,T,F)) %>%
      ggplot(aes(x = -log2_FC,y = -log10(p.val), colour = ETC, label = ID))+
      geom_point(data = . %>% subset(ETC ==F))+theme_bw()+
      geom_point(data = . %>% subset(ETC ==T))+
      ggrepel::geom_text_repel(data = . %>% subset(ETC == T), max.overlaps =  20)+
      ggtitle("24hrs release vs No release")+
      annotate("text", x = c(-1.5,1.5), y=0, label = rev(c("Release 24Hrs","Etop 3hrs")))
  ggsave(here::here("Output","Figures", 'ETC_Proteins_volcano_24_T0.pdf'), width = 20, height = 20 ) 
  
  Volcano_DFs$dmso_vs_x0h_diff %>% 
      mutate(ETC = if_else(Single_Uniprot %in% Complexes_mito_proteins,T,F)) %>%
      ggplot(aes(x = -log2_FC,y = -log10(p.val), colour = ETC, label = ID))+
      geom_point(data = . %>% subset(ETC ==F))+theme_bw()+
      geom_point(data = . %>% subset(ETC ==T))+
      ggrepel::geom_text_repel(data = . %>% subset(ETC == T), max.overlaps =  20)+
      ggtitle("No release vs DMSO")+
      annotate("text", x = c(-1.5,1.5), y=0, label = c("DMSO","Etop 3hrs"))
  ggsave(here::here("Output","Figures", 'ETC_Proteins_volcano_T0_DMSO.pdf'), width = 20, height = 20 ) 
  