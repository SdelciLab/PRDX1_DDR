install.packages('pacman')
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("EnrichmentBrowser")
BiocManager::install("clusterProfiler")
BiocManager::install("seqinr")
BiocManager::install("org.Hs.eg.db")
pacman::p_load(clusterProfiler,org.Hs.eg.db,'here','readxl',
               tidyverse,data.table)
#Figure 1 
kegg_genes <-  EnrichmentBrowser::getGenesets("hsa",db  = "kegg")
gs <- EnrichmentBrowser::getGenesets("hsa")

HUMAN_9606 <- read_tsv(here::here("Datasets","Raw", "HUMAN_9606_idmapping.dat"),
                       col_names = FALSE) %>% set_names(c("Uniprot","Type","ID"))
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
go_KEGG <- purrr::map_dfr(.x = list(kegg_genes,gs),~
                              .x %>% unlist() %>% enframe(name = "pathway",
                                                          value = "KEGG") %>% 
                              # mutate(pathway = str_remove_all(pathway,"_[:graph:]*$")) %>% 
                              separate(pathway,into = c("id","name"),sep = "_",extra = "merge") %>% 
                              left_join(Human_hsa))
fwrite(go_KEGG, here::here("Datasets","Processed","go_KEGG.csv"))

#these are the essential genes as classified by MageckFlute publication in 
#https://www.bioconductor.org/packages/release/bioc/html/MAGeCKFlute.html
mageckflute_essential <- readxl::read_xls(here::here("Datasets","Raw","Mageck_flute_essential_625.xls"))[-1,2] %>% unlist()

convert_KEGG_into_Symbol <- function(x){
    # x <- "163/476/483/481/476/1213/1175/160/161/5566/1211"
    subset(Human_hsa, KEGG %in% ( x %>% str_split("/",simplify = T))) %>% pull(ID) %>% unique()
}
core_enrichment_KEGG <- openxlsx::read.xlsx(here::here("Datasets","Processed","KEGG-Multi-Omic_reduced.xlsx")) %>% 
    mutate(Genes = purrr::map(core_enrichment,convert_KEGG_into_Symbol))
Selection_dirs <- list.dirs(here::here("Datasets","Processed","Amandine_screen")) %>% str_subset("Selection") %>% 
    set_names(.,purrr::map_chr(.x =., ~ str_match(.x,"MAGeCKFlute_([:print:]*)_") %>% .[,2]))

hits <- purrr::imap_dfr(.x = Selection_dirs, 
                        ~list.files(.x) %>% str_subset("squareview_data_fix") %>% here::here(.x,.) %>% fread() %>% 
                            .[!(Label == ""),.(Gene,Diff,Depmap)] %>% mutate(Comparison = .y)
) 
#write_tsv(hits, here::here("Datasets","Processed","Amandine_screen","cell_cycle_hits.tsv"))

conditions <- c("D10_D10_etop",
                #"D10etophighlow",
                "D14_D14_etop")
All_interesting_genes <- high_low_hits %>% as.character()
list_of_comp <- list()
for(condition in conditions){
    condition <- conditions[2]
    Interesting_genes <- list.files(Selection_dirs[condition]) %>% str_subset("^Data_ScatterView_TreatvsCtrl_fixed_names.txt") %>% 
        here::here(Selection_dirs[condition],.) %>% fread() %>% subset(group != "none") %>% pull(Gene)
    # d14_genes <- list.files(Selection_dirs[condition]) %>% str_subset("^Data_ScatterView_TreatvsCtrl.txt") %>% 
    #     here::here(Selection_dirs[condition],.) %>% fread() 
    # subset(d14_genes,group =="bottom") %>% pull(Gene)
    # 
    # All_interesting_genes <- c(All_interesting_genes,Interesting_genes)
    # All_interesting_genes %>% cat(sep = "\n")
    Unnorm_genes <- fread(here::here("Datasets","Processed","Amandine_screen",glue::glue("{condition}_cnv_ntsg.mle.gene_summary.txt")))[,.(Gene,`dmso|beta`,`plx|beta`)]
    essential_genes <- readxl::read_xls(here::here("Datasets","Raw","NIHMS1058492-supplement-Supplementary_Data_2.xls"),skip = 1)[,2] %>% unlist  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6862721/#SD1
    non_essential <- readxl::read_xls(here::here("Datasets","Raw","NIHMS1058492-supplement-Supplementary_Data_1.xls"),skip = 0)[,1] %>% unlist
    Unnorm_genes %>% 
        mutate(essentiality = case_when(
            Gene %in% essential_genes~ "Essential",
            Gene %in% non_essential~ "Non-essential",
            TRUE~"Other"
        )) %>% #View()
        ggplot(aes(x = `dmso|beta`, y= `plx|beta`, colour = essentiality))+
        geom_point(data = . %>% subset(essentiality =="Other"), size = 2)+
        geom_point(data = . %>% subset(essentiality !="Other"), size = 2)+theme_bw()+
        scale_colour_manual(values = c(Essential = scales::muted("red"),
                                       `Non-essential` = scales::muted("blue"),
                                       Other = "grey70"))+
        ggtitle("Essential genes are more essential in both dmso and etop condition at d14")
    ggsave(here::here("Output","Figures","D14_dmso_etop_essentiality_mageckflute.pdf"))
    norm_genes <- list.files(Selection_dirs[condition]) %>% str_subset("squareview_data_fix") %>% here::here(Selection_dirs[condition],.) %>% fread() %>% 
        .[,.(Gene,dmso,plx,Depmap)]%>% 
        mutate(Etop_vs_DMS0 = plx-dmso) %>% 
        arrange(-Etop_vs_DMS0) %>% 
        mutate(Gene = as.factor(Gene),
               Rank = order(Etop_vs_DMS0),
               Significant = case_when(Etop_vs_DMS0<0 & (Gene %in% Interesting_genes)~"Low",
                                       Etop_vs_DMS0>0&  (Gene %in% Interesting_genes)~"High",
                                       TRUE~NA_character_)) 
    norm_genes$Significant %>% table()
    
    hits_d14 <- norm_genes %>% subset(!is.na(Significant)) %>% pull(Gene) 
    
    set.seed(1234)
    KEGG_pathways <- read.csv(here::here("Datasets","Raw","KEGG_genes.csv"))
    norm_genes %>%  
        mutate(Random_index= sample(Rank),
               ID = if_else(!is.na(Significant),as.character(Gene),NA_character_),
               KEGG = case_when(
                   Gene %in% unlist(subset(core_enrichment_KEGG, 
                                           Description == "Chemical carcinogenesis - reactive oxygen species" & comparison == "D14_BP") %>%
                                        pull(Genes))~"Chemical carcinogenesis - reactive oxygen species",
                   Gene %in% unlist(subset(core_enrichment_BP,
                                           Description == "pyrimidine deoxyribonucleotide metabolic process" & comparison == "D10_high_BP") %>%
                                        pull(Genes))~"pyrimidine deoxyribonucleotide metabolic process",
                   
                   
                   Gene %in% unlist(subset(core_enrichment_KEGG,
                                           Description == "Citrate cycle (TCA cycle)" & comparison == "D14_BP") %>%
                                        pull(Genes))~"Citrate cycle (TCA cycle)",
                   # Gene %in% unlist(subset(core_enrichment_KEGG,
                   #                         Description == "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis" & comparison == "D14_BP") %>%
                   #                      pull(Genes))~"Glycosylphosphatidylinositol (GPI)-anchor biosynthesis",
                   TRUE~NA_character_
               ),
               Significant_alpha = if_else(!is.na(Significant)|!is.na(KEGG),T,F)) %>% 
        ggplot(aes(x= Random_index, y= Etop_vs_DMS0,label = ID,colour =KEGG , alpha = Significant_alpha))+
        geom_point(size = 3)+theme_bw()+
        # annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=-0.1645, alpha=0.05, fill=scales::muted("red")) +
        # annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.1645 , ymax=Inf, alpha=0.05, fill=scales::muted("blue")) +
        theme( panel.grid.major = element_blank(),legend.position = "none",text = element_text(size = 25),
               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
        # geom_point(data = . %>% subset(is.na(Significant)), alpha = 0.1)+
        ggrepel::geom_text_repel(data = . %>% subset(!is.na(Significant) & abs(Etop_vs_DMS0)>0.24), aes(size = abs(Etop_vs_DMS0)),max.overlaps = 7)+
        lims(y = c(-0.4,0.4))+
        geom_abline(slope  = 0,intercept = 0.1645, alpha = 0.5,linetype="dashed")+
        geom_abline(slope  = 0,intercept = -0.1645, alpha = 0.5,linetype="dashed")+
        scale_colour_manual(values = 
                                c("Chemical carcinogenesis - reactive oxygen species" = "#A95049" ,
                                  "Citrate cycle (TCA cycle)" ="#363471",
                                  "pyrimidine deoxyribonucleotide metabolic process" = "#1EBCD8"),
        )+
        scale_size_continuous(range = c(7, 15)) +
        
        ggtitle("D14 Difference between Etop and DMSO Annotating Significant KEGG terms")
    
    
    ggsave(here::here("Output","Figures","D14_rank_plot_annotated_illustrator.pdf"), height = 10, width = 20)
    norm_genes %>%  
        mutate(Random_index= sample(Rank),
               ID = if_else(!is.na(Significant),as.character(Gene),NA_character_),
               KEGG = case_when(
                   Gene %in% (go_KEGG %>% subset(id == "hsa05208") %>% pull(ID))~
                       "Chemical carcinogenesis - reactive oxygen species",
                   Gene %in% (go_KEGG %>% subset(id == "hsa00020") %>% pull(ID))~
                       "Citrate_cycle_(TCA_cycle)",
                   # Gene %in% unlist(subset(core_enrichment_KEGG,
                   #                         Description == "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis" & comparison == "D14_BP") %>%
                   #                      pull(Genes))~"Glycosylphosphatidylinositol (GPI)-anchor biosynthesis",
                   TRUE~NA_character_
               ),
               Significant_alpha = if_else(!is.na(Significant)|!is.na(KEGG),T,F)) %>% 
        ggplot(aes(x= Random_index, y= Etop_vs_DMS0,label = ID,colour =KEGG , alpha = Significant_alpha))+
        geom_point()+theme_bw()+
        annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=-0.1645, alpha=0.05, fill=scales::muted("red")) +
        annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.1645 , ymax=Inf, alpha=0.05, fill=scales::muted("blue")) +
        theme( panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
        # geom_point(data = . %>% subset(is.na(Significant)), alpha = 0.1)+
        ggrepel::geom_text_repel(aes(size = abs(Etop_vs_DMS0)))+
        lims(y = c(-0.4,0.4))+
        geom_abline(slope  = 0,intercept = 0.1645, alpha = 0.5,linetype="dashed")+
        geom_abline(slope  = 0,intercept = -0.1645, alpha = 0.5,linetype="dashed")+
        
        ggtitle("D14 Difference between Etop and DMSO Annotating Significant KEGG terms")
    ggsave(here::here("Output","Figures","D14_rank_plot_annotated_all_genes.pdf"))
    norm_genes$Significant %>% table()
    top5_hits_D5 <- norm_genes %>% mutate(abs_Etop_vs_DMS0 = abs(Etop_vs_DMS0)) %>%
        arrange(-abs_Etop_vs_DMS0) %>%  head(17) %>% 
        pull(Etop_vs_DMS0, Gene) 
    
    Hits <- c(top5_hits_D5,top5_hits)
    Hits %>% names %>% cat(sep = "\n")
    Essential_genes_correction <- rbind(Unnorm_genes %>% set_names(c("Gene","DMSO_D14","Etop_D14")) %>% 
                                            mutate(Type = "Unnormalised"),
                                        norm_genes %>% dplyr::select(Gene,dmso,plx) %>% set_names(c("Gene","DMSO_D14","Etop_D14")) %>% 
                                            mutate(Type = "Cell_cycle_Norm")) %>% subset(Gene %in% mageckflute_essential)
    ggplot(Essential_genes_correction,aes(x = DMSO_D14, y = Etop_D14, colour = Type))+
        geom_point(alpha = 1)+
        geom_smooth(method="lm", alpha = 0.5)+ theme_bw()+ theme(panel.grid.major = element_blank(), 
                                                                 panel.grid.minor = element_blank(),
                                                                 panel.background = element_blank())+
        geom_abline(intercept = 0, slope = 1)+
        lims(x= c(-1.8,0),
             y = c(-1.8,0))+
        ggtitle("Cell Cycle Normalisation For Day 14 Screen",
                subtitle = "After normalisation, Essential genes are overall equally depleted in Etop as in DMSO against the library")
    ggsave(here::here("Output","Figures","D14_cell_cycle_norm_Essential_Genes.pdf"))
    list_of_comp[[condition]] <- norm_genes
    ggplot(norm_genes,aes(x = Etop_vs_DMS0, y = Rank, label = Gene, colour = Significant))+
        geom_point(data = . %>% subset(is.na(Significant)))+
        geom_point(data = . %>% subset(!is.na(Significant)))+
        ggrepel::geom_label_repel(data = . %>% subset(!is.na(Significant) & !between(Etop_vs_DMS0,-0.26,0.25 )),max.overlaps = 14,
                                  size = 3)+
        scale_colour_manual(values = c(High = "red", Low = "blue"))+
        theme_bw()+
        ggtitle(condition,
                subtitle = "High is positive, and Low is negative selection for Etop_vs_DMSO")
    ggsave(filename = here::here(here::here("Output","Figures"), glue::glue("CNV_corrected_{condition}_Etop_vs_DMSO.pdf")))
    gene_list <- norm_genes %>% 
        left_join(Human_hsa, by = c("Gene" ="ID")) %>% 
        arrange(-Etop_vs_DMS0)
    Go_list <- gene_list %>% 
        pull(Etop_vs_DMS0,Gene)
    D14_BP <- gseGO(geneList     =Go_list,
                    OrgDb        = org.Hs.eg.db,
                    ont          = "BP",
                    keyType = "SYMBOL",
                    #nPerm        = 1000,
                    minGSSize    = 10,
                    maxGSSize    = 500,
                    pvalueCutoff = 1,
                    #qvalueCutoff = 1,
                    verbose      = FALSE)
    ego3 <- gseGO(geneList     =Go_list,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "BP",
                  keyType = "SYMBOL",
                  #nPerm        = 1000,
                  minGSSize    = 10,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)
    if((ego3@result %>% nrow)>0){
        enrichplot::ridgeplot(ego3 %>%  clusterProfiler::simplify(), showCategory = 68)+
            ggtitle(glue::glue("{condition} Etop_vs_DMSO  Gene BP-GSEA"))
        ggsave(here::here(output_folder, glue::glue("{condition} Etop_vs_DMSO  BP-GSEA.pdf")), height = 20, width  = 15)
    }
    
    ego3 <- gseGO(geneList     =Go_list,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "MF",
                  keyType = "SYMBOL",
                  #nPerm        = 1000,
                  minGSSize    = 10,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)
    
    if((ego3@result %>% nrow)>0){
        enrichplot::ridgeplot(ego3 %>% clusterProfiler::simplify(), showCategory = 68)+
            ggtitle(glue::glue("{condition} Etop_vs_DMSO  Gene  MF-GSEA"))
        ggsave(here::here(output_folder, glue::glue("{condition} Etop_vs_DMSO  MF-GSEA.pdf")), height = 20, width  = 15)
    }
    
    KEGG_list <- gene_list%>% dplyr::select(Etop_vs_DMS0,KEGG) %>% 
        na.omit()%>%  pull(Etop_vs_DMS0,KEGG)
    
    D5_KEGG <- gseKEGG(geneList     = KEGG_list,
                       organism     = 'hsa',
                       minGSSize    = 10,
                       pvalueCutoff = 0.05,
                       verbose      = FALSE)
    terms_to_reduce <- D5_KEGG@result %>% subset(`p.adjust`<0.05 ) %>% 
        dplyr::select(core_enrichment,Description) %>% subset(!duplicated(Description)) %>% 
        pull(core_enrichment,Description)  %>% purrr::map(str_split,"/") %>% flatten()
    Interesting_genes_chemical <- c(Enrichment_KEGG_D14 %>% subset(Category == "Oxphos_entrez") %>% pull(ID),
                                    "ASS1","SLC38A2")
    norm_genes %>% ggplot(aes(x = dmso, y = plx,label = Gene))+
        geom_point(data = . %>% subset(!(Gene %in% Interesting_genes_chemical)), colour = "grey70")+
        
        geom_point(data = . %>% subset((Gene %in% Interesting_genes_chemical)), colour = "red")+
        ggrepel::geom_label_repel(data = . %>% subset((Gene %in% Interesting_genes_chemical)),max.overlaps = 100)+
        ggtitle("D14 Screen")
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
            index <- index[JaccIndex == row_max][,top_similar:= fifelse(str_detect(Pathway1,"Inosit"),Pathway2,Pathway1)]
            
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
    reduced_terms <- Jaccard_index_list(terms_to_reduce,max_jacc = 0.5)
    D5_KEGG@result <- D5_KEGG@result %>% subset(Description %in% (reduced_terms %>% names()) )
    
    enrichplot::dotplot(D5_KEGG)+ggtitle("Enrichent D14 Etop vs D14 DMSO")
    if((D5_KEGG@result %>% nrow)>0){
        
        enrichplot::ridgeplot(D5_KEGG, showCategory = 100)+
            ggtitle(glue::glue("{condition} Etop_vs_DMSO  Gene  KEGG "))
        ggsave(here::here(output_folder, glue::glue("{condition} Etop_vs_DMSO  KEGG.pdf")), height = 20, width  = 15)}
    
    
    mkk2 <- gseMKEGG(geneList = KEGG_list,
                     organism = 'hsa',
                     minGSSize = 10,
                     pvalueCutoff = 1)
    
    if((mkk2@result %>% nrow)>0){
        
        enrichplot::ridgeplot(mkk2, showCategory = 68)+
            ggtitle(glue::glue("{condition} Etop_vs_DMSO  Gene MKEGG "))
        mkk2@result %>% 
            ggplot(aes(x = NES, y = -log10(pvalue),label = Description, colour = NES))+
            geom_point()+
            ggrepel::geom_label_repel(data = . %>% subset(p.adjust<0.05))+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                           panel.background = element_blank())+
            scale_colour_gradient2()+ ggtitle("MKEGG enrichment D14_etop vs D14_DMSO",
                                              subtitle = "So negative pathways are more essential in Etop, this is opposite to IMPDH2 which is significant but less essential in Etop")
        glimpse()
        dotplot(mkk2, showCategory = 68)+
            ggtitle(glue::glue("{condition} Etop_vs_DMSO  Gene MKEGG "))
        ggsave(here::here(output_folder, glue::glue("{condition} Etop_vs_DMSO  MKEGG.pdf")), height = 20, width  = 15)
        
        Comparison <- Unnorm_genes[norm_genes, on = "Gene"][,Essentiality := fcase(
            between(Depmap,-10,-0.5,),"Very Essential",
            between(Depmap,-0.5,-0.25,),"Essential",
            between(Depmap,-0.25,0.2),"Not Essential",
            between(Depmap,0.2,10),"Growth") ]
        Comparison %>% ggplot(aes(x = `plx|beta`- `dmso|beta`, y= plx-dmso, colour = Depmap))+
            geom_point(alpha = 0.5)+geom_abline(yintercept=0, slope=1)+
            scale_colour_gradient2()+facet_wrap("Essentiality")+
            labs(x = "Etop - DMSO Unnormalised",
                 y = "Etop - DMSO normalised")+
            ggtitle(condition)
        Comparison[,lfc := plx-dmso]
        ggsave(here::here("Output",glue::glue("{condition} normalisation MAGeCKflute cell cycle.pdf")))
    }}
columns_interest <- c("Gene", "library" ,'10d', "10d-etop")
counts <- fread(here::here("Datasets","Processed","Amandine_screen", "Amandine.count.txt"), select = columns_interest)
median_normalised <- counts[,-1] %>% as.matrix() %>% median_normalization() %>% as.data.table()
median_normalised[,Gene := counts$Gene]
counts <- median_normalised[, lapply(.SD, log2), by = Gene][, lapply(.SD, median), by = Gene]
essentiality <- norm_genes[,Essentiality := fcase(
    between(Depmap,-10,-0.7,),"Very Essential",
    between(Depmap,-0.7,-0.25,),"Partially Essential",
    between(Depmap,-0.25,0.2),"Not Essential",
    between(Depmap,0.2,10),"Growth") ][,.(Gene,Essentiality,Depmap)]
counts <- essentiality[counts, on = "Gene"]
counts[,`:=`(DMSO_d10_vs_library = `10d`-library,
             Etop_d10_vs_library = `10d-etop`-library)]
ggplot(counts,aes(x = DMSO_d10_vs_library, y = Etop_d10_vs_library, colour = Depmap ))+
    geom_point()+facet_wrap("Essentiality")+
    geom_abline(yintercept=0, slope=1)+
    theme_bw()+
    lims(x= c(-2,0.75), y = c(-2,0.75))+
    scale_colour_gradient2(mid = "grey90")+
    ggtitle("Very Essential Genes are less depleted in Etop_10_vs_d0, thant DMSO_d10_vs_d0")
ggsave(filename = here::here(here::here("Output","Figures"), glue::glue("Cell_cycle_correction_D10.pdf")))

columns_interest <- c("Gene", "library" ,'14d', "14d-etop")
counts <- fread(here::here("Datasets","Processed","Amandine_screen", "Amandine.count.txt"), select = columns_interest)
median_normalised <- counts[,-1] %>% as.matrix() %>% median_normalization() %>% as.data.table()
median_normalised[,Gene := counts$Gene]
counts <- median_normalised[, lapply(.SD, log2), by = Gene][, lapply(.SD, median), by = Gene]
essentiality <- norm_genes[,Essentiality := fcase(
    between(Depmap,-10,-0.7,),"Very Essential",
    between(Depmap,-0.7,-0.25,),"Partially Essential",
    between(Depmap,-0.25,0.2),"Not Essential",
    between(Depmap,0.2,10),"Growth") ][,.(Gene,Essentiality,Depmap)]
counts <- essentiality[counts, on = "Gene"]
counts[,`:=`(DMSO_d14_vs_library = `14d`-library,
             Etop_d14_vs_library = `14d-etop`-library)]
ggplot(counts,aes(x = DMSO_d14_vs_library, y = Etop_d14_vs_library, colour = Depmap ))+
    geom_point()+facet_wrap("Essentiality")+
    geom_abline(yintercept=0, slope=1)+
    theme_bw()+
    lims(x= c(-2,0.75), y = c(-2,0.75))+
    scale_colour_gradient2(mid = "grey90")+
    ggtitle("Very Essential Genes are less depleted in Etop_14_vs_d0, thant DMSO_d14_vs_d0")
ggsave(filename = here::here(here::here("Output","Figures"), glue::glue("Cell_cycle_correction_D14.pdf")))
Combined <- list_of_comp %>% purrr::imap(.x = ., ~set_names(.x,paste0(colnames(.x),.y))) %>% 
    reduce(full_join,by = c("GeneD10_D10_etop" =
                                "GeneD14_D14_etop" )) %>% 
    mutate(Essentiality = case_when(
        DepmapD14_D14_etop<(-1)~"Essential_genes",
        DepmapD14_D14_etop>(0)~"Non_essential_genes",
        TRUE~NA_character_))
ggplot(Combined,aes(x = dmsoD10_D10_etop,y = dmsoD14_D14_etop, colour = Essentiality))+
    geom_point(data =. %>%  subset(is.na(Essentiality)),size = 2)+
    geom_point(data =. %>%  subset(!is.na(Essentiality)),size = 2)+
    theme_bw()+theme(#panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_colour_manual(values = c("Essential_genes" = "#B1172B",
                                   "Non_essential_genes" = "#2066AC"))+
    ggtitle("Essential Genes are depleted in both D10 and D14 DMSO more than none Essential Genes")
ggsave(here::here("Output","Figures","Supp_2B_Essentiality_Depletetion.pdf"))

Complexes <- read.delim(here::here("Datasets","Raw","coreComplexes.txt")) %>% 
    dplyr::select(ComplexName,Organism,`subunits.Gene.name.`) %>% subset(Organism == "Human") %>% 
    separate_rows(`subunits.Gene.name.`,sep =';') %>% dplyr::rename(Genes =`subunits.Gene.name.`,
                                                                    Complex = ComplexName) %>% dplyr::select(Genes,Complex)
Complexes_mito <- read.delim(here::here("Datasets","Raw","allComplexes.txt")) %>% dplyr::select(`Approved.symbol`,`Group.name`) %>% 
    dplyr::rename(Genes =`Approved.symbol`,
                  Complex =`Group.name` )
Complexes <- rbind(Complexes_mito,Complexes) %>% distinct(Genes,.keep_all = T)
Complexes <- Complexes_mito

chemical_carc_complex <- core_enrichment_KEGG %>% subset(Description == "Chemical carcinogenesis - reactive oxygen species" & Significant == T) %>% 
    dplyr::select(comparison, Genes) %>% unnest(Genes)
chemical_carc_complex <- chemical_carc_complex %>% left_join(Complexes)
#chemical_carc_complex <- chemical_carc_complex %>% mutate(Complex = if_else(is.na(Complex),"No Complex",Complex))
Genes_in_complexes <- chemical_carc_complex %>% na.omit() %>% pull(Genes)
complex_graph_edges <- data.table()
for(i in Genes_in_complexes){
    tmp_complex <-   chemical_carc_complex %>% 
        subset(Genes == i) %>% pull(Complex)
    complex_graph_edges <- rbind(complex_graph_edges,
                                 data.table(from = i,
                                            to =chemical_carc_complex %>% 
                                                subset(Complex == tmp_complex) %>% pull(Genes) ))
    
}
complex_graph_edges <- complex_graph_edges[from>to] 
for(i in Genes_in_complexes){
    
    complex_graph_edges <- rbind(complex_graph_edges,
                                 data.table(from = i,
                                            to =chemical_carc_complex %>% 
                                                subset(Genes == i) %>% pull(comparison) ))
    
}
complex_graph_edges <- complex_graph_edges%>% 
    mutate(colour = if_else(str_detect(from,"vs")|
                                str_detect(to,"vs"),"grey40","grey80"))
library("RColorBrewer")
nodes_df <- data.table(id = unique(c(complex_graph_edges$from,
                                     complex_graph_edges$to)) ) %>% 
    left_join(unique(chemical_carc_complex %>% dplyr::select(Genes,Complex)),by = c("id" = "Genes")) %>% 
    mutate(Complex = if_else(is.na(Complex),"Comparison",Complex),
           value = if_else(Complex == "Comparison",15,7))
colour_pallete <- data.frame(color.background = c(brewer.pal(n = 4, name = 'Set3'),
                                                  brewer.pal(n = length(unique(nodes_df$Complex))-4, name = 'Set2')
))%>% mutate(Complex = unique(nodes_df$Complex))

nodes_df <- left_join(nodes_df,colour_pallete)
nodes_df <- left_join(nodes_df,nodes_df %>% dplyr::select(Complex,id) %>% 
                          group_by(Complex) %>% 
                          top_n(50) %>%
                          dplyr::mutate(label = id)) %>% 
    ungroup %>% mutate(label = if_else(Complex == "No Complex",id,label))
library(igraph)



complex_graph_edges <- distinct(complex_graph_edges)
net <- graph_from_data_frame(d=complex_graph_edges, vertices=nodes_df, directed=F) 
E(net)$colour <-as.numeric(if_else(str_detect(complex_graph_edges$from,"_")|
                                       str_detect(complex_graph_edges$to,"_"),0.8,0.5))

set.seed(1234)
ggraph(net, layout="dh") +
    geom_edge_fan(aes(alpha = colour),colour = "grey20", width=0.8) + 
    
    scale_edge_alpha_manual(values = c(0.03,0.08))+
    geom_node_point(#shape = 21,
        # fill=V(net)$color.background, 
        size=V(net)$value,
        
        
        colour =  V(net)$color.background,
        #stroke = 1
    ) +
    geom_node_text(label = V(net)$label, size=5, 
                   color="gray30", repel=T) +
    theme_void()
##### average enrichment per omic###
chemical_carc_complex <- core_enrichment_KEGG %>% subset(Description == "Chemical carcinogenesis - reactive oxygen species" & Significant == T) %>% 
    dplyr::select(comparison, Genes) %>% unnest(Genes)
chemical_carc_complex <- chemical_carc_complex %>% left_join(Complexes_mito) #%>% na.omit()

Proteomics <- Volcano_DFs$dmso_vs_x24h_diff %>% 
    mutate(EtopvsDMSO = scale(rank(-log2_FC)),
           comparison  = "dmso_vs_x24h_diff") %>% 
    dplyr::select(EtopvsDMSO,ID,comparison) %>% 
    dplyr::rename( Genes= ID) 
Screen = norm_genes   %>%
    mutate(EtopvsDMSO = -scale(rank(Etop_vs_DMS0)),
           comparison  = "D14_BP") %>%
    dplyr::select(EtopvsDMSO,Gene,comparison)%>% 
    dplyr::rename( Genes= Gene ) 


Mito_complex_genes <- rbind(inner_join(chemical_carc_complex,Proteomics),
                            inner_join(chemical_carc_complex,Screen))%>% 
    mutate(Complex = str_remove_all(Complex,"oxidoreductase [:graph:]*"))
openxlsx::write.xlsx(Mito_complex_genes, here::here("Datasets","Processed","Mito_complex_genes_enrichment.xlsx"), overwrite = T)
Mito_complex_genes <- rbind(inner_join(chemical_carc_complex,Proteomics),
                            inner_join(chemical_carc_complex,Screen))%>% 
    mutate(Complex = str_remove_all(Complex,"oxidoreductase [:graph:]*"))
openxlsx::write.xlsx(Mito_complex_genes, here::here("Datasets","Processed","Mito_complex_genes_enrichment.xlsx"))

Mito_complex_genes <-     Mito_complex_genes %>%     mutate(tally_add = 1) %>% 
    group_by(comparison, Complex) %>% 
    dplyr::summarise(median_Complex_rank = median(EtopvsDMSO, na.rm = T),
                     Count = sum(tally_add)) %>% 
    left_join(Complexes_mito %>% 
                  mutate(Complex = str_remove_all(Complex,"oxidoreductase [:graph:]*"))%>%
                  dplyr::count(Complex)) %>% 
    mutate(portion_of_complex =Count/n )
rbind(Mito_complex_genes,
      data.table(Type = "Error",
                 portion_of_complex = 0,
                 Complex = "error",
                 median_Complex_rank =  0)) %>% 
    ggplot(aes(x = rank(median_Complex_rank), y = median_Complex_rank, colour = median_Complex_rank,label = Complex
    ))+
    geom_point(aes(size = portion_of_complex))+ theme_bw()+
    ggrepel::geom_text_repel()+
    scale_colour_gradient2(low ="#053061" ,mid = "white", high = "#670A1F")
ggsave(here::here("Output","Figures","ETC_plot.pdf"))

