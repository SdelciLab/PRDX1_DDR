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
high_low_hits = read_tsv(here::here("Datasets","Processed","cell_cycle_hits.tsv"))

conditions <- c("D14_D14_etop")
All_interesting_genes <- high_low_hits$Gene %>% as.character()
list_of_comp <- list()
# for(condition in conditions){
    condition <- conditions[1]
    Interesting_genes <- fread(here::here("Datasets","Processed","Data_ScatterView_TreatvsCtrl_fixed_names_D14.txt")) %>% 
        subset(group != "none") %>% pull(Gene)
    Unnorm_genes <- fread(here::here("Datasets","Processed",glue::glue("{condition}_cnv_ntsg.mle.gene_summary.txt")))[,.(Gene,`dmso|beta`,`plx|beta`)]
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
    ggsave(here::here("Output","D14_dmso_etop_essentiality_mageckflute.pdf"))
    norm_genes <- fread(here::here("Datasets","Processed","D14_D14_etop_squareview_data_fixed_names.txt")) %>% 
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
    
    
    ggsave(here::here("Output","D14_rank_plot_annotated_illustrator.pdf"), height = 10, width = 20)
    norm_genes$Significant %>% table()
    top5_hits_D5 <- norm_genes %>% mutate(abs_Etop_vs_DMS0 = abs(Etop_vs_DMS0)) %>%
        arrange(-abs_Etop_vs_DMS0) %>%  head(17) %>% 
        pull(Etop_vs_DMS0, Gene) 
    
    Hits <- c(top5_hits_D5)
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
    ggsave(here::here("Output","D14_cell_cycle_norm_Essential_Genes.pdf"))

    
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
     
        ggsave(here::here("Output",glue::glue("normalisation MAGeCKflute cell cycle.pdf")))

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

Complexes_mito <- read.delim(here::here("Datasets","Raw","allComplexes.txt")) %>% dplyr::select(`Approved.symbol`,`Group.name`) %>% 
    dplyr::rename(Genes =`Approved.symbol`,
                  Complex =`Group.name` )

Mito_complex_genes <- openxlsx::read.xlsx(here::here("Datasets","Processed","Mito_complex_genes_enrichment.xlsx"))
Mito_complex_genes <-     Mito_complex_genes %>%    na.omit() |> 
    mutate(tally_add = 1) %>% 
    group_by(comparison, Complex) %>% 
    dplyr::summarise(median_Complex_rank = median(as.numeric(EtopvsDMSO), na.rm = T),
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
ggsave(here::here("Output","ETC_plot.pdf"))

