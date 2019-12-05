library(tidyverse)
library(stringr)
library(KEGGREST)
library(AnnotationForge)

#' Title
#'
#' @param f_emapper_anno eggnog-mapper annotation result
#' @param author Who is the creator of this package? like "xxx <xxx@xxx.xx>"
#' @param tax_id The Taxonomy ID that represents your organism. (NCBI has a nice online browser for finding the one you need)
#' @param genus Single string indicating the genus
#' @param species Single string indicating the species
#'
#' @return OrgDb name
#' @export
#'
#' @examples
#' 
makeOrgPackageFromEmapper <- function(f_emapper_anno, 
                                      author, 
                                      tax_id = "0", 
                                      genus = "default", 
                                      species = "default") {
  
  # read emapper result
  emapper <- read_delim(f_emapper_anno,
                        "\t", escape_double = FALSE, trim_ws = TRUE)
  
  # extract gene name from emapper
  gene_info <- emapper %>%
    dplyr::select(GID = query_name, GENENAME = `eggNOG annot`) %>%
    na.omit()
  
  # extract go annotation from emapper
  gos <- emapper %>%
    dplyr::select(query_name, GO_terms) %>%
    na.omit()
  
  gene2go = data.frame(GID = character(),
                       GO = character(),
                       EVIDENCE = character())
  
  df_temp <- list()
  for (row in 1:nrow(gos)) {
    the_gid <- gos[row, "query_name"][[1]]
    the_gos <- str_split(gos[row,"GO_terms"], ",", simplify = FALSE)[[1]]
    
    df_temp[[row]] <- data_frame(GID = rep(the_gid, length(the_gos)),
                                 GO = the_gos,
                                 EVIDENCE = rep("IEA", length(the_gos)))
  }
  gene2go <- bind_rows(df_temp)
  
  # extract kegg pathway annotation from emapper
  gene2ko <- emapper %>%
    dplyr::select(GID = query_name, Ko = KEGG_KOs) %>%
    na.omit()
  
  load(file = "kegg_info.RData")
  gene2pathway <- gene2ko %>% left_join(ko2pathway, by = "Ko") %>% 
    dplyr::select(GID, Pathway) %>%
    na.omit()
  
  # extract COG annotation from emapper -------------------------------------
  cog_info <- read_delim("cog_funclass.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
  
  cogs <- emapper %>%
    dplyr::select(query_name, COG = "COG cat") %>%
    na.omit()
  
  gene2cog = data.frame(GID = character(),
                        COG = character())
  
  df_temp <- list()
  for (row in 1:nrow(cogs)) {
    the_gid <- cogs[row, "query_name"][[1]]
    the_cogs <- str_trim(str_split(cogs[row,"COG"], ",", simplify = FALSE)[[1]])
    
    df_temp[[row]] <- data_frame(GID = rep(the_gid, length(the_cogs)),
                                 COG = the_cogs)
  }
  gene2cog <- bind_rows(df_temp)
  
  gene2cog <- gene2cog %>% left_join(cog_info, by = "COG")
  
  save(gene_info, gene2go, gene2pathway, gene2cog, file = "gene_annotation.RData")

  
  # make OrgDb
  makeOrgPackage(gene_info=gene_info,
                 go=gene2go,
                 ko=gene2ko,
                 pathway=gene2pathway,
                 cog=gene2cog,
                 maintainer=author,
                 author=author,
                 outputDir=".",
                 tax_id=tax_id,
                 genus=genus,
                 species=species,
                 goTable="go",
                 version="1.0")
  
  my_orgdb <- str_c("org.", str_to_upper(str_sub(genus, 1, 1)) , species, ".eg.db", sep = "")
  return(my_orgdb)
}

###################################################

my_orgdb <- makeOrgPackageFromEmapper("genome_eggnog.emapper2.annotations", 
                                      "test <test@genek.tv>", 
                                      tax_id = "0000", 
                                      genus = "G", 
                                      species = "E")
my_orgdb <- makeOrgPackageFromEmapper("effectome_eggnog.emapper.annotations", 
                                      "test <test@genek.tv>", 
                                      tax_id = "0000", 
                                      genus = "E", 
                                      species = "F")
my_orgdb <- makeOrgPackageFromEmapper("contraction_extension.emapper.annotations", 
                                      "test <test@genek.tv>", 
                                      tax_id = "0000", 
                                      genus = "C", 
                                      species = "E")

dir.create("R_library_contraction_extension")
install.packages(my_orgdb, repos = NULL, lib = "R_library_contraction_extension")


dir.create("R_library_genome")
install.packages(my_orgdb, repos = NULL, lib = "R_library_genome")

install.packages("org.GE.eg.db", repos = NULL, lib = "R_library_genome")

library(my_orgdb, character.only = TRUE, lib.loc = "R_library_genome")

library("org.CE.eg.db", character.only = TRUE, lib.loc = "R_library_contraction_extension")

##################功能注释###########################

# load("gene_annotation.RData")
library(clusterProfiler)
# groupGO(gene, OrgDb, keyType = "ENTREZID", ont = "CC", level = 2,
  readable = FALSE)

list <- as.character(read.table("effectome_Biotoph_specialst_gene_id.txt", quote="\"", comment.char="")$V1)

group_go_bp <- groupGO(gene     = list,
                   OrgDb    = my_orgdb,
                   keyType  = "GID",
                   ont      = "BP",
                   level    = 3,
                   readable = FALSE)
  
 group_go_bp <- as.data.frame(group_go_bp)
 group_go_bp$GO_Class <- "Biological Process"
  
 group_go_cc <- groupGO(gene     = list,
                   OrgDb    = my_orgdb,
                   keyType  = "GID",
                   ont      = "CC",
                   level    = 3,
                   readable = FALSE)
  
 group_go_cc <- as.data.frame(group_go_cc)
 group_go_cc$GO_Class <- "Cellular Component"
  
  group_go_mf <- groupGO(gene     = list,
                   OrgDb    = my_orgdb,
                   keyType  = "GID",
                   ont      = "MF",
                   level    = 3,
                   readable = FALSE)
 group_go_mf <- as.data.frame(group_go_mf)
 group_go_mf$GO_Class <- "Molecular Function"
  
 group_go_all <- rbind(group_go_bp, group_go_cc, group_go_mf)
  
write.csv(group_go_all,"effectome_Biotoph_specialst_group_go_all.csv")





####################### encrichGO ################################
list <- as.character(read.table("Gc_contraction_id.txt", quote="\"", comment.char="")$V1)
# load("gene_annotation.RData")
library(clusterProfiler)

enrich_go_bp <- enrichGO(gene = list,
                OrgDb = my_orgdb ,
                ont = "BP",
				keyType  = "GID",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01) 
				
enrich_go_bp <- as.data.frame(enrich_go_bp)
enrich_go_bp$GO_Class <- "Biological Process"

enrich_go_cc <- enrichGO(gene = list,
                OrgDb = my_orgdb,
                ont = "CC",
				keyType  = "GID",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01) 
				
enrich_go_cc <- as.data.frame(enrich_go_cc)
enrich_go_cc$GO_Class <- "Cellular Component"


enrich_go_mf <- enrichGO(gene = list,
                OrgDb = my_orgdb ,
                ont = "MF",
				keyType  = "GID",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01) 
				
enrich_go_mf <- as.data.frame(enrich_go_mf)
enrich_go_mf$GO_Class <- "Molecular Function"

enrich_go_all <- rbind(enrich_go_bp, enrich_go_cc, enrich_go_mf)
write.csv(enrich_go_all,"Gc_contraction_enrich_go_all.csv")

######################作图################################

p<-ggplot(go_all) + 
  geom_bar(aes(x = Description, 
                     y = Count,
                     fill = GO_Class),
                 stat="identity",width=0.8,position="dodge") + facet_wrap(~GO_Class, scales = "free_x",) + 
  labs(title = "GO function classification", y = "Number of genes") +
  theme_set(theme_bw())+
  theme(axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=10),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=15),
        axis.text.x = element_text(vjust=1,hjust=1,angle=45,size=6))
ggsave("Fo_go.pdf",width = 10, height = 7)


qvalue<- go_all$qvalue
p<-ggplot(go_all,aes(GeneRatio,Description)) + 
       geom_point(aes(size = Count,color = -1*log10(qvalue))) + 
	   scale_colour_gradient(low="green",high="red")+
  facet_wrap(~GO_Class, scales = "free_x",) + 
  labs(color = -1*log10(qvalue),size = "Gene number",title = "GO function classification", y = "Number of genes") +
  theme_set(theme_bw())+
  theme(axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=10),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=15),
        axis.text.x = element_text(vjust=1,hjust=1,angle=45,size=6))
ggsave("Fo_go.pdf",width = 10, height = 7)


############ COG ############################################
emapper <- read_delim("genome_eggnog.emapper2.annotations",
                        "\t", escape_double = FALSE, trim_ws = TRUE)
cog_info <- read_delim("cog_funclass.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
  
  cogs <- emapper %>%
    dplyr::select(query_name, COG = "COG cat") %>%
    na.omit()
  
  gene2cog = data.frame(GID = character(),
                        COG = character())
  
  df_temp <- list()
  for (row in 1:nrow(cogs)) {
    the_gid <- cogs[row, "query_name"][[1]]
    the_cogs <- str_trim(str_split(cogs[row,"COG"], ",", simplify = FALSE)[[1]])
    
    df_temp[[row]] <- data_frame(GID = rep(the_gid, length(the_cogs)),
                                 COG = the_cogs)
  }
  gene2cog <- bind_rows(df_temp)
  
  gene2cog <- gene2cog %>% left_join(cog_info, by = "COG")


gene2cog$COG_Name = paste("(", gene2cog$COG, ")", gene2cog$COG_Name, sep = " ")

write.csv(gene2cog,"genome_gene2cog.csv")

p <- ggplot(data = gene2cog) + 
  geom_bar(aes(x = COG, 
               fill = COG_Name)) +
  labs(title = "COG/KOG Function Classification ", 
       x = "",
       y = "Number of genes") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.key.size=unit(1,"line"),
        legend.text = element_text(size = 7.5)) +
  guides(fill=guide_legend(ncol=1))
ggsave(paste(argv$out_dir, "cog.pdf", sep = "/"), p, width = 16, height = 7)
write.csv(gene2cog, paste(argv$out_dir, "genome_gene2cog.csv"))


library(ggplot2)
############# ggplot cog geom_bar ####################

ggplot(data=cog_go_core,mapping=aes(x=Class,y=Number,fill=COG_Name))+
  geom_bar(stat="identity",width=0.8,position="dodge")+
  xlab("")+
  ylab("Number of genes")+
  theme_set(theme_bw())+
  theme(axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(family="Times New Roman",face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(family="Times New Roman",face = "bold",colour="black",size=15),
        axis.text.x = element_text(vjust=1,hjust=1,angle=0,size=11))
ggsave("Fo_COG.pdf",width = 10, height = 7)

##################### KEGG ############################ 

library(purrr)
library(tidyverse)
library(clusterProfiler)
				
pathway2gene <- AnnotationDbi::select(org.GE.eg.db, 
                                        keys = keys(org.GE.eg.db), 
                                        columns = c("Pathway","Ko")) %>%
    na.omit() %>%
    dplyr::select(Pathway, GID)

		
load("kegg_info.RData")

ekp <- enricher(list, 
                  TERM2GENE = pathway2gene, 
                  TERM2NAME = pathway2name, 
                  pvalueCutoff = 1, 
                  qvalueCutoff = 1,
                  pAdjustMethod = "BH",
                  minGSSize = 1)

write.csv(ekp,"effectome_Biotoph_specialst_kegg.csv")
				  
ekp_results <- as.data.frame(ekp)
  
  barplot(ekp, showCategory=20,color="pvalue",
          font.size=10)
  dotplot(ekp)
  
  emapplot(ekp)