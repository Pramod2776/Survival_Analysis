## required libraries
require(survival)
require(survminer)
require(ggplot2)
require("biomaRt")
library("broom")
library("dplyr")

## Data Processing
# TCGA-SKCM gene expression and clinical datasets
expr <- data.table::fread("./skcm_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt",
                   header = T, check.names = F)%>%
  as.data.frame()
rownames(expr) <- make.unique(expr$Hugo_Symbol)
expr <- expr[,3:dim(expr)[2]]

clinical <- read.delim("./skcm_tcga_pan_can_atlas_2018/data_clinical_patient.txt",
                       header = T, check.names = F, skip = 4)
clinical <- clinical[, c(1, 5,6,30, 31)]

samples <- read.delim("./skcm_tcga_pan_can_atlas_2018/data_clinical_sample.txt", skip = 4,
                      header = T, stringsAsFactors = F)
samples <- samples[,1:2]

tmp1 <- clinical[match(samples$PATIENT_ID, clinical$PATIENT_ID), ]
metadat <- cbind(samples, tmp1)
table(samples$PATIENT_ID == tmp1$PATIENT_ID)
metadat <- metadat[,2:7]

metadat$status <- rep(0, dim(metadat)[1])
metadat$status[metadat$OS_STATUS == "1:DECEASED"] <- 1

# Compile expression of selected genes
genes = readxl::read_xlsx("./Composite_File.xlsx", sheet = "Sel_Genes")$Genes
genes_unique = stringr::str_to_title(unique(genes))
homology_table = readxl::read_xlsx("./HGNC_AllianceHomology.xlsx")%>%
  as.data.frame()
homology_table_subset = homology_table %>%
  dplyr::select(Marker_Symbol, HGNC_ID)%>%
  dplyr::filter(Marker_Symbol %in% genes_unique)


homology_table_subset_columns = homology_table_subset %>%
  mutate("HGNC_ID" = strsplit(as.character(homology_table_subset$HGNC_ID), "\\|")) %>% 
  tidyr::unnest("HGNC_ID")%>%
  as.data.frame()

## check any genes missing
setdiff_genes = setdiff(genes_unique, homology_table_subset_columns$Marker_Symbol)

## extract attributes from biomaRT
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
get_genesymbol = getBM(attributes = c( "hgnc_symbol", "hgnc_id"), filters = "hgnc_id", 
      values = homology_table_subset_columns$HGNC_ID, mart = mart)  

colnames(homology_table_subset_columns) = c("Marker_Symbol", "hgnc_id")
homology_table_subset_survival = homology_table_subset_columns %>%
  dplyr::left_join(get_genesymbol, by = "hgnc_id")

genes_expr <- homology_table_subset_survival$hgnc_symbol
expr.sel_surv = expr %>%
  dplyr::filter(rownames(.) %in% genes_expr)

expr.sel_surv <- t(expr.sel_surv)

tmp <- expr.sel_surv[match(metadat$SAMPLE_ID, rownames(expr.sel_surv)), ]
table(rownames(tmp) == metadat$SAMPLE_ID)

metadat <- cbind(metadat, tmp)

# Survival analysis and plots
genes <- colnames(metadat)[8:ncol(metadat)]
res.cut <- surv_cutpoint(metadat, time = "OS_MONTHS", event = "status", variables = genes)
res.cat <- surv_categorize(res.cut)
head(res.cat, n=20)

formulae <- list()
formulae <- lapply(genes, function(x) as.formula(paste0("Surv(`OS_MONTHS`, status) ~ ", x)))
fits <- surv_fit(formulae, data = res.cat)

g = ggsurvplot_list(fit = fits, 
               data = res.cat,
               ####### Format Title #######
               title = "",
               subtitle = "",
               font.title = c(16, "black"),
               ggtheme = theme_bw() + theme(plot.title = element_text(hjust = 0.5))+ # theme_classic will give a white background with no lines on the plot
                 theme(plot.subtitle = element_text(hjust = 0.5, size = 16))+
                 theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), 
                       panel.grid.minor.x = element_blank(), 
                       panel.grid.minor.y = element_blank()),
               ####### Format Axes #######
               xlab="Months elapsed", # changes xlabel,
               ylab = "Survival Probability",
               font.x=c(18), # changes x axis labels
               font.y=c(18), # changes y axis labels
               font.xtickslab=c(14,"plain"), # changes the tick label on x axis
               font.ytickslab=c(14,"plain"),
               ####### Format Curve Lines #######
               palette = c("#bd0026","#1d91c0","grey50"),
               ####### Censor Details ########
               censor = T, # logical value. If TRUE, censors will be drawn,
               censor.shape="|",
               censor.size = 5,
               ####### Confidence Intervals ########
               conf.int = F, # To Remove conf intervals use "FALSE"
               conf.int.fill = "purple", # fill color to be used for confidence interval
               surv.median.line = "none", # allowed values include one of c("none", "hv", "h", "v"). v: vertical, h:horizontal
                ######## Risk Table #######
                risk.table = T, # Adds Risk Table
                risk.table.height = 0.2,
                pval = T,
                pval.size = 5,
                pval.coord = c(3,1)) # Adjusts the height of the risk table (default is 0.25)

pdf("survplot.pdf")
print(g)
dev.off()
dev.set(dev.next())
while (!is.null(dev.list()))  dev.off()


