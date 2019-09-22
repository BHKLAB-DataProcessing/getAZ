library(PharmacoGx)
options(stringsAsFactors=FALSE)


matchToIDTable <- function(ids,tbl, column, returnColumn="unique.cellid") {
  sapply(ids, function(x) {
                          myx <- grep(paste0("((///)|^)",Hmisc::escapeRegex(x),"((///)|$)"), tbl[,column])
                          if(length(myx) > 1){
                            stop("Something went wrong in curating ids, we have multiple matches")
                          }
        if(length(myx) == 0){return(NA_character_)}
                          return(tbl[myx, returnColumn])
                        })
}


myDirPrefix <- "pfs"

geneMap <- read.csv(file.path(myDirPrefix, "downAnnotations/annot_ensembl_all_genes.csv"))

curationCell <- readRDS(file.path(myDirPrefix,"normalizedAZData/curationCell.rds"))
curationDrug <- readRDS(file.path(myDirPrefix,"normalizedAZData/curationDrug.rds"))


profiles <- get(load(file.path(myDirPrefix, "azProfiles/profiles.RData")))

sensitivity_info <- readRDS(file.path(myDirPrefix, "normalizedAZData/AZ_sens_info.rds"))
Az_raw_sensitivity <- readRDS(file.path(myDirPrefix, "normalizedAZData/AZ_sens_raw.rds"))

profiles <- profiles[rownames(sensitivity_info),]

# profiles <- cbind("auc_recomputed"=recomputed$AUC/100, "ic50_recomputed"=recomputed$IC50, )
# rownames(profiles) <- rownames(Az_raw_sensitivity)
library(Biobase)
load(file.path(myDirPrefix, "gdscU219normalized/GDSC_U219_ENSG.RData"))


message("Loading RNA Data")

load(file.path(myDirPrefix, "gdscU219normalized/GDSC_U219_ENSG.RData"))

cell.all <- read.csv(file.path(myDirPrefix, "downAnnotations/cell_annotation_all.csv"))


rna.cellid <- matchToIDTable(ids=phenoData(cgp.u219.ensg)$Characteristics.cell.line., tbl=cell.all, column = "CGP_EMTAB3610.cellid", returnColumn="unique.cellid")


annot <- geneMap
rownames(annot) <- annot$gene_id
gdsc.u219.ensg <- cgp.u219.ensg
annotation(gdsc.u219.ensg) <- "rna"
ensemblIds <- sapply(strsplit(rownames(exprs(gdsc.u219.ensg)), "_"), function (x) { return (x[[1]]) }) 
fData(gdsc.u219.ensg) <- data.frame("Probe"=rownames(exprs(gdsc.u219.ensg)), 
                          "EnsemblGeneId"=ensemblIds,
                          "EntrezGeneId"=annot[ensemblIds, "EntrezGene.ID"],
                          "Symbol"=annot[ensemblIds, "gene_name"],
                          "GeneBioType"=annot[ensemblIds, "gene_biotype"],
                          "BEST"=TRUE)
rownames(fData(gdsc.u219.ensg)) <- rownames(exprs(gdsc.u219.ensg))
pData(gdsc.u219.ensg)[,"batchid"] <- NA
pData(gdsc.u219.ensg)[,"cellid"] <- rna.cellid

cells <- curationCell[,"unique.cellid"]

myx <- pData(gdsc.u219.ensg)[,"cellid"] %in% cells

gdsc.u219.ensg <- gdsc.u219.ensg[,myx]



xx <- which(is.na(curationCell$unique.tissueid))

## WHY was the following done?

curationCell$unique.tissueid[xx] <- curationCell$disease_type[xx]

curationCell$tissueid <- curationCell$unique.tissueid
curationTissue <- curationCell[,c("tissueid", "unique.tissueid")]




AZ <- PharmacoSet(name="AZ",
                  molecularProfiles=list("rna"=gdsc.u219.ensg),
                  cell=curationCell,
                  drug=curationDrug,
                  sensitivityInfo=sensitivity_info,
                  sensitivityRaw=Az_raw_sensitivity,
                  sensitivityProfiles=profiles,
                  sensitivityN=NULL,
                  curationCell=curationCell,
                  curationDrug=curationDrug,
                  curationTissue=curationTissue, datasetType="sensitivity")

saveRDS(AZ, file="/pfs/out/AZ.rds")


