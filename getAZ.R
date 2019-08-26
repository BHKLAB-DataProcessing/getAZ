
options(stringsAsFactors=FALSE)

myDirPrefix <- "/pfs"

geneMap <- read.csv(file.path(myDirPrefix, "downAnnotations/annot_ensembl_all_genes.csv"))


profiles <- get(load("/pfs/AZrecomp/profiles.rds"))


# profiles <- cbind("auc_recomputed"=recomputed$AUC/100, "ic50_recomputed"=recomputed$IC50, )
# rownames(profiles) <- rownames(Az_raw_sensitivity)
library(Biobase)
load("/pfs/gdscU219normalized/GDSC1000.RData")


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


pp <- pData(GDSC1000@molecularProfiles$rna)
samples <- rownames(pp)[which(pp$cellid %in% curationCell$unique.cellid)]
rna <- GDSC1000@molecularProfiles$rna[,samples]
xx <- which(is.na(curationCell$unique.tissueid))
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

save(AZ, file="AZ.RData")

