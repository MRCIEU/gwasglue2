

# create summarySet 1
data1 <- ieugwasr::associations(variants = "5:74132993-75132993", id = "ukb-d-I9_IHD")
meta1 <- create_metadata(metadata = ieugwasr::gwasinfo("ukb-d-I9_IHD"))
sumset1 <- create_summaryset(data1, metadata=meta1)
getMetadata(sumset1)


# create SummarySet2
data2 <- ieugwasr::associations(variants = "5:74132993-75132993", id = "finn-b-I9_CHD")
meta2 <- create_metadata(metadata = ieugwasr::gwasinfo("finn-b-I9_CHD"))
meta2 

sumset2 <- create_summaryset(data2, metadata=meta2, qc = FALSE)

# There were warnings
sumset2 <- create_summaryset(data2, metadata=meta2, qc = TRUE)
# problematic variants removed
getMetadata(sumset2)
#  There is no same_size information in metadata. add ncontrol to sample_size
sumset2 = addToMetadata(sumset2, sample_size=getMetadata(sumset2)$ncontrol +getMetadata(sumset2)$ncase)
getMetadata(sumset2)


# create dataset
summarysets <- list(sumset1, sumset2)
dataset <-  create_dataset(summarysets, harmonise = TRUE, tolerance = 0.08, action = 1)%>%
  # meta-analysis to create a new summary set
  meta_analysis(.) 


bed_ref <- "data/ld/EUR"
summarysets <- list(sumset1, sumset2)
dt <-  create_dataset(summarysets, harmonise = TRUE, tolerance = 0.08, action = 1) %>% 
    # harmonise dataset against LD matrix
    harmonise_ld(., bfile = bed_ref , plink_bin = "plink")
#################################
# LOCAL FILES
# GWAS summary data
library(susieR)
library(hyprcoloc)
d1 <- dplyr::as_tibble(read.table(system.file("tests", "ieu-a-2_TopHits_sumdata.txt", package="gwasglue2")))
d2 <- dplyr::as_tibble(read.table(system.file("tests", "ieu-a-7_sumdata_ieu-0-7TopHits.txt", package="gwasglue2")))
data <- list(d1,d2)

# get metadata and create metadata objects
m1 <- read.table(system.file("tests", "ieu-a-2_metadata.txt", package="gwasglue2"))
m2 <- read.table(system.file("tests", "ieu-a-7_metadata.txt", package="gwasglue2"))
metadata <- list(create_metadata(m1), create_metadata(m2))

bed_ref <- "data/ld/EUR"
# create dataset using create_summaryset() and create_dataset()
  dataset <- lapply(seq_along(data), function(i){
     # create summarysets
    dt <- create_summaryset(data[[i]], metadata=metadata[[i]])
  }) %>%
    # create dataset
    create_dataset(., harmonise = TRUE, tolerance = 0.08, action = 1) %>% 
    # harmonise dataset against LD matrix
    harmonise_ld(., bfile = bed_ref , plink_bin = "plink")

plot_gwasglue(dataset, type="manhattan", title = "tophits")

ntraits <- getLength(dataset)
dataset_marginalised <- lapply(1:ntraits, function(trait)
  {
    # Takes in 1 SS
    # Outputs 1 DS (with at least 1 SS)
    ds <- susie_rss(R = getLDMatrix(dataset), summaryset = getSummarySet(dataset, trait))
})
dataset_marginalised <- merge_datasets(dataset_marginalised)

res_dataset <- hyprcoloc(dataset = dataset) 
#hyprcoloc is using gwasglue2 DataSet class object as input
print(res_dataset)


res_dataset_marginalised <- hyprcoloc(dataset = dataset_marginalised)
print(res_dataset_marginalised)


library(VariantAnnotation)


vcffile <- "data/vcf/IEU-a-2.vcf.gz"
vcf <- VariantAnnotation::readVcf(vcffile)
class(vcf)
vcf
VariantAnnotation::samples(header(vcf))
SummarizedExperiment::rowRanges(vcf)
x=gwasvcf::vcf_to_granges(vcf)%>% dplyr::as_tibble()


data1_vcf <- gwasvcf::query_gwas(vcf, chrompos=c("1:1097291-1099437"))%>% 
gwasvcf::vcf_to_granges() %>% dplyr::as_tibble()
dim(data1_vcf)
data1_ieu <- ieugwasr::associations(variants = "1:1097291-1099437", id = "ieu-a-2")
dim(data1_ieu)
ss1_vcf=create_summaryset(data1_vcf, type="vcf")
ss1_vcf=create_summaryset_from_gwasvcf(data1_vcf)
ss1_ieu=create_summaryset(data1_ieu, type="tibble")


data2_vcf <- gwasvcf::query_gwas(vcf, chrompos=c("5:74132993-75132993"))%>% 
gwasvcf::vcf_to_granges() %>% dplyr::as_tibble()
dim(data2_vcf)
ss2_vcf=create_summaryset(data2_vcf, type="vcf",build="GRCh37")

sumset_vcf <- gwasvcf::query_gwas(vcf, chrompos=c("5:74132993-75132993"))%>%           gwasvcf::gwasvcf_to_summaryset()

a=ss2_vcf %>% 
  liftover(., to="GRCh38")

data1 <- ieugwasr::associations(variants = "5:74132993-75132993", id = "ukb-d-I9_IHD")
 
 sumset1 <- create_summaryset(data1, build ="GRCh37")
 getMetadata(sumset1)

sumset1_lift <- download_chainfile(from = "GRCh37", to = "GRCh38") %>%

                liftover(sumset1,.)

 

sumset2_lift <- download_chainfile(from = "GRCh38", to = "GRCh37") %>%

                liftover(sumset1_lift,., to = "GRCh37")

# create_summaryset %>%

#     add_metadata %>%

#     convert_build(., chain)

 


 
# query_gwas %>%

#     gwasvcf_to_SummarySet()

 

# (variantAnnotation::read_vcf)

 

# gwasvcf::gwasvcf_to_SummarySet(vcf, ...) {

#     # get metadata from vcf

#     md <- vcf@metadata

 

#     # get summary data

#     tab <- vcf %>% gwasvcf_to_tibble()

 

#     # create summary set

#     ss <- create_summaryset_from_tibble(tab, ...)

 

#     # add metadata

#     ss <- ss %>% add_metadata(md)

#     return(ss)

# }

 

# ieugwasr::ieugwasr_to_SummarySet(dat, ...){

 

# }
 
