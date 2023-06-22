#' modified from ieugwasr::ld_matrix_local ()
#' Get LD matrix using local plink binary and reference dataset
#'
#' @param variants List of variants (in plink 'set range' format). 
#' @param bfile Path to bed/bim/fam ld reference panel
#' @param plink_bin Specify path to plink binary. Default = `NULL`. 
#' @importFrom utils read.table write.table
#'
#' @return data frame
ld_matrix_local <- function(variants, bfile, plink_bin)
{
	# Make textfile
	shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
	fn <- tempfile()
	write.table(data.frame(variants), file=fn, row.names=F, col.names=F, quote=F)

	
	fun1 <- paste0(
		shQuote(plink_bin, type=shell),
		" --bfile ", shQuote(bfile, type=shell),
		" --extract range ", shQuote(fn, type=shell), 
		" --make-just-bim ", 
		" --keep-allele-order ",
		" --out ", shQuote(fn, type=shell)
	)
	system(fun1)

	bim <- read.table(paste0(fn, ".bim"), stringsAsFactors=FALSE)

	fun2 <- paste0(
		shQuote(plink_bin, type=shell),
		" --bfile ", shQuote(bfile, type=shell),
		" --extract range ", shQuote(fn, type=shell), 
		" --r square ", 
		" --keep-allele-order ",
		" --out ", shQuote(fn, type=shell)
	)
	system(fun2)
	res <- read.table(paste0(fn, ".ld"), header=FALSE) %>% as.matrix

	# standardise LD matrix. When alleles are not in alphabetical order flip and multiply col/row by -1
	alleles <- bim[,5:6]
	alleles_sorted <- t(apply(alleles,1,sort))
	flip <- alleles[,1] !=alleles_sorted[,1]
	res[flip,flip] <- res[flip,flip] * -1
	diag(res) <- 1

	rownames(res) <- colnames(res) <- paste(create_variantid_plink(bim), bim$V5, bim$V6, sep="|")
	
	return(res)
}





create_variantid_plink <-function(bim_file) {
  bim <- bim_file
  nvariants <- dim(bim)[1] 
  variantid <- lapply(1:nvariants, function(i){
    x <- sort(c(bim[i,]$V5,bim[i,]$V6))
    if (nchar(x[1]) > 10 || nchar(x[2]) <= 10){
      id <- paste0(bim[i,]$V1,":", bim[i,]$V4,"_#",digest::digest(x[1],algo= "murmur32"),"_",x[2]) 
    }
    if (nchar(x[1]) <= 10 || nchar(x[2]) > 10){
      id <- paste0(bim[i,]$V1,":", bim[i,]$V4,"_",x[1],"_#",digest::digest(x[2],algo= "murmur32")) 
    }
    if (nchar(x[1]) > 10 || nchar(x[2]) > 10){
      id <- paste0(bim[i,]$V1,":", bim[i,]$V4,"_#",digest::digest(x[1],algo= "murmur32"),"_#",digest::digest(x[2],algo= "murmur32")) 
    } else {
      id <- paste0(bim[i,]$V1,":", bim[i,]$V4,"_",x[1],"_",x[2]) 
    }
  }) 

  variantid <- unlist(variantid)
  return(variantid)
}



#' Harmonise LD matrix against summary data (now it just looks for overlapef variants. The harmonisation is done in ld_matrix_local)
#' harmonise_ld_dat() is based TwoSampleMR::harmonise_ld_dat()
#'
#  LD matrix returns with variantid_ea_oa identifiers. Make sure that they are oriented to the same effect allele as the summary dataset. Summary dataset can be dat1 dataset or harmonised dartaset
#
#' @param x harmonised dataset
#' @param ld Output from ld_matrix
#' @return List of dataset and harmonised LD matrix
harmonise_ld_dat <- function(x, ld){
	snpnames <- do.call(rbind, strsplit(rownames(ld), split="\\|"))
	i1 <- snpnames[,1] %in% x$variantid
	ld <- ld[i1,i1]
	snpnames <- snpnames[i1,]
	i2 <- x$variantid %in% snpnames[,1]
	x <- x[i2,]
	# stopifnot(all(snpnames[,1] == x$variantid))
	x$ea <- as.character(x$ea)
	x$nea <- as.character(x$nea)
	# Set1 x and ld alleles match
	snpnames <- data.frame(snpnames, stringsAsFactors=FALSE)
	snpnames <- merge(subset(x, select=c(variantid, ea, nea)), snpnames, by.x="variantid", by.y="X1")
	snpnames <- snpnames[match(x$variantid, snpnames$variantid),]
	snpnames$keep <- (snpnames$X2 == snpnames$ea & snpnames$X3 == snpnames$nea) |
		(snpnames$X3 == snpnames$ea & snpnames$X2 == snpnames$nea)

	# What happens if everything is gone?
	if(nrow(x) == 0)
	{
		message(" - none of the SNPs could be aligned to the LD reference panel")
		return(NULL)
	}

	if(any(!snpnames$keep))
	{
		message(" - the following SNPs could not be aligned to the LD reference panel: \n- ", paste(subset(snpnames, !keep)$variantid, collapse="\n - "))
	}


	# snpnames$flip1 <- snpnames$X2 != snpnames$ea
	# x <- subset(x, variantid %in% snpnames$variantid)
	# temp1 <- x$ea[snpnames$flip1]
	# temp2 <- x$nea[snpnames$flip1]
	# x$beta[snpnames$flip1] <- x$beta[snpnames$flip1] * -1
	# x$ea[snpnames$flip1] <- temp2
	# x$nea[snpnames$flip1] <- temp1

	rownames(ld) <- snpnames$variantid
	colnames(ld) <- snpnames$variantid

	if(any(!snpnames$keep))
	{
		message("Removing ", sum(!snpnames$keep), " variants due to harmonisation issues")
		ld <- ld[snpnames$keep, snpnames$keep]
		x <- x[snpnames$keep, ]

	}
	return(list(x,ld))
}

