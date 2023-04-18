# WIP
# harmonise_ld_dat() is based TwoSampleMR::harmonise_ld_dat()


#Harmonise LD matrix against summary data
#
# LD matrix returns with rsid_ea_oa identifiers. Make sure that they are oriented to the same effect allele as the summary dataset. Summary dataset can be dat1 dataset or harmonised dartaset
#
# @param x harmonised dataset
# @param ld Output from ld_matrix
# @return List of dataset and harmonised LD matrix
harmonise_ld_dat <- function(x, ld){
	snpnames <- do.call(rbind, strsplit(rownames(ld), split="_"))
	i1 <- snpnames[,1] %in% x$rsid
	ld <- ld[i1,i1]
	snpnames <- snpnames[i1,]
	i2 <- x$rsid %in% snpnames[,1]
	x <- x[i2,]
	# stopifnot(all(snpnames[,1] == x$rsid))
	x$ea <- as.character(x$ea)
	x$nea <- as.character(x$nea)
	# Set1 x and ld alleles match
	snpnames <- data.frame(snpnames, stringsAsFactors=FALSE)
	snpnames <- merge(subset(x, select=c(rsid, ea, nea)), snpnames, by.x="rsid", by.y="X1")
	snpnames <- snpnames[match(x$rsid, snpnames$rsid),]
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
		message(" - the following SNPs could not be aligned to the LD reference panel: \n- ", paste(subset(snpnames, !keep)$rsid, collapse="\n - "))
	}


	snpnames$flip1 <- snpnames$X2 != snpnames$ea
	x <- subset(x, rsid %in% snpnames$rsid)
	temp1 <- x$ea[snpnames$flip1]
	temp2 <- x$nea[snpnames$flip1]
	x$beta[snpnames$flip1] <- x$beta[snpnames$flip1] * -1
	x$ea[snpnames$flip1] <- temp2
	x$nea[snpnames$flip1] <- temp1

	rownames(ld) <- snpnames$rsid
	colnames(ld) <- snpnames$rsid

	if(any(!snpnames$keep))
	{
		message("Removing ", sum(!snpnames$keep), " variants due to harmonisation issues")
		ld <- ld[snpnames$keep, snpnames$keep]
		x <- x[snpnames$keep, ]

	}
	return(list(x,ld))
}

