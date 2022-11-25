# WIP
# harmonise_ld_dat() is based TwoSampleMR::harmonise_ld_dat()


#' Harmonise LD matrix against summary data
#'
#' LD matrix returns with rsid_ea_oa identifiers. Make sure that they are oriented to the same effect allele as the summary dataset. Summary dataset can be dat1 dataset or harmonised dartaset
#'
#' @param x harmonised dataset
#' @param ld Output from ld_matrix
#' @return List of dataset and harmonised LD matrix
#'
#'
harmonise_ld_dat <- function(x, ld)
{
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



# ##################################################
# # WIP

# #' Harmonise dat1 and dat2 for multivariable MR
# #'
# #' @md
# #' @param dat1 Output from [`mv_extract_dat1s`].
# #' @param dat2 Output from `extract_dat2a(dat1$rsid, id_output)`.
# #' @param harmonise_strictness See the `action` option of [`harmonise_data`]. The default is `2`.
# #'
# #' @export
# #' @return List of vectors and matrices required for mv analysis.
# #' \describe{
# #' \item{dat1_beta}{a matrix of beta coefficients, in which rows correspond to SNPs and columns correspond to dat1s.}
# #' \item{dat1_se}{is the same as `dat1_beta`, but for standard errors.}
# #' \item{dat1_pval}{the same as `dat1_beta`, but for p-values.}
# #' \item{expname}{A data frame with two variables, `id.dat1` and `dat1` which are character strings.}
# #' \item{dat2_beta}{an array of effects for the dat2, corresponding to the SNPs in dat1_beta.}
# #' \item{dat2_se}{an array of standard errors for the dat2.}
# #' \item{dat2_pval}{an array of p-values for the dat2.}
# #' \item{outname}{A data frame with two variables, `id.dat2` and `dat2` which are character strings.}
# #' }
# #'
# mv_harmonise_data <- function(dat1, dat2, harmonise_strictness=2)
# {

#   stopifnot(all(c("rsid", "id", "ea", "beta", "se", "p") %in% names(dat1)))
#   nexp <- length(unique(dat1$id))
#   stopifnot(nexp > 1)
#   tab <- table(dat1$rsid)
#   keepsnp <- names(tab)[tab == nexp]
#   dat1 <- subset(dat1, rsid %in% keepsnp)


#   dat1_mat <- reshape2::dcast(dat1, rsid ~ id, value.var="beta")


#   # Get dat2 data
#   dat <- harmonise_data(subset(dat1, id == dat1$id[1]), dat2, action=harmonise_strictness) #TODO strictness, harmonise data
#   dat <- subset(dat, mr_keep) #TODO mrkeep
#   dat$rsid <- as.character(dat$rsid)

#   dat1_beta <- reshape2::dcast(dat1, rsid ~ id, value.var="beta")
#   dat1_beta <- subset(dat1_beta, rsid %in% dat$rsid)
#   dat1_beta$rsid <- as.character(dat1_beta$rsid)

#   dat1_pval <- reshape2::dcast(dat1, rsid ~ id, value.var="p")
#   dat1_pval <- subset(dat1_pval, rsid %in% dat$rsid)
#   dat1_pval$rsid <- as.character(dat1_pval$rsid)

#   dat1_se <- reshape2::dcast(dat1, rsid ~ id, value.var="se")
#   dat1_se <- subset(dat1_se, rsid %in% dat$rsid)
#   dat1_se$rsid <- as.character(dat1_se$rsid)

#   index <- match(dat1_beta$rsid, dat$rsid)
#   dat <- dat[index, ]
#   stopifnot(all(dat$rsid == dat1_beta$rsid))

#   dat1_beta <- as.matrix(dat1_beta[,-1])
#   dat1_pval <- as.matrix(dat1_pval[,-1])
#   dat1_se <- as.matrix(dat1_se[,-1])

#   rownames(dat1_beta) <- dat$rsid
#   rownames(dat1_pval) <- dat$rsid
#   rownames(dat1_se) <- dat$rsid

#   dat2_beta <- dat$beta
#   dat2_se <- dat$se
#   dat2_pval <- dat$pval

#   return(list(dat1,dat))
# }

