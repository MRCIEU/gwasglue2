#' Harmonise the alleles and effects between two summary sets
#' This is based o the harmonise function of TwoSampleMR



# Changes:
# 1:instead of having one dataset (dat), we have the variables needed
# 2:
# 3: change .exposure and .outcome columns names in dat
# 4: columns names is the same as in IEUopenGWAS


#' @param action Level of strictness in dealing with SNPs.
#' * `action = 1`: Assume all alleles are coded on the forward strand, i.e. do not attempt to flip alleles
#' * `action = 2`: Try to infer positive strand alleles, using allele frequencies for palindromes (default, conservative);
#' * `action = 3`: Correct strand for non-palindromic SNPs, and drop all palindromic SNPs from the analysis (more conservative).

#'
harmonise <- function(rsid, A1, A2, B1, B2, betaA, betaB, fA, fB, tolerance, action)
{

  orig_SNP<-rsid

  SNP_index<-sapply(1:length(rsid),function(i)sum(rsid[1:i]==rsid[i]))
  rsid <- paste0(rsid,"_",SNP_index)


  i22 <- !is.na(A1) & !is.na(A2) & !is.na(B1) & !is.na(B2)
  i21 <- !is.na(A1) & !is.na(A2) & !is.na(B1) & is.na(B2)
  i12 <- !is.na(A1) & is.na(A2) & !is.na(B1) & !is.na(B2)
  i11 <- !is.na(A1) & is.na(A2) & !is.na(B1) & is.na(B2)

  d22 <- harmonise_22(rsid[i22], A1[i22], A2[i22], B1[i22], B2[i22], betaA[i22], betaB[i22], fA[i22], fB[i22], tolerance, action)
  d21 <- harmonise_21(rsid[i21], A1[i21], A2[i21], B1[i21], betaA[i21], betaB[i21], fA[i21], fB[i21], tolerance, action)
  d12 <- harmonise_12(rsid[i12], A1[i12], B1[i12], B2[i12], betaA[i12], betaB[i12], fA[i12], fB[i12], tolerance, action)
  d11 <- harmonise_11(rsid[i11], A1[i11], B1[i11], betaA[i11], betaB[i11], fA[i11], fB[i11], tolerance, action)

  # jlog <- dplyr::bind_rows(
  # 	as.data.frame(attr(d22, "log"), stringsAsFactors=FALSE),
  # 	as.data.frame(attr(d21, "log"), stringsAsFactors=FALSE),
  # 	as.data.frame(attr(d12, "log"), stringsAsFactors=FALSE),
  # 	as.data.frame(attr(d11, "log"), stringsAsFactors=FALSE)
  # )
  # jlog <- cbind(data.frame(id.exposure=dat$id.exposure[1], id.outcome=dat$id.outcome[1], stringsAsFactors=FALSE), jlog)

  d <- rbind(d21, d22, d12, d11)

  d$rsid <- orig_SNP
  d$keep <- TRUE

  if(action == 3)
  {
    d$keep[d$palindromic | d$remove | d$ambiguous] <- FALSE
    if(any(d$palindromic))
    {
      message("Removing the following SNPs for being palindromic:\n", paste(d$rsid[d$palindromic], collapse=", "))
    }
    if(any(d$remove))
    {
      message("Removing the following SNPs for incompatible alleles:\n", paste(d$rsid[d$remove], collapse=", "))
    }
    # jlog[['incompatible_alleles']] <- sum(d$remove)
    if(any(d$ambiguous & !d$palindromic)) #TODO I think it should be | instead of & or just any(d$ambiguous)
    {
      message("Removing the following SNPs for having incompatible allele frequencies:\n", paste(d$rsid[d$ambiguous], collapse=", "))
    }
    # jlog[['ambiguous_alleles']] <- sum(d$ambiguous)
  }
  if(action == 2)
  {
    d$keep[d$remove | d$ambiguous] <- FALSE
    if(any(d$remove))
    {
      message("Removing the following SNPs for incompatible alleles:\n", paste(d$rsid[d$remove], collapse=", "))
    }
    # jlog[['incompatible_alleles']] <- sum(d$remove)
    if(any(d$ambiguous))
    {
      message("Removing the following SNPs for being palindromic with intermediate allele frequencies:\n", paste(d$rsid[d$ambiguous], collapse=", "))
    }
    # jlog[['ambiguous_alleles']] <- sum(d$ambiguous)
  }
  if(action == 1)
  {
    d$keep[d$remove] <- FALSE
    if(any(d$remove))
    {
      message("Removing the following SNPs for incompatible alleles:\n", paste(d$rsid[d$remove], collapse=", "))
    }else{message("No alleles were removed")}
    # jlog[['incompatible_alleles']] <- sum(d$remove)
  }

  # if(nrow(d1) != nrow(d))
  # {
  # 	message("Removing the following SNPs due to harmonising issues:\n",
  # 		paste(subset(d, ! SNP %in% d1$SNP)$SNP, collapse="\n")
  # 	)
  # }


  # attr(d, "log") <- jlog
  return(d)
}


check_palindromic <- function(A1, A2)
{
	(A1 == "T" & A2 == "A") |
	(A1 == "A" & A2 == "T") |
	(A1 == "G" & A2 == "C") |
	(A1 == "C" & A2 == "G")
}


flip_alleles <- function(x)
{
	x <- toupper(x)
	x <- gsub("C", "g", x)
	x <- gsub("G", "c", x)
	x <- gsub("A", "t", x)
	x <- gsub("T", "a", x)
	return(toupper(x))
}


recode_indels_22 <- function(A1, A2, B1, B2)
{

	ncA1 <- nchar(A1)
	ncA2 <- nchar(A2)
	ncB1 <- nchar(B1)
	ncB2 <- nchar(B2)


	i1 <- ncA1 > ncA2 & B1 == "I" & B2 == "D"
	B1[i1] <- A1[i1]
	B2[i1] <- A2[i1]

	i1 <- ncA1 < ncA2 & B1 == "I" & B2 == "D"
	B1[i1] <- A2[i1]
	B2[i1] <- A1[i1]

	i1 <- ncA1 > ncA2 & B1 == "D" & B2 == "I"
	B1[i1] <- A2[i1]
	B2[i1] <- A1[i1]

	i1 <- ncA1 < ncA2 & B1 == "D" & B2 == "I"
	B1[i1] <- A1[i1]
	B2[i1] <- A2[i1]


	i1 <- ncB1 > ncB2 & A1 == "I" & A2 == "D"
	A1[i1] <- B1[i1]
	A2[i1] <- B2[i1]

	i1 <- ncB1 < ncB2 & A1 == "I" & A2 == "D"
	A2[i1] <- B1[i1]
	A1[i1] <- B2[i1]

	i1 <- ncB1 > ncB2 & A1 == "D" & A2 == "I"
	A2[i1] <- B1[i1]
	A1[i1] <- B2[i1]

	i1 <- ncB1 < ncB2 & A1 == "D" & A2 == "I"
	A1[i1] <- B1[i1]
	A2[i1] <- B2[i1]

	keep <- rep(TRUE, length(A1))
	keep[(ncA1 > 1 & ncA1 == ncA2 & (B1 == "D" | B1 == "I"))] <- FALSE
	keep[(ncB1 > 1 & ncB1 == ncB2 & (A1 == "D" | A1 == "I"))] <- FALSE
	keep[A1 == A2] <- FALSE
	keep[B1 == B2] <- FALSE

	return(data.frame(A1=A1, A2=A2, B1=B1, B2=B2, keep=keep, stringsAsFactors=FALSE))
}


recode_indels_21 <- function(A1, A2, B1)
{
	ncA1 <- nchar(A1)
	ncA2 <- nchar(A2)
	ncB1 <- nchar(B1)

	B2 <- rep(NA, length(B1))

	i1 <- ncA1 > ncA2 & B1 == "I"
	B1[i1] <- A1[i1]
	B2[i1] <- A2[i1]

	i1 <- ncA1 < ncA2 & B1 == "I"
	B1[i1] <- A2[i1]
	B2[i1] <- A1[i1]

	i1 <- ncA1 > ncA2 & B1 == "D"
	B1[i1] <- A2[i1]
	B2[i1] <- A1[i1]

	i1 <- ncA1 < ncA2 & B1 == "D"
	B1[i1] <- A1[i1]
	B2[i1] <- A2[i1]

	keep <- rep(TRUE, length(A1))
	keep[A1 == "I" & A2 == "D"] <- FALSE
	keep[A1 == "D" & A2 == "I"] <- FALSE
	keep[(ncA1 > 1 & ncA1 == ncA2 & (B1 == "D" | B1 == "I"))] <- FALSE
	keep[A1 == A2] <- FALSE

	return(data.frame(A1=A1, A2=A2, B1=B1, B2=B2, keep=keep, stringsAsFactors=FALSE))
}


recode_indels_12 <- function(A1, B1, B2)
{
	ncA1 <- nchar(A1)
	ncB1 <- nchar(B1)
	ncB2 <- nchar(B2)

	A2 <- rep(NA, length(A1))

	i1 <- ncB1 > ncB2 & A1 == "I"
	A1[i1] <- B1[i1]
	A2[i1] <- B2[i1]

	i1 <- ncB1 < ncB2 & A1 == "I"
	A2[i1] <- B1[i1]
	A1[i1] <- B2[i1]

	i1 <- ncB1 > ncB2 & A1 == "D"
	A2[i1] <- B1[i1]
	A1[i1] <- B2[i1]

	i1 <- ncB1 < ncB2 & A1 == "D"
	A1[i1] <- B1[i1]
	A2[i1] <- B2[i1]

	keep <- rep(TRUE, length(A1))
	keep[B1 == "I" & B2 == "D"] <- FALSE
	keep[B1 == "D" & B2 == "I"] <- FALSE
	keep[(ncB1 > 1 & ncB1 == ncB2 & (A1 == "D" | A1 == "I"))] <- FALSE
	keep[B1 == B2] <- FALSE

	return(data.frame(A1=A1, A2=A2, B1=B1, B2=B2, keep=keep, stringsAsFactors=FALSE))
}



harmonise_22 <- function(SNP, A1, A2, B1, B2, betaA, betaB, fA, fB, tolerance, action)
{
	if(length(SNP) == 0) return(data.frame())
	jlog <- list()
	jlog[['alleles']] <- "2-2"

	indel_index <- nchar(A1) > 1 | nchar(A2) > 1 | A1 == "D" | A1 == "I"
	temp <- recode_indels_22(A1[indel_index], A2[indel_index], B1[indel_index], B2[indel_index])

	A1[indel_index] <- temp$A1
	A2[indel_index] <- temp$A2
	B1[indel_index] <- temp$B1
	B2[indel_index] <- temp$B2


	# Find SNPs with alleles that match in A and B
	status1 <- (A1 == B1) & (A2 == B2)
	to_swap <- (A1 == B2) & (A2 == B1)
	jlog[['switched_alleles']] <- sum(to_swap)


	# If B's alleles are the wrong way round then swap
	Btemp <- B1[to_swap]
	B1[to_swap] <- B2[to_swap]
	B2[to_swap] <- Btemp
	betaB[to_swap] <- betaB[to_swap] * -1
	fB[to_swap] <- 1 - fB[to_swap]

	# Check again
	status1 <- (A1 == B1) & (A2 == B2)
	palindromic <- check_palindromic(A1, A2)

	# If NOT palindromic and alleles DON'T match then try flipping
	i <- !palindromic & !status1
	B1[i] <- flip_alleles(B1[i])
	B2[i] <- flip_alleles(B2[i])
	status1 <- (A1 == B1) & (A2 == B2)
	jlog[['flipped_alleles_basic']] <- sum(i)

	# If still NOT palindromic and alleles DON'T match then try swapping
	i <- !palindromic & !status1
	to_swap <- (A1 == B2) & (A2 == B1)
	Btemp <- B1[to_swap]
	B1[to_swap] <- B2[to_swap]
	B2[to_swap] <- Btemp
	betaB[to_swap] <- betaB[to_swap] * -1
	fB[to_swap] <- 1 - fB[to_swap]

	# Any SNPs left with unmatching alleles need to be removed
	status1 <- (A1 == B1) & (A2 == B2)
	remove <- !status1
	remove[indel_index][!temp$keep] <- TRUE

	# Now deal with palindromic SNPs
	# If the frequency is within tolerance then they are ambiguous and need to be flagged
	minf <- 0.5 - tolerance
	maxf <- 0.5 + tolerance
	tempfA <- fA
	tempfB <- fB
	tempfA[is.na(tempfA)] <- 0.5
	tempfB[is.na(tempfB)] <- 0.5
	ambiguousA <- tempfA > minf & tempfA < maxf
	ambiguousB <- tempfB > minf & tempfB < maxf

	# If action = 2 (flip alleles based on frequency) then flip and swap
	if(action == 2)
	{
		status2 <- ((tempfA < 0.5 & tempfB > 0.5) | (tempfA > 0.5 & tempfB < 0.5)) & palindromic
		to_swap <- status2 & !remove
		betaB[to_swap] <- betaB[to_swap] * -1
		fB[to_swap] <- 1 - fB[to_swap]
		if(!is.null(jlog))
		{
			jlog[['flipped_alleles_palindrome']] <- sum(to_swap)
		}
	} else {
		if(!is.null(jlog))
		{
			jlog[['flipped_alleles_palindrome']] <- 0
		}
	}
	#TODO
	d <- data.frame(rsid=SNP,A1=A1, A2=A2, B1=B1, B2=B2, betaA=betaA, betaB=betaB, fA=fA, fB=fB, remove=remove, palindromic=palindromic, ambiguous=(ambiguousA|ambiguousB) & palindromic)
	# d <- data.frame(SNP=SNP, effect_allele.exposure=A1, other_allele.exposure=A2, effect_allele.outcome=B1, other_allele.outcome=B2, beta.exposure=betaA, beta.outcome=betaB, eaf.exposure=fA, eaf.outcome=fB, remove=remove, palindromic=palindromic, ambiguous=(ambiguousA|ambiguousB) & palindromic)
	attr(d, "log") <- jlog
	return(d)
}

harmonise_21 <- function(SNP, A1, A2, B1, betaA, betaB, fA, fB, tolerance, action)
{
	if(length(SNP) == 0) return(data.frame())
	jlog <- list()
	jlog[['alleles']] <- "2-1"

	n <- length(A1)
	B2 <- rep(NA, n)
	ambiguous <- rep(FALSE, n)
	palindromic <- check_palindromic(A1, A2)
	remove <- palindromic


	indel_index <- nchar(A1) > 1 | nchar(A2) > 1 | A1 == "D" | A1 == "I"
	temp <- recode_indels_21(A1[indel_index], A2[indel_index], B1[indel_index])

	A1[indel_index] <- temp$A1
	A2[indel_index] <- temp$A2
	B1[indel_index] <- temp$B1
	B2[indel_index] <- temp$B2
	remove[indel_index][!temp$keep] <- TRUE


	status1 <- A1 == B1
	minf <- 0.5 - tolerance
	maxf <- 0.5 + tolerance

	tempfA <- fA
	tempfB <- fB
	tempfA[is.na(tempfA)] <- 0.5
	tempfB[is.na(tempfB)] <- 0.5

	freq_similar1 <- (tempfA < minf & tempfB < minf) | (tempfA > maxf & tempfB > maxf)
	ambiguous[status1 & !freq_similar1] <- TRUE

	B2[status1] <- A2[status1]

	to_swap <- A2 == B1
	jlog[['switched_alleles']] <- sum(to_swap)
	freq_similar2 <- (tempfA < minf & tempfB > maxf) | (tempfA > maxf & tempfB < minf)

	ambiguous[to_swap & !freq_similar2] <- TRUE
	B2[to_swap] <- B1[to_swap]
	B1[to_swap] <- A1[to_swap]
	betaB[to_swap] <- betaB[to_swap] * -1
	fB[to_swap] <- 1 - fB[to_swap]

	to_flip <- A1 != B1 & A2 != B1
	jlog[['flipped_alleles_no_oa']] <- sum(to_flip)

	ambiguous[to_flip] <- TRUE

	B1[to_flip] <- flip_alleles(B1[to_flip])
	status1 <- A1 == B1
	B2[status1] <- A2[status1]

	to_swap <- A2 == B1
	B2[to_swap] <- B1[to_swap]
	B1[to_swap] <- A1[to_swap]
	betaB[to_swap] <- betaB[to_swap] * -1
	fB[to_swap] <- 1 - fB[to_swap]

	#TODO
	d <- data.frame(rsid=SNP,A1=A1, A2=A2, B1=B1, B2=B2, betaA=betaA, betaB=betaB, fA=fA, fB=fB, remove=remove, palindromic=palindromic, ambiguous=(ambiguousA|ambiguousB) & palindromic)
# d <- data.frame(SNP=SNP, effect_allele.exposure=A1, other_allele.exposure=A2, effect_allele.outcome=B1, other_allele.outcome=B2, beta.exposure=betaA, beta.outcome=betaB, eaf.exposure=fA, eaf.outcome=fB, remove=remove, palindromic=palindromic, ambiguous=ambiguous | palindromic)

	attr(d, "log") <- jlog
	return(d)

}

harmonise_12 <- function(SNP, A1, B1, B2, betaA, betaB, fA, fB, tolerance, action)
{
	if(length(SNP) == 0) return(data.frame())
	jlog <- list()
	jlog[['alleles']] <- "1-2"

	n <- length(A1)
	A2 <- rep(NA, n)
	ambiguous <- rep(FALSE, n)
	palindromic <- check_palindromic(B1, B2)
	remove <- palindromic

	indel_index <- nchar(B1) > 1 | nchar(B2) > 1 | B1 == "D" | B1 == "I"
	temp <- recode_indels_21(A1[indel_index], B1[indel_index], B2[indel_index])

	A1[indel_index] <- temp$A1
	A2[indel_index] <- temp$A2
	B1[indel_index] <- temp$B1
	B2[indel_index] <- temp$B2
	remove[indel_index][!temp$keep] <- TRUE

	status1 <- A1 == B1
	minf <- 0.5 - tolerance
	maxf <- 0.5 + tolerance

	tempfA <- fA
	tempfB <- fB
	tempfA[is.na(tempfA)] <- 0.5
	tempfB[is.na(tempfB)] <- 0.5

	freq_similar1 <- (tempfA < minf & tempfB < minf) | (tempfA > maxf & tempfB > maxf)
	ambiguous[status1 & !freq_similar1] <- TRUE

	A2[status1] <- B2[status1]

	to_swap <- A1 == B2
	jlog[['switched_alleles']] <- sum(to_swap)

	freq_similar2 <- (tempfA < minf & tempfB > maxf) | (tempfA > maxf & tempfB < minf)

	ambiguous[to_swap & !freq_similar2] <- TRUE
	A2[to_swap] <- A1[to_swap]
	A1[to_swap] <- B1[to_swap]
	betaA[to_swap] <- betaA[to_swap] * -1
	fA[to_swap] <- 1 - fA[to_swap]

	to_flip <- A1 != B1 & A1 != B2
	jlog[['flipped_alleles_no_oa']] <- sum(to_flip)

	ambiguous[to_flip] <- TRUE

	A1[to_flip] <- flip_alleles(A1[to_flip])
	status1 <- A1 == B1
	A2[status1] <- B2[status1]

	to_swap <- B2 == A1
	B2[to_swap] <- B1[to_swap]
	B1[to_swap] <- A1[to_swap]
	betaB[to_swap] <- betaB[to_swap] * -1
	fB[to_swap] <- 1 - fB[to_swap]

	#TODO
	d <- data.frame(rsid=SNP,A1=A1, A2=A2, B1=B1, B2=B2, betaA=betaA, betaB=betaB, fA=fA, fB=fB, remove=remove, palindromic=palindromic, ambiguous=(ambiguousA|ambiguousB) & palindromic)
	# d <- data.frame(SNP=SNP, effect_allele.exposure=A1, other_allele.exposure=A2, effect_allele.outcome=B1, other_allele.outcome=B2, beta.exposure=betaA, beta.outcome=betaB, eaf.exposure=fA, eaf.outcome=fB, remove=remove, palindromic=palindromic, ambiguous=ambiguous | palindromic)
	attr(d, "log") <- jlog

	return(d)

}


harmonise_11 <- function(SNP, A1, B1, betaA, betaB, fA, fB, tolerance, action)
{
	if(length(SNP) == 0) return(data.frame())
	jlog <- list()
	jlog[['alleles']] <- "1-1"

	n <- length(A1)
	A2 <- rep(NA, n)
	B2 <- rep(NA, n)
	ambiguous <- rep(FALSE, n)
	palindromic <- FALSE

	status1 <- A1 == B1
	remove <- !status1

	minf <- 0.5 - tolerance
	maxf <- 0.5 + tolerance

	tempfA <- fA
	tempfB <- fB
	tempfA[is.na(tempfA)] <- 0.5
	tempfB[is.na(tempfB)] <- 0.5

	freq_similar1 <- (tempfA < minf & tempfB < minf) | (tempfA > maxf & tempfB > maxf)
	ambiguous[status1 & !freq_similar1] <- TRUE

	#TODO
	d <- data.frame(rsid=SNP,A1=A1, A2=A2, B1=B1, B2=B2, betaA=betaA, betaB=betaB, fA=fA, fB=fB, remove=remove, palindromic=palindromic, ambiguous=(ambiguousA|ambiguousB) & palindromic)
	# d <- data.frame(SNP=SNP, effect_allele.exposure=A1, other_allele.exposure=A2, effect_allele.outcome=B1, other_allele.outcome=B2, beta.exposure=betaA, beta.outcome=betaB, eaf.exposure=fA, eaf.outcome=fB, remove=remove, palindromic=palindromic, ambiguous=ambiguous | palindromic)

	attr(d, "log") <- jlog
	return(d)
}

###############################################
# The code below is not being used

harmonise_data <- function(dat1, dat2, action=2)
{
  stopifnot(all(action %in% 1:3))
  res.tab <- merge(dat2, dat1, by="rsid") #TODO if I merge the data there is going to be columns with the same name
  ncombinations <- length(unique(res.tab$id.outcome))
  if(length(action) == 1)
  {
    action <- rep(action, ncombinations)
  } else if(length(action) != ncombinations) {
    stop("Action argument must be of length 1 (where the same action will be used for all outcomes), or number of unique id.outcome values (where each outcome will use a different action value)")
  }

  res.tab <- harmonise_cleanup_variables(res.tab)
  if (nrow(res.tab) == 0) {
    # NOTE: currently the columns are not in the same order as they would
    # be should there be overlapping variants
    return(res.tab)
  }

  d <- data.frame(id.outcome=unique(res.tab$id.outcome), action=action)
  res.tab <- merge(res.tab, d, by="id.outcome")

  combs <- subset(res.tab, !duplicated(paste(id.exposure, id.outcome)), select=c(id.exposure, id.outcome))

  fix.tab <- list()
  mr_cols <- c("beta.exposure", "beta.outcome", "se.exposure", "se.outcome")
  for(i in 1:nrow(combs))
  {
    x <- subset(res.tab, id.exposure == combs$id.exposure[i] & id.outcome == combs$id.outcome[i])
    message("Harmonising ", x$exposure[1], " (", x$id.exposure[1], ") and ", x$outcome[1], " (", x$id.outcome[1], ")")
    x <- harmonise(x, 0.08, x$action[1])
    attr(x, "log")[['candidate_variants']] <- sum(dat1$id.exposure == x$id.exposure[1])
    attr(x, "log")[['variants_absent_from_reference']] <- sum(dat1$id.exposure == x$id.exposure[1]) - nrow(x)

    x$mr_keep[apply(x[, mr_cols], 1, function(y) any(is.na(y)))] <- FALSE
    attr(x, "log")[["total_variants"]] <- nrow(x)
    attr(x, "log")[["total_variants_for_mr"]] <- sum(x$mr_keep)
    attr(x, "log")[["proxy_variants"]] <- ifelse(is.null(x$proxy.outcome), 0, sum(x$proxy.outcome, na.rm=TRUE))
    fix.tab[[i]] <- x
  }
  # fix.tab <- plyr::dlply(res.tab, c("id.exposure", "id.outcome"), function(x)
  # {
  # 	x <- plyr::mutate(x)
  # 	message("Harmonising ", x$exposure[1], " (", x$id.exposure[1], ") and ", x$outcome[1], " (", x$id.outcome[1], ")")

  # 	return(x)
  # })

  jlog <- dplyr::bind_rows(lapply(fix.tab, function(x) attr(x, "log")))
  fix.tab <- dplyr::bind_rows(fix.tab)
  attr(fix.tab, "log") <- jlog

  # fix.tab <- harmonise_make_snp_effects_positive(fix.tab)
  if(!"samplesize.outcome" %in% names(fix.tab))
  {
    fix.tab$samplesize.outcome <- NA
  }


  return(fix.tab)
}


harmonise_cleanup_variables <- function(res.tab)
{
  res.tab$beta.exposure <- as.numeric(res.tab$beta.exposure)
  res.tab$beta.outcome <- as.numeric(res.tab$beta.outcome)
  res.tab$eaf.exposure <- as.numeric(res.tab$eaf.exposure)
  res.tab$eaf.outcome[res.tab$eaf.outcome=="NR"] <- NA
  res.tab$eaf.outcome[res.tab$eaf.outcome=="NR "] <- NA
  res.tab$eaf.outcome <- as.numeric(res.tab$eaf.outcome)

  # convert alleles to upper case
  res.tab$effect_allele.exposure <- toupper(res.tab$effect_allele.exposure)
  res.tab$other_allele.exposure <- toupper(res.tab$other_allele.exposure)
  res.tab$effect_allele.outcome <- toupper(res.tab$effect_allele.outcome)
  res.tab$other_allele.outcome <- toupper(res.tab$other_allele.outcome)
  res.tab$other_allele.outcome[res.tab$other_allele.outcome == ""] <- NA

  return(res.tab)
}


harmonise_make_snp_effects_positive <- function(res.tab)
{
  # code SNP effect so that effect allele is the allele that increases the trait
  #res.tab[,c("SNP","beta.exposure","beta.outcome","effect_allele.exposure","other_allele.exposure","effect_allele","other_allele","eaf.exposure","eaf.outcome","info_s1","RSQ_s2","RSQ_s3","info_s4","q_p.value")]

  pos.change <- which(res.tab$beta.exposure < 0)
  res.tab$eaf.exposure[pos.change] <- 1 - res.tab$eaf.exposure[pos.change]
  res.tab$beta.exposure[pos.change] <- res.tab$beta.exposure[pos.change] * -1
  eff.allele.change <- res.tab$effect_allele.exposure[pos.change]
  oth.allele.change <- res.tab$other_allele.exposure[pos.change]
  res.tab$effect_allele.exposure[pos.change] <- oth.allele.change
  res.tab$other_allele.exposure[pos.change] <- eff.allele.change

  return(res.tab)
}




check_required_columns <- function(dat, type="exposure")
{
	required_columns <- c(
		"SNP",
		paste0(c("id.", "", "beta.", "se.", "effect_allele.", "other_allele."), type)
	)
	index <- required_columns %in% names(dat)
	if(!all(index))
	{
		stop("The following required columns are missing from ", type, ": ", paste(required_columns[!index], collapse=", "))
	}
	return(NULL)
}
