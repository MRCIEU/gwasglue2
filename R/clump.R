#' Wrapper for clump function using local plink binary and ld reference dataset (based on `ieugwasr::ld_clump_local()`)
#'
#' @param data Dataframe. Must have a variant name column (`variant`) and pval column called `p`. 
#' If `id` is present then clumping will be done per unique id.
#' @param clump_kb Clumping kb window. Default is very strict, `10000`
#' @param clump_r2 Clumping r2 threshold. Default is very strict, `0.001`
#' @param clump_p Clumping sig level for index variants. Default = `1` (i.e. no threshold)
#' @param bfile If this is provided then will use the API. Default = `NULL`
#' @param plink_bin Specify path to plink binary. Default = `NULL`. 
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @export
#' @return data frame of clumped variants
ld_clump <- function(data, clump_kb=10000, clump_r2=0.001, clump_p=1, bfile=NULL, plink_bin=NULL)
{

	# Make textfile
	shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
	fn <- tempfile()
	write.table(data.frame(SNP=data[["rsid"]], P=data[["p"]]), file=fn, row.names=F, col.names=T, quote=F)

	fun2 <- paste0(
		shQuote(plink_bin, type=shell),
		" --bfile ", shQuote(bfile, type=shell),
		" --clump ", shQuote(fn, type=shell), 
		" --clump-p1 ", clump_p, 
		" --clump-r2 ", clump_r2, 
		" --clump-kb ", clump_kb, 
		" --out ", shQuote(fn, type=shell)
	)
	system(fun2)
	res <- read.table(paste(fn, ".clumped", sep=""), header=T)
	unlink(paste(fn, "*", sep=""))
	y <- subset(data, !data[["rsid"]] %in% res[["SNP"]])
	if(nrow(y) > 0)
	{
		message("Removing ", length(y[["rsid"]]), " of ", nrow(data), " variants due to LD with other variants or absence from LD reference panel")
	}
	return(subset(data, data[["rsid"]] %in% res[["SNP"]]))
}

#' Extract top hits variants from a data frame
#' @param data GWAS summary data. Dataframe.
#' @param pval P-value threshold. Default = `5e-08`
#' @param clump Logical. If `TRUE` then clump the data. Default = `TRUE`
#' @param r2 Clumping r2 threshold. Default = `0.001`
#' @param kb Clumping kb window. Default = `10000`
#' @param variant_col Name of the variant column. Default = `rsid`
#' @param pval_col Name of the p-value column. Default = `p`
#' @param bfile It corresponds to the path and prefix of the plink files used to build the LD correlation matrix. 
#' @param plink_bin Path to the plink executable
get_tophits_from_data <- function(data, pval = 5e-08, clump = TRUE, r2 = 0.001, kb = 10000, variant_col="rsid", pval_col = "p", bfile=NULL, plink_bin=NULL){

	if(clump)
	{
		return(ld_clump(data, clump_kb=kb, clump_r2=r2, clump_p=pval, bfile=bfile, plink_bin=plink_bin))
	}
	else
	{
		return(subset(data, data[[pval_col]] < pval))
	}


}
