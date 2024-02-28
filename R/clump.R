#' Wrapper for clump function using local plink binary and ld reference dataset (based on ieugwasr::ld_clump_local)
#'
#' @param dat Dataframe. Must have a variant name column (`variant`) and pval column called `pval`. 
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
ld_clump <- function(dat, clump_kb=10000, clump_r2=0.001, clump_p=1, bfile=NULL, plink_bin=NULL)
{

	# Make textfile
	shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
	fn <- tempfile()
	write.table(data.frame(SNP=dat[["rsid"]], P=dat[["pval"]]), file=fn, row.names=F, col.names=T, quote=F)

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
	y <- subset(dat, !dat[["rsid"]] %in% res[["SNP"]])
	if(nrow(y) > 0)
	{
		message("Removing ", length(y[["rsid"]]), " of ", nrow(dat), " variants due to LD with other variants or absence from LD reference panel")
	}
	return(subset(dat, dat[["rsid"]] %in% res[["SNP"]]))
}
