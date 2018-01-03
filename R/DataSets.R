#' sHiCDatum
#' 
#' Sample of HiC data from MES cells depleated of CTCF using auxin. Data obtained
#' from GEO, accession GSE98671, sample GSM2644947. These data are in sparse HiC format.
#' @format A sparseHiCdatum \describe{
#'     \item{sHicDatum@sampleName}{Here, auxin_rep1.}
#'     \item{sHicDatum@resolutionNamedList}{List of sparse contact matrices, indexed by resolution.}
#'     \item{sHicDatum@metaData}{Reference genome.}
#' }
"sHicDatum"

#' Auxin_rep1
#' 
#' Sample of HiC data from MES cells depleated of CTCF using auxin. Data obtained
#' from GEO, accession GSE98671, sample GSM2644947.
#' @format A list of two sparse matrices of class dtCMatrix \describe{ 
#'   \item{chr1}{Sample of data from chromosome one.} 
#'   \item{chr2}{Sample of data from chromosome two.} 
#'   }
#' 
"Auxin_rep1"

#' Auxin_rep2
#' 
#' Sample of HiC data from MES cells depleated of CTCF using auxin. Data obtained
#' from GEO, accession GSE98671, sample GSM2644948.
#' @format A list of two sparse matrices of class dtCMatrix \describe{ 
#'   \item{chr1}{Sample of data from chromosome one.}
#'   \item{chr2}{Sample of data from chromosome two.} 
#'   }
"Auxin_rep2"

#' Control_rep1
#' 
#' Sample of HiC data from control MES cells. Data obtained
#' from GEO, accession GSE98671, sample GSM2644945. 
#' @format A list of two sparse matrices of class dtCMatrix \describe{ 
#'   \item{chr1}{Sample of data from chromosome one.}
#'   \item{chr2}{Sample of data from chromosome two.} 
#'   }
"Control_rep1"

#' Control_rep2
#' 
#' Sample of HiC data from control MES cells. Data obtained
#' from GEO, accession GSE98671, sample GSM2644946. 
#' @format A list of two sparse matrices of class dtCMatrix \describe{ 
#'   \item{chr1}{Sample of data from chromosome one.}
#'   \item{chr2}{Sample of data from chromosome two.} 
#'   }
"Control_rep2"

#' Normalized HiC Experiment
#' 
#' A \code{DCSexp} containing pre-processed and normalized HiC data for control
#' (\code{\link{Control_rep1}}, \code{\link{Control_rep2}}) and auxin treated
#' (\code{\link{Auxin_rep1}}, \code{\link{Auxin_rep2}}) MES cells.
#' @format A \code{DCSexp} with two \code{DCSchr}s \describe{
#'     \item{\code{Z@Data[["1"]]}}{\code{DCSchr} with HiC data for chromosome 1.}
#'     \item{\code{Z@Data[["2"]]}}{\code{DCSchr} with HiC data for chromosome 2.}
#' }
"Z"

#' Border Scores
#' 
#' An \code{FSexp} containing border scores calculated using the normalized HiC
#' data in \code{\link{Z}}.
#' @format A \code{FSexp} with two \code{FSchr}s \describe{
#'     \item{\code{B@Data[["1"]]}}{\code{FSchr} with border scores for chromosome 1.}
#'     \item{\code{B@Data[["2"]]}}{\code{FSchr} with border scores  for chromosome 2.} 
#' }
"B"
