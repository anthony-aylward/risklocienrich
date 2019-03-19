#===============================================================================
# risklocienrich.R
#===============================================================================

# Imports ======================================================================

#' @import BSgenome.Hsapiens.UCSC.hg19
#' @import GenomicRanges
#' @import gwascat
#' @import haploR
#' @import regioneR




# Functions ====================================================================

#' @title Set Difference of GRanges
#'
#' @description Set Difference of GRanges
#'
#' @param x list of GRanges
#' @param ignore_strand ignore strand if TRUE
#' @return GRanges, the set difference
setdiff_granges <- function(x, ignore_strand = FALSE) {
  if (length(x) > 1) {
    diffs <- list()
    for (i in 1:length(x)) {
      xi <- c(x[i], x[-i])
      diffs[[names(x)[i]]] <- Reduce(
        function(x, y) setdiff(x, y, ignore.strand = ignore_strand),
        xi
      )
    }
    diffs
  } else {
    x[[1]]
  }
}

#' @title SNPs per locus
#'
#' @description Collect SNPs at trait associated loci based on LD with lead SNP
#'
#' @param trait character, trait to consider
#' @param ld_threshold numeric, LD threshold for inclusion in a locus
#' @param population 1KGP population for LD reference
#' @return list of data frames giving haploR output per locus
snps_per_locus <- function(trait, ld_threshold = 0.5, population = "EUR") {
  if (population %in% unique(unlist(pop_list))) {
    lapply(
      mcols(subsetByTraits(ebicat37, tr = trait))[["SNPS"]],
      queryHaploreg,
      ldThresh = ld_threshold,
      ldPop = population
    )
  } else {
    lapply(
      mcols(subsetByTraits(ebicat37, tr = trait))[["SNPS"]],
      queryHaploreg,
      ldThresh = ld_threshold - 1,
      ldPop = pop_list[[population]]
    )
  }
}

#' @title SNPS to boundaries
#'
#' @description conver a vector of SNPS to boundaries
#'
#' @param snps vector of SNPs at a locus
#' @param population 1KGP population for LD reference
#' @return \describe{
#'   \item{chr}{chromosome of the locus}
#'   \item{start}{start position of the locus}
#'   \item{end}{end position of the locus}
#' }
snps_to_boundaries <- function(snps, ld_threshold = 0.5, population = "EUR") {
  if (population %in% unique(unlist(pop_list))) {
    positions <- start(
      ranges(ebicat37[intersect(getRsids(ebicat37), snps[["rsID"]])])
    )
  } else {
    ldmat <- LDlink.LDmatrix(snps[["rsID"]], population = population)
    lead_snp <- snps[["rsID"]][snps[["is_query_snp"]] == 1]
    rsids <- ldmat[["matrix.r2"]][
      ldmat[["matrix.r2"]][,lead_snp] > ld_threshold,
      RS_number
    ]
    positions <- start(
      ranges(ebicat37[intersect(getRsids(ebicat37), rsids)])
    )
  }
  start <- min(positions)
  end <- max(positions)
  list(chr = unique(snps[["chr"]]), start = start, end = end)
}

#' @title Risk Loci Enrichment Test
#'
#' @description test regions for enrichment with disease risk loci
#'
#' @param traits character, list of disease traits
#' @param regions GRanges, regions to test
#' @param population 1KGP population
#' @param ld_threshold LD threshold for definition of locus boundaries
#' @param ntimes number of permutations
#' @return \describe{
#'   \item{perm_test}{result of the permutation test}
#'   \item{loci}{GRanges object containing loci}
#'   \item{regions}{GRanges object containing regions}
#' }
#' @export
risk_loci_enrichment_test <- function(
    traits,
    regions,
    population = "EUR",
    ld_threshold = 0.5,
    ntimes = 50
) {
  loci_by_trait <- lapply(
    setNames(traits, traits),
    function(trait) {
      loci_boundaries <- lapply(
        snps_per_locus(
          trait,
          ld_threshold = ld_threshold,
          population = population
        ),
        snps_to_boundaries,
        population = population
      )
      chr <- as.character(
        sapply(loci_boundaries, function(locus) locus[["chr"]])
      )
      start <- as.integer(
        sapply(loci_boundaries, function(locus) locus[["start"]])
      )
      end <- as.integer(
        sapply(loci_boundaries, function(locus) locus[["end"]])
      )
      valid <- !is.na(start) & !is.na(end) & chr %in% valid_chromosomes
      if (startsWith(as.character(seqnames(regions))[1], "chr")) {
        chr <- paste("chr", chr, sep = "")
      }
      loci <- GRanges(
        seqnames = Rle(chr[valid]),
        ranges = IRanges(start = start[valid], end = end[valid])
      )
      union(loci, loci)
    }
  )
  loci_by_trait <- setdiff_granges(loci_by_trait, ignore_strand = TRUE)
  genome_and_mask <- getGenomeAndMask("hg19")
  lapply(
    setNames(traits, traits),
    function(trait) {
      list(
        perm_test = overlapPermTest(
          A = loci_by_trait[[trait]],
          B = regions,
          ntimes = ntimes,
          genome = genome_and_mask[["genome"]],
          mask = genome_and_mask[["mask"]]
        ),
        loci = loci_by_trait[[trait]],
        regions = regions
      )
    }
  )
}
