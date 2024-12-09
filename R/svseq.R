

#' Extract Fusion Sequence or Coordinates
#'
#' @param fusion fusion string (e.g BCR::chr22:23632645:chr9:133729466::ABL1). General form is `upstreamGene:upstreamChrom:upstreamPos-downstreamChrom:downstreamPos::downstreamGene`
#' @param ranges A set of genomic ranges representing the derivative chromosome
#' @param sv genomic breakpoint, start position, end position,
#' @return the genomic DNA fusion sequence from Gene1 start to Gene2 stop coordinates.
#' @export
#'
#' @examples
#' svseq(
#'   fusion = "BCR::chr22:23632645|chr9:133729466::ABL1",
#'   inversion = FALSE,  # Is the first gene inverted in its new position?
#'   ref = "hg38",
#'   return = "sequence"
#' )
svseq <- function(fusion, inversion=FALSE, ref = "hg38", return = c("sequence", "coords")){

  ref <- rlang::arg_match(ref)
  return <- rlang::arg_match(return)

  if(ref != "hg38"){
   cli::cli_abort("Only the hg38 reference genome is supported by svseq at the present")
  }

  # Reference Data (hg38)
  rlang::check_installed("TxDb.Hsapiens.UCSC.hg38.knownGene", reason = "to look up gene coordinates")
  rlang::check_installed("org.Hs.eg.db", reason = "to map gene symbol to entrez IDs")
  rlang::check_installed("BSgenome.Hsapiens.UCSC.hg38", reason = "extract the derived chromosoome sequence")
  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  db = org.Hs.eg.db::org.Hs.eg.db
  genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

  # Genomic ranges object that represents the derivative fusion genomic sequence.
  coords <- fusion_string_to_coords(fusion = fusion, inversion = inversion, txdb = txdb, db = db)

  if(return == "coords"){
    return(coords)
  }

  # Extract Sequence
  seq <- Biostrings::getSeq(genome, coords)

  # Collapse into one string
  singleseq <- unlist(seq)

  return(singleseq)
}

#' Export Genomic Coordinates as BED File
#'
#' This function exports a set of genomic coordinates to a BED file, suitable for downstream analysis or visualization.
#'
#' @param coords A `GenomicRanges` object representing genomic coordinates.
#' @param outfile A string specifying the output file path for the BED file.
#'
#' @return Invisible `NULL`. The BED file is written to the specified path.
#' @export
#'
#' @examples
#' export_bed(coords = example_sv, outfile = "example.bed")
export_bed <- function(coords, outfile){
  rlang::check_installed(pkg = "rtracklayer", reason = "export genomic ranges object as BED")
  rtracklayer::export.bed(coords, con = outfile)
}

#' Export DNA Sequence as FASTA File
#'
#' This function exports a DNA sequence to a FASTA file, enabling further bioinformatics analyses.
#'
#' @param sequence A string representing the DNA sequence to export.
#' @param outfile A string specifying the output file path for the FASTA file.
#' @param fasta_name An optional string to name the sequence in the FASTA file. Defaults to "fusion_sequence" if not provided.
#'
#' @return Invisible `NULL`. The FASTA file is written to the specified path.
#' @export
#'
#' @examples
#' export_fasta(sequence = "AGCTTGCA", outfile = "sequence.fasta", fasta_name = "example_sequence")
export_fasta <- function(sequence, outfile, fasta_name = NULL){
  rlang::check_installed(pkg = "Biostrings", reason = "export genomic ranges object as FASTA")
  seqset <- Biostrings::DNAStringSet(sequence)
  fasta_name <- if(is.null(fasta_name)) "fusion_sequence" else fasta_name
  names(seqset) <- fasta_name
  Biostrings::writeXStringSet(x = seqset, append = FALSE, filepath = outfile, format = "fasta")
}

fusion_string_to_coords <- function(fusion, inversion=FALSE, txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, db = org.Hs.eg.db::org.Hs.eg.db){

  components <- stringr::str_match(fusion, "^(.+)::(.+):(.+)\\|(.+):(.+)::(.+)$")
  assertions::assert(ncol(components) == 7, msg = "Unexpected fusion string. Should take the form {.emph GENE1::chrom1:position1|chrom2:position2::GENE2}")

  # Check if parsing worked well
  assertions::assert_no_missing(
  components,
  msg = "fusion string is not in the correct format. Should take the form {.emph GENE1::chrom1:position1|chrom2:position2::GENE2}")

  # Assign to variables
  gene1 <- components[2]
  chrom1 <- components[3]
  pos1 <- as.numeric(components[4])
  chrom2 <- components[5]
  pos2 <- as.numeric(components[6])
  gene2 <- components[7]

  # Print out info
  cli::cli_alert("Upstream Gene: {gene1}")
  cli::cli_alert("Upstream Gene SV: {chrom1}:{pos1}")
  cli::cli_alert("Upstream Gene Inverted: {inversion}")
  cli::cli_alert("Downstream Gene: {gene2}")
  cli::cli_alert("Downstream Gene SV: {chrom2}:{pos2}")


  # Fetch Start Coordinates
  gene1_coords <- lookup_gene_start_and_end(gene1, txdb=txdb, db=db)
  gene2_coords <- lookup_gene_start_and_end(gene2, txdb=txdb, db=db)

  # Adjust gene based on SV pos
  gene1_fixed <- set_end_position_in_read_dir(gene1_coords, pos1)
  gene2_fixed <- set_start_position_in_read_dir(gene2_coords, pos2)

  # Invert gene1 strand if its an inversion
  if(inversion){
    gene1_fixed <- BiocGenerics::invertStrand(gene1_fixed)
  }

  # Return Fixed Coordinates
  sv_coords <- c(gene1_fixed, gene2_fixed)
  names(sv_coords) <- c(paste0("Seq1 (", gene1,")"), paste0("Seq2 (", gene2,")"))
  return(sv_coords)
}

# Set the first position in read direction (start if + strand, end if -ve strand
set_start_position_in_read_dir <- function(range, newpos){
  if(as.character(GenomicRanges::strand(range)) == "+"){
    GenomicRanges::start(range) <- newpos
  }
  else {
    GenomicRanges::end(range) <- newpos
  }

  return(range)
}

set_end_position_in_read_dir <- function(range, newpos){
  if(as.character(GenomicRanges::strand(range)) == "+"){
    GenomicRanges::end(range) <- newpos
  }
  else{
    GenomicRanges::start(range) <- newpos
  }

  return(range)
}

lookup_gene_start_and_end <- function(gene, txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, db = org.Hs.eg.db::org.Hs.eg.db){

  # Grab Entrez ID
  entrez <- gene_symbol_to_entrez_id(gene, db=db)
  assertions::assert(length(entrez) == 1, msg = "More than 1 entrez ID corresponding to gene symbol {gene}")

  # Grab Gene Coords
  gene_coords <- GenomicFeatures::genes(txdb, filter = list(gene_id = entrez))

  # Assertions
  nrow(gene_coords)
  assertions::assert(length(gene_coords) == 1, msg = "More than 1 set of gene coordinates for gene symbol {gene}")

  # Return
  return(gene_coords)
}


gene_symbol_to_entrez_id <- function(gene, db = org.Hs.eg.db::org.Hs.eg.db){
  AnnotationDbi::select(db,
                        keys = gene,
                        keytype = "SYMBOL",
                        columns = c("ENTREZID"))[["ENTREZID"]]
}


example_sv <- GenomicRanges::GRanges(seqnames = c("chr11"), ranges = IRanges::IRanges(start = 63764705, end = 65663542))
