library(rentrez)
library(tidyverse)
library(biomaRt)

#====================================================================================================================#
#====================================================================================================================#

LiteratureRetrieval <- function(term, db = "pubmed", retmax = 250) {
  # Validate inputs
  if (!is.character(term) || length(term) != 1) {
    stop("The 'term' parameter must be a single string.")
  }
  if (!is.character(db) || length(db) != 1) {
    stop("The 'db' parameter must be a single string.")
  }
  if (!is.numeric(retmax) || retmax <= 0) {
    stop("The 'retmax' parameter must be a positive number.")
  }
  
  # Retrieve literature results
  literature_results <- entrez_search(db = db, term = term, retmax = retmax)
  
  # Check if any literature was found
  if (length(literature_results$ids) == 0) {
    stop("No literature found for the specified search terms.")
  }
  
  # Retrieve summaries for the literature results
  literature_summary <- entrez_summary(db = db, id = literature_results$ids)
  
  # Create a tibble to store the literature information
  literature_table <- tibble(
    Uid = extract_from_esummary(literature_summary, "uid", simplify = TRUE),
    Publication.Date = str_extract(extract_from_esummary(literature_summary, "pubdate", simplify = TRUE), "\\d{4}"),
    Journal = extract_from_esummary(literature_summary, "fulljournalname", simplify = TRUE),
    ISSN = extract_from_esummary(literature_summary, "issn", simplify = TRUE),
    ESSN = extract_from_esummary(literature_summary, "essn", simplify = TRUE),
    Title = extract_from_esummary(literature_summary, "title", simplify = TRUE),
    Electronic.Location = extract_from_esummary(literature_summary, "elocationid", simplify = TRUE),
    Authors = sapply(literature_summary, function(x) paste(x$authors$name, collapse = ", ")),
    Publication.Type = sapply(literature_summary, function(x) paste(x$pubtype, collapse = ", "))
  )
  
  return(literature_table)
}

#====================================================================================================================#
#====================================================================================================================#

GenomeRetrieval <- function(organism_name) {
  # Search for the taxonomy ID of the given organism
  taxonomy <- entrez_search(db = "taxonomy", term = organism_name, retmax = 1)
  taxonomy_id <- taxonomy$ids
  
  if (length(taxonomy$ids) == 0) {
    stop("No Taxonomy ID found for the specified organism.")
  }
  
  # Link the taxonomy ID to the genome database
  link_taxonomy_genome <- entrez_link(dbfrom = "taxonomy",
                                      id = taxonomy_id,
                                      db = "genome")
  
  genome_ids <- link_taxonomy_genome$links$taxonomy_genome
  
  # Link the genome IDs to the nucleotide database and retrieve web history
  link_genome_sequences <- entrez_link(dbfrom = "genome",
                                       id = genome_ids,
                                       db = "nuccore",
                                       cmd = "neighbor_history")
  
  web_history <- link_genome_sequences$web_histories
  
  # Fetch summaries of the genomic sequences
  genome_summary <- entrez_summary(db = "nuccore",
                                   web_history = web_history$genome_nuccore,
                                   retmax = 500)
  
  # Create a table with the relevant details
  genome_table <- tibble(
    Uid = extract_from_esummary(genome_summary, "uid", simplify = TRUE),
    Date = extract_from_esummary(genome_summary, "createdate", simplify = TRUE),
    Accession = extract_from_esummary(genome_summary, "caption", simplify = TRUE),
    Acc_version = extract_from_esummary(genome_summary, "accessionversion", simplify = TRUE),
    Title = extract_from_esummary(genome_summary, "title", simplify = TRUE),
    Seq.Length = extract_from_esummary(genome_summary, "slen", simplify = TRUE),
    Biomolecule = extract_from_esummary(genome_summary, "biomol", simplify = TRUE),
    Database = extract_from_esummary(genome_summary, "sourcedb", simplify = TRUE),
    Organism = extract_from_esummary(genome_summary, "organism", simplify = TRUE)
  )
  
  return(genome_table)
}

#====================================================================================================================#
#====================================================================================================================#

GenomeFetcher <- function(genome_ids, fasta_file) {
  
  # Fetch genome sequences from NCBI
  genome_sequence_fetch <- entrez_fetch(db = "nuccore", id = genome_ids, rettype = "fasta")
  
  # Write the fetched sequences to the specified fasta file
  writeLines(genome_sequence_fetch, fasta_file)
  
  return(paste("Sequences written into", fasta_file))
}

#====================================================================================================================#
#====================================================================================================================#

#GeneRetrieval <- function(gene_name,
#                          organism_name,
#                          retmax = 500) {
  
  # Construct the search term
#  search_term <- paste0(gene_name, "[TITL] gene AND ", organism_name, "[ORGN]")
  
  # Search for genes in the specified organism
#  genes <- entrez_search(db = "gene", term = search_term, retmax = retmax, use_history = TRUE)
  
  # Check if any genes were found
# if (genes$count == 0) {
#    stop("No genes found for the specified search term.")
#  }
  
  # Fetch summaries of the genes
#  genes_summary <- entrez_summary(db = "gene", web_history = genes$web_history)
  
  # Create a table with the relevant details
#  genes_table <- tibble(
#    Uid = extract_from_esummary(genes_summary, "uid", simplify = TRUE),
#    Chr.Acc.Number = sapply(genes_summary, function(x) x$genomicinfo$chraccver),
#    Name = extract_from_esummary(genes_summary, "name", simplify = TRUE),
#    Description = extract_from_esummary(genes_summary, "description", simplify = TRUE),
#    Organism = sapply(genes_summary, function(x) x$organism$scientificname),
#    Chromossome = sapply(genes_summary, function(x) x$genomicinfo$chrloc),
#    Start = sapply(genes_summary, function(x) x$genomicinfo$chrstart),
#    Stop = sapply(genes_summary, function(x) x$genomicinfo$chrstop),
#    Summary = extract_from_esummary(genes_summary, "summary", simplify = TRUE)
#  )
  
#  genes_table <- genes_table %>%
#    unnest(cols = c(Chr.Acc.Number, Chromossome, Start, Stop))
  
#  return(genes_table)
#}


GeneRetrieval <- function(gene_name,
                          organism_name,
                          retmax = 2000,
                          batch_size = 500) {
  
  # Construct the search term
  search_term <- paste0(gene_name, "[TITL] gene AND ", organism_name, "[ORGN]")
  
  # Search for genes in the specified organism
  genes <- entrez_search(db = "gene", term = search_term, retmax = retmax, use_history = TRUE)
  
  # Check if any genes were found
  if (genes$count == 0) {
    stop("No genes found for the specified search term.")
  }
  
  # Split the web history into batches
  web_hist <- genes$web_history
  batches <- seq(0, genes$count - 1, by = batch_size)
  
  # Function to fetch and process summaries for a batch
  fetch_summaries <- function(start) {
    summaries <- entrez_summary(db = "gene", web_history = web_hist, retstart = start, retmax = batch_size)
    
    tibble(
      Uid = extract_from_esummary(summaries, "uid", simplify = TRUE),
      Chr.Acc.Number = sapply(summaries, function(x) x$genomicinfo$chraccver),
      Name = extract_from_esummary(summaries, "name", simplify = TRUE),
      Description = extract_from_esummary(summaries, "description", simplify = TRUE),
      Organism = sapply(summaries, function(x) x$organism$scientificname),
      Chromosome = sapply(summaries, function(x) x$genomicinfo$chrloc),
      Start = sapply(summaries, function(x) x$genomicinfo$chrstart),
      Stop = sapply(summaries, function(x) x$genomicinfo$chrstop),
      Summary = extract_from_esummary(summaries, "summary", simplify = TRUE)
    )
  }
  
  # Fetch and process all batches
  all_genes <- map_dfr(batches, fetch_summaries)
  
  all_genes <- all_genes %>%
    unnest(cols = c(Chr.Acc.Number, Chromosome, Start, Stop))
  
  return(all_genes)
}


#====================================================================================================================#
#====================================================================================================================#

mRNARetrieval <- function(gene_uids) {
  # Link gene IDs to the nuccore database
  search_links <- entrez_link(dbfrom = "gene", id = gene_uids, db = "nuccore", cmd = "neighbor_history")
  
  # Check if any RNA sequences were found
  if (length(search_links$web_histories$gene_nuccore_refseqrna) == 0) {
    stop("No RNA sequences found for the specified gene UIDs.")
  }
  
  # Fetch summaries of the RNA sequences
  RNA_summary <- entrez_summary(db = "nuccore", web_history = search_links$web_histories$gene_nuccore_refseqrna)
  
  # Extract details from the summaries
  RNA_table <- tibble(
    Uid = extract_from_esummary(RNA_summary, "uid", simplify = TRUE),
    Accession = extract_from_esummary(RNA_summary, "caption", simplify = TRUE),
    Acc_version = extract_from_esummary(RNA_summary, "accessionversion", simplify = TRUE),
    Title = extract_from_esummary(RNA_summary, "title", simplify = TRUE),
    Seq_Length = extract_from_esummary(RNA_summary, "slen", simplify = TRUE),
    Biomolecule = extract_from_esummary(RNA_summary, "biomol", simplify = TRUE),
    Database = extract_from_esummary(RNA_summary, "sourcedb", simplify = TRUE),
    Organism = extract_from_esummary(RNA_summary, "organism", simplify = TRUE)
  )
  
  # Return the RNA details table
  return(RNA_table)
}

#====================================================================================================================#
#====================================================================================================================#

mRNAFetcher <- function(mRNA_ids, fasta_file) {
  
  # Fetch genome sequences from NCBI
  mRNA_sequence_fetch <- entrez_fetch(db = "nuccore",
                                        id = mRNA_ids,
                                        rettype = "fasta")
  
  # Write the fetched sequences to the specified fasta file
  writeLines(mRNA_sequence_fetch, fasta_file)
  
  return(paste("mRNA Sequences written into", fasta_file))
}

#====================================================================================================================#
#====================================================================================================================#

ProteinRetrieval <- function(gene_uids) {
  # Link RNA IDs to the protein database
  search_links <- entrez_link(dbfrom = "gene", id = gene_uids, db = "protein", cmd = "neighbor_history")
  
  # Check if any protein sequences were found
  if (length(search_links$web_histories$gene_protein_refseq) == 0) {
    stop("No protein sequences found for the specified RNA UIDs.")
  }
  
  # Fetch summaries of the protein sequences
  Protein_summary <- entrez_summary(db = "protein", web_history = search_links$web_histories$gene_protein_refseq)
  
  # Extract details from the summaries
  Protein_tibble <- tibble(
    Uid = extract_from_esummary(Protein_summary, "uid", simplify = TRUE),
    Accession = extract_from_esummary(Protein_summary, "caption", simplify = TRUE),
    Title = extract_from_esummary(Protein_summary, "title", simplify = TRUE),
    Seq_Length = extract_from_esummary(Protein_summary, "slen", simplify = TRUE),
    Biomolecule = extract_from_esummary(Protein_summary, "moltype", simplify = TRUE),
    Database = extract_from_esummary(Protein_summary, "sourcedb", simplify = TRUE),
    Organism = extract_from_esummary(Protein_summary, "organism", simplify = TRUE)
  )
  
  # Return the protein details table
  return(Protein_tibble)
}

#====================================================================================================================#
#====================================================================================================================#

ProteinFetcher <- function(protein_ids, fasta_file) {
  
  # Fetch genome sequences from NCBI
  Protein_sequence_fetch <- entrez_fetch(db = "protein",
                                      id = protein_ids,
                                      rettype = "fasta")
  
  # Write the fetched sequences to the specified fasta file
  writeLines(Protein_sequence_fetch, fasta_file)
  
  return(paste("Protein Sequences written into", fasta_file))
}

#====================================================================================================================#
#====================================================================================================================#

EnsemblRetrieval <- function(biomart_name,
                             dataset_name,
                             entrez_gene_ids,
                             host = "https://fungi.ensembl.org",
                             retrieve_gene = TRUE,
                             retrieve_mRNA = TRUE,
                             retrieve_protein = TRUE) {
  
  # Helper function to filter invalid characters
  filter_valid_sequences <- function(sequences, valid_chars) {
    gsub(paste0("[^", valid_chars, "]"), "", sequences)
  }
  
  # Connect to the specified Ensembl mart
  ensembl_mart <- useMart(host = host, biomart = biomart_name)
  
  # Select the dataset
  dataset <- useDataset(dataset_name, mart = ensembl_mart)
  
  # Get the Ensembl gene IDs
  results <- getBM(
    attributes = c("ensembl_gene_id", "entrezgene_id"),
    filters = "entrezgene_id",
    values = entrez_gene_ids,
    mart = dataset
  )
  
  if (nrow(results) == 0) {
    stop("No results found for the provided Entrez gene IDs.")
  }
  
  gene_tibble <- tibble()
  mRNA_tibble <- tibble()
  protein_tibble <- tibble()
  
  # Helper function to retrieve and process sequences
  retrieve_sequences <- function(attributes, output_file, column_order, description_field, sequence_type) {
    sequences <- getBM(
      attributes = attributes,
      filters = "ensembl_gene_id",
      values = results$ensembl_gene_id,
      mart = dataset
    )
    
    sequence_tibble <- as_tibble(sequences) %>%
      dplyr::select(all_of(column_order)) %>%
      arrange(ensembl_gene_id)
    
    # Ensure sequences are correctly formatted
    valid_sequences <- sequences[[description_field]]
    valid_sequences <- valid_sequences[!is.na(valid_sequences) & valid_sequences != ""]
    
    if (sequence_type == "DNA") {
      valid_sequences <- filter_valid_sequences(valid_sequences, "ACGTNacgtn")
      fasta_sequences <- DNAStringSet(setNames(valid_sequences,
                                               paste(sequences$ensembl_gene_id,
                                                     sep = "_")))
    } else if (sequence_type == "RNA") {
      valid_sequences <- filter_valid_sequences(valid_sequences, "ACGUNacgun")
      fasta_sequences <- RNAStringSet(setNames(valid_sequences,
                                               paste(sequences$ensembl_gene_id,
                                                     sep = "_")))
    } else if (sequence_type == "AA") {
      valid_sequences <- filter_valid_sequences(valid_sequences, "ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy*")
      fasta_sequences <- AAStringSet(setNames(valid_sequences,
                                              paste(sequences$ensembl_gene_id,
                                                    sep = "_")))
    } else {
      stop("Invalid sequence type specified. Must be 'DNA', 'RNA', or 'AA'.")
    }
    
    writeXStringSet(fasta_sequences, filepath = output_file)
    
    return(sequence_tibble)
  }
  
  # Retrieve gene sequences
  if (retrieve_gene) {
    gene_tibble <- retrieve_sequences(
      attributes = c("ensembl_gene_id",
                     "entrezgene_id",
                     "description",
                     "gene_exon_intron"),
      output_file = "ensembl_gene_sequence.fasta",
      column_order = c("ensembl_gene_id", "entrezgene_id", "description"),
      description_field = "gene_exon_intron",
      sequence_type = "DNA"
    )
    cat("Gene sequences have been successfully retrieved and saved to FASTA files.\n")
  }
  
  # Retrieve mRNA sequences
  if (retrieve_mRNA) {
    mRNA_tibble <- retrieve_sequences(
      attributes = c("ensembl_gene_id",
                     "entrezgene_id",
                     "ensembl_transcript_id",
                     "description",
                     "transcript_exon_intron"),
      output_file = "ensembl_mRNA_sequence.fasta",
      column_order = c("ensembl_gene_id", "entrezgene_id", "ensembl_transcript_id", "description"),
      description_field = "transcript_exon_intron",
      sequence_type = "RNA"
    )
    cat("mRNA sequences have been successfully retrieved and saved to FASTA files.\n")
  }
  
  # Retrieve protein sequences
  if (retrieve_protein) {
    protein_tibble <- retrieve_sequences(
      attributes = c("ensembl_gene_id",
                     "ensembl_peptide_id",
                     "description",
                     "peptide"),
      output_file = "ensembl_protein_sequence.fasta",
      column_order = c("ensembl_gene_id", "ensembl_peptide_id", "description"),
      description_field = "peptide",
      sequence_type = "AA"
    )
    cat("Protein sequences have been successfully retrieved and saved to FASTA files.\n")
  }
  
  if (!retrieve_gene && !retrieve_mRNA && !retrieve_protein) {
    stop("No sequence type specified. Please select at least one sequence type to retrieve.")
  }
  
  return(list(gene_sequences = gene_tibble, mRNA_sequences = mRNA_tibble, protein_sequences = protein_tibble))
}

#====================================================================================================================#
#====================================================================================================================#

EnsemblInfo <- function(kingdom = "vertebrates") {
  # Define the base URLs for different kingdoms
  ensembl_urls <- list(
    vertebrates = "https://www.ensembl.org",
    plants = "https://plants.ensembl.org",
    fungi = "https://fungi.ensembl.org",
    protists = "https://protists.ensembl.org",
    metazoa = "https://metazoa.ensembl.org"
  )
  
  # Check if the specified kingdom is valid
  if (!kingdom %in% names(ensembl_urls)) {
    stop("Invalid kingdom specified. Valid options are: vertebrates, plants, fungi, protists, metazoa.")
  }
  
  # Get the BioMart list for the specified kingdom
  ensembl_url <- ensembl_urls[[kingdom]]
  bioMarts <- listMarts(host = ensembl_url)
  
  # Function to get datasets for a chosen BioMart
  getDatasetsForBioMart <- function(bioMartName) {
    selectedMart <- useMart(biomart = bioMartName, host = ensembl_url)
    datasets <- listDatasets(selectedMart)
    return(datasets)
  }
  
  # Create a tibble to store the results
  results <- tibble()
  
  # Loop through each BioMart and get its datasets
  for (i in 1:nrow(bioMarts)) {
    bioMartName <- bioMarts[i, "biomart"]
    datasets <- getDatasetsForBioMart(bioMartName)
    
    # Create a tibble for the current BioMart
    bioMartTibble <- tibble(
      biomart = bioMartName,
      dataset = datasets$dataset,
      host = ensembl_url,
      description = datasets$description
    )
    
    # Append to the results tibble
    results <- bind_rows(results, bioMartTibble)
  }
  
  return(results)
}


