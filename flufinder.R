## Function 1: upload_fasta()

upload_fasta <- function(testfasta.txt) {
  #Opening seqinr library for handling FASTA files; make sure you have seqinr installed
  library(seqinr)
  
  #Reading the fasta file
  read.fasta(testfasta.txt, seqtype = "AA", as.string = TRUE,
             set.attributes = FALSE)
}

upload_fasta("testfasta.txt")

upload_fasta("testfasta.txt")$A


#2:
trypsinize <- function(proteins) {
  
  # Load stringr for string manipulation
  library(stringr)
  
  # "(?<=R|K)" is a regex "lookbehind" — it splits the string
  # at positions RIGHT AFTER an R or K, but keeps the R/K attached
  # to the left side. This is how Trypsin works: it cuts the
  # peptide bond on the C-terminal side of R and K.
  lapply(proteins, str_split_1, pattern="(?<=R|K)")
}


## Function 3: split_peptides()

split_peptides <- function(peptides) {
  # Opening stringr for simple string manipulation
  library(stringr)
  #Splitting peptides into individual amino acids using str_split; generates a list of lists of amino acids for each peptide
  lapply(peptides, str_split, pattern="")
}

peptides <- list(A=c("LVK","LHHIIFESMLK","DMQR","R","HR","VW"),
                 B=c("ADEFQGSMQK","IEACWQSYDVQF"),
                 C=c("MINEPFSWR","LEFHLSER","K","YDEIM"))

split_peptides(peptides)

split_peptides(peptides)$B


#4:
splitpeptides_to_masses <- function(aa) {
  
  # Mass table: monoisotopic mass (in Daltons) for each of the 
  # 20 standard amino acids
  aa_masses <- c(A=71.037, R=156.101, N=114.042, D=115.026, C=103.009, 
                 Q=128.058, E=129.042, G=57.021, H=137.058, I=113.084, 
                 L=113.084, K=128.094, M=131.040, F=147.068, P=97.052, 
                 S=87.032, T=101.047, W=186.079, Y=163.063, V=99.068)
  
  # Loop over each protein's list of split peptides
  # For each peptide (a character vector of amino acid letters),
  # use the letters to look up masses from aa_masses, then sum
  peptide_masses <- aa
  
  for(i in 1:length(aa)) {
    peptide_masses[[i]] <- lapply(aa[[i]], 
                                  function(x) sum(aa_masses[x]))
  }
  
  # Unlist the inner lists so each protein has a simple numeric vector
  lapply(peptide_masses, unlist)
}


## Function 5: count_matching_masses()

count_matching_masses <- function(protein_masses, sample, masses_list)
{
  #Virus masses is a list of masses for each protein so we use sapply to iterate over the list; sum (of TRUEs) is used to count the number of times a mass in the sample is found (%in%) among the masses of each of the proteins (virus_masses); note that masses are converted into strings (as.character) because %in% is not very reliable with numbers
  
  df <- as.data.frame(sapply(protein_masses, function (x)
    sum(as.character(sample) %in% as.character(x))))
  # Adding peptide_counts as the column name of the counts column
  
  names(df) <- "peptide_counts"
  return(df) 
}

masses_list <- list(
  A = c(340.246, 1348.728, 530.225),
  B = c(1121.476, 1469.624),
  C = c(1160.540, 1011.511, 651.255)
)

sample <- c(340.246, 530.225, 1348.728)

count_matching_masses(masses_list, sample)

# Create a new sample vector (randomly picked from masses_list)
new_sample <- c(340.246, 1121.476, 651.255)

# Run the function
count_matching_masses(masses_list, new_sample)


#6:
ggbarplot <- function(peptide_counts_table) {
  
  library(ggplot2)
  
  # geom_col draws bars with height = peptide_counts
  # rownames() extracts the flu strain names for the x-axis
  ggplot(peptide_counts_table) +
    aes(rownames(peptide_counts_table), peptide_counts) +
    geom_col(fill="steelblue", width=0.5) +
    theme_bw() +
    labs(x="Flu Strain", y="Peptide Counts")
}