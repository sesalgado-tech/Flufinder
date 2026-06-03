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