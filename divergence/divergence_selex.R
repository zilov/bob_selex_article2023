library('DiffLogo')

reverse_complement_matrix <- function(motif) {
  # Check if the matrix has row names "A", "C", "G", "T"
  if (!all(rownames(motif) %in% c("A", "C", "G", "T"))) {
    stop("Row names must be 'A', 'C', 'G', 'T'")
  }
  
  # Reverse the columns
  motif_rev <- motif[, rev(seq_len(ncol(motif)))]
  
  # Create an empty matrix to store the reverse complement
  motif_complement <- matrix(0, nrow=4, ncol=ncol(motif))
  
  # Take the complement while keeping original row names
  motif_complement[rownames(motif) == "A", ] <- motif_rev[rownames(motif) == "T", ]
  motif_complement[rownames(motif) == "T", ] <- motif_rev[rownames(motif) == "A", ]
  motif_complement[rownames(motif) == "C", ] <- motif_rev[rownames(motif) == "G", ]
  motif_complement[rownames(motif) == "G", ] <- motif_rev[rownames(motif) == "C", ]
  
  # Assign the row names
  rownames(motif_complement) <- rownames(motif)
  
  return(motif_complement)
}


meme_to_pwm <- function(file_path) {
  # Read the file into a data frame
  motif_probs <- read.csv(file_path, sep=" ", header=FALSE)
  
  # Set the column names to "A", "C", "G", "T"
  colnames(motif_probs) <- c("A", "C", "G", "T")
  
  # Transpose the matrix
  motif_probs_transposed <- t(as.matrix(motif_probs))
  
  # Set the column names to be the positions
  colnames(motif_probs_transposed) <- c(1:(ncol(motif_probs_transposed)))
  
  return(motif_probs_transposed)
}


o1_forv = read.csv("./data/O1_1-ATTTGCATAHNNNNNN_freqs.txt", sep=" ", header=FALSE)
colnames(o1_forv) = c("A", "C", "G", "T")
o1_forv = t(o1_forv)
colnames(o1_forv) = c(1:(ncol(o1_forv)))

o1_forv_logo = seqLogo::seqLogo(pwm = o1_forv)

o1_rev = read.csv("./data/O1_2-ATGCAAATBNVNNNNN_freqs_rev.txt", sep=" ", header=FALSE)
colnames(o1_rev) = c("A", "C", "G", "T")
o1_rev = t(o1_rev)
colnames(o1_rev) = c(1:(ncol(o1_rev)))
o1_rev_logo = seqLogo::seqLogo(pwm = o1_rev)
o1_rev_fixed = reverse_complement_matrix(o1_rev)
o1_rev_fixed_logo = seqLogo::seqLogo(pwm = o1_rev_fixed)

o1_motif = (o1_forv + o1_rev_fixed)/2
o1_motif_logo = seqLogo::seqLogo(pwm = o1_motif)


o1b_forv = meme_to_pwm("./data/O1B_1-ATTTGCATAHNNNNNN_freqs.txt")
o1b_forv_logo = seqLogo::seqLogo(pwm = o1b_forv)

o1b_rev = meme_to_pwm("./data/O1B_2-ATGCAAATNNNNNNNN_freqs_rev.txt")
o1b_rev_fixed = reverse_complement_matrix(o1b_rev)
o1b_rev_fixed_logo = seqLogo::seqLogo(pwm = o1b_rev_fixed)

o1b_motif = (o1b_forv + o1b_rev_fixed)/2
o1b_motif_logo = seqLogo::seqLogo(pwm = o1b_motif)

o2_forv = meme_to_pwm("./data/O2_1-ATTTGCATAHNNNNNN_freqs.txt")
o2_forv_logo = seqLogo::seqLogo(pwm = o2_forv)

o2_rev = meme_to_pwm("./data/O2_3-TATGCAAATBARNNNN_freqs_rev.txt")
o2_rev_fixed = reverse_complement_matrix(o2_rev)
o2_rev_fixed_logo = seqLogo::seqLogo(pwm = o2_rev_fixed)

o2_motif = (o2_forv + o2_rev_fixed)/2
o2_motif_logo = seqLogo::seqLogo(pwm = o2_motif)


o2b_forv = meme_to_pwm("./data/O2B_1-DATTTGCATAHNNNNN_freqs.txt")
o2b_forv_logo = seqLogo::seqLogo(pwm = o2b_forv)

o2b_rev = meme_to_pwm("./data/O2B_2-ATGCAAATNHNNNNNN_freqs_rev.txt")
o2b_rev_fixed = reverse_complement_matrix(o2b_rev)
o2b_rev_fixed_logo = seqLogo::seqLogo(pwm = o2b_rev_fixed)

o2b_motif = (o2b_forv + o2b_rev_fixed)/2
o2b_motif_logo = seqLogo::seqLogo(pwm = o2b_motif)

o1_motif = reverse_complement_matrix(o1_motif)
o1b_motif = reverse_complement_matrix(o1b_motif)

o2_motif = reverse_complement_matrix(o2_motif)
o2b_motif = reverse_complement_matrix(o2b_motif)

o1_motif_logo = seqLogo::seqLogo(pwm = o1_motif)
o1b_motif_logo = seqLogo::seqLogo(pwm = o1b_motif)
o2_motif_logo = seqLogo::seqLogo(pwm = o2_motif)
o2b_motif_logo = seqLogo::seqLogo(pwm = o2b_motif)

o1_o1b_diff = diffLogoFromPwm(o1_motif, o1b_motif)
o2_o2b_diff = diffLogoFromPwm(o2_motif, o2b_motif)

