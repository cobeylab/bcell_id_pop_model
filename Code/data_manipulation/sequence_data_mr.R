#This file produces two .txt files, WIDE_OUTPUT AND LONG_OUTPUT
#which are text versions of wide and long form dataframes

#Neither are particularly useful for human reading, but I threw them together anticipating
#needing machine readable versions of that data sometime

WD <- "~/Desktop/CobeyLab/bcell_id_pop_model/Data/reidms"

WIDE_OUTPUT <- "~/Desktop/CobeyLab/bcell_id_pop_model/Data/transformations/wide_form_sequence_data.txt"
LONG_OUTPUT <- "~/Desktop/CobeyLab/bcell_id_pop_model/Data/transformations/long_form_sequence_data.txt"


setwd(WD)

library('reshape2')

get_sequences <- function(filename){
  raw_lines <- read.table(file = filename, header = FALSE,
                          stringsAsFactors = FALSE)
  
  sequence_vec <- raw_lines[which(as.numeric(row.names(raw_lines)) %% 2 == 0),1]
  
  return(sequence_vec)
}

get_names <- function(filename){
  raw_lines <- read.table(file = filename, header = FALSE,
                          stringsAsFactors = FALSE)
  raw_names <- raw_lines[which(as.numeric(row.names(raw_lines)) %% 2 != 0),1]
  split_names <- unlist(strsplit(raw_names, split=c('>'), fixed=TRUE))
  names <- split_names[which(split_names != '')]
  names <- unlist(strsplit(names,split='-Hchain', fixed=TRUE))
}

names <- get_names("Ab sequences aminoacids H.txt")
aa_h_chain <- get_sequences("Ab sequences aminoacids H.txt")
aa_k_chain <- get_sequences("Ab sequences aminoacids K.txt")
nt_h_sequence <- get_sequences("Ab sequences nucleotide H.txt")
nt_k_sequence <- get_sequences("Ab sequences nucleotide K.txt")

sequence_df_wide <- data.frame(names, aa_h_chain, aa_k_chain,
                               nt_h_sequence, nt_k_sequence, 
                               stringsAsFactors = FALSE)

sequence_df_long <- m_df <- melt(sequence_df_wide, id=c('names')) #melts into long form data.frame

write.table(sequence_df_wide, file=WIDE_OUTPUT)
write.table(sequence_df_long, file=LONG_OUTPUT)