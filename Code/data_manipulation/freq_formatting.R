#This file takes the individual data for B cell frequencies and builds a .txt file
#containing the averages and their difference from the S12 data

WD <- "~/Desktop/CobeyLab/bcell_id_pop_model/Data/transformations"

MLN_OUTPUT <- "~/Desktop/CobeyLab/bcell_id_pop_model/Data/transformations/mln_bcell_freq_means.txt"
SPLN_OUTPUT <- "~/Desktop/CobeyLab/bcell_id_pop_model/Data/transformations/spln_bcell_freq_means.txt"

setwd(WD)

b_freq <- read.csv('GC_B_Freq.csv',header=TRUE,stringsAsFactors = FALSE)

viruses <- c('PR8','D4 Sa B9','D4 Sb C7','D4 Ca1 1B1',
             'D4 Ca2 D6','D4 Cb A1','SEQ12')
days <- c(14,21,28)
day_names <- c('14 dpi','21 dpi','28 dpi')
col_names <- c('14 dpi','21 dpi','28 dpi','14 dpi - S12','21dpi - S12', '28dpi - S12')
locations <- c('GC MLN','GC SPLN')

for (location in locations){
  
  loc_mat <- matrix(nrow=length(viruses),ncol=length(days))
  
  for (i in 1:length(days)){
    day_vec <- c()
    for (j in 1:length(viruses)){
      spec_mean <- mean(b_freq[which(b_freq['virus']==viruses[j] & b_freq['dpi']==days[i] & b_freq['loc']==location),'freq'])
      day_vec[j] <- spec_mean
    }
    loc_mat[,i] <- day_vec
  }
  
  df <- as.data.frame(loc_mat,row.names=viruses)
  colnames(df) <- day_names
  
  diff_mat <- matrix(nrow=length(viruses),ncol=length(days))
  for (i in 1:length(days)){
    diff_vec <- c()
    for (j in 1:length(viruses)){
      spec_mean <- mean(b_freq[which(b_freq['virus']==viruses[j] & b_freq['dpi']==days[i] & b_freq['loc']==location),'freq'])
      diff_mean <- spec_mean - df['SEQ12',i]
      diff_vec[j] <- diff_mean
    }
    diff_mat[,i] <- diff_vec
  }
  
  comb_mat <- matrix(c(loc_mat,diff_mat), nrow=length(viruses), ncol=2*length(days))
  final_df <- as.data.frame(comb_mat,row.names=viruses)
  colnames(final_df) <- col_names
}

if (location == 'GC MLN') {
  write.table(final_df, file=MLN_OUTPUT)
} else if (location == 'GC SPLN')(
  write.table(final_df, file=SPLN_OUTPUT)
)





