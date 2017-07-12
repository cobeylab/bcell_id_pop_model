# Modeling B Cell Immunodominance and Population Dynamics
This repository contains all of the work related to B-Cell Immunodominance and population dynamical models

**Literature**:
* Smith (1999)
  - from PNAS
* Skowronski and Smith (2017)
  - from J. Infectious Diseases
* Zarnitsyna and Antia (2015)
  - from Phil. Trans.
* Zarnitsyna and Antia (2016)
  - from PLOS Pathogens
* Ndifon (2015)
  - from J. of the Royal Society Interface
* Wang and Chakraborty (2015)
  - from Cell
* Angeletti and Yewdell (2017)
  - from Nature
  - Note: This paper is useful for its *data*

**Code**:
* Basic Implementation
  * Z+A_3_epitope.R
     - Zarnitsyna et al. 2016 3 epitope model
  * log_growth_test.R
    - Code for scalable logistic growth with competition
  * log_growth_fixed_R
    - Code to examine how fixed rs for each population impacts dynamics
  * comp_matrix.csv & fixed_r.csv
    - .csv files called in log_growth_test and log_growth_fixed_R respectively
* Data Manipulation
  * sequence_data_mr.R
    - Code which produces data files aggregating the sequence data for each Ab
    - Produced long and wide_form sequence data
  * freq_formatting.R
    - takes GC_B_Freq.csv and produces the means and differences from S12
    - Produced mln and spln_freq_means.txt

**Data**
* reidms.zip
  * Just a zip file that makes the directory reidms
* reidms (The original data as sent by D.A.)
  * numbers.xlsx
  * Ab sequences --- .txt
* transformations
  * GC_B_Freq.csv
  * long_form_sequence_data.txt
  * wide_form_sequence_data.txt
  * mln_bcell_freq_means.txt
  * spln_bcell_freq_means.txt
