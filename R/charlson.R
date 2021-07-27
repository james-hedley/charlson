# Project: R function to create a variable for Charlson Comorbidity Index (CCI),
#          based on data from the NSW Admitted Patient Data Collection (APDC)
# Adapted from: 'Charlson_code.SAS' created by Judy Simpson
# Changes: converting code from SAS to R
# Created by: James Hedley
# Date created: 22nd July 2021
# Last updated: 22nd July 2021


# Create an example dataset of diagnosis codes
library('tidyr')
library('dplyr')
library('comorbidity')


set.seed(12345)
example_icdcodes_data <- data.frame(id=rep(seq(1:200),5)) %>%
  as_tibble() %>%
  arrange(id) %>%
  group_by(id) %>%
  mutate(episode=seq_len(n())) %>%
  ungroup() %>%
  slice(rep(1:n(),each=51)) %>%
  group_by(id,episode) %>%
  mutate(diagnosis_code_seq=seq(0:50)-1) %>%
  ungroup() %>%
  arrange(id,episode,diagnosis_code_seq) %>%
  mutate(diagnosis_code=comorbidity::sample_diag(n())) %>%
  pivot_wider(names_from=diagnosis_code_seq,
              names_prefix='diagnosis_code',
              values_from=diagnosis_code) %>%
  rename(diagnosis_codeP=diagnosis_code0)


# Define a function that creates a flag variable based on the ICD codes specified ----
icdcode_flag <- function(data=data,
                         newvar=flag,
                         diagnosis_code_vars=paste0('diagnosis_code',c('P',seq(1:50))),
                         icdcodes=NULL,
                         regex=NULL) {

  # Inputs:
  # data - the dataset where a flag variable should be created
  # diagnosis_code_vars - a character vector of the existing variables containing diagnosis codes that should be checked
  # newvar - a new variable which will be coded as 1 if any of the diagnosis codes are flagged, otherwise 0
  # icdcodes - a character vector of ICD codes to be flagged. These must be an exact match in order to be flagged.
  #            e.g. 'I20.0' will not be flagged if icdcodes='I20'
  # regex - a character string containing a regular expression to be used to flag ICD codes


  # Load required libraries
  library('tidyr')
  library('dplyr')


  # Enquote variable names
  newvar <- enquo(newvar)


  # Only one of icdcodes or regex may be specified. If not the case, return an error message
  if (is.null(icdcodes) & is.null(regex)) {
    return(paste('ERROR: One of icdcodes or regex must be specified'))
  }

  if (!is.null(icdcodes) & !is.null(regex)) {
    return(paste('ERROR: Only one of icdcodes OR regex may be specified'))
  }


  # Create a flag variable based on specified list of ICD codes
  if (!is.null(icdcodes) & is.null(regex)) {
    data <- data %>%
      mutate(across(.cols=all_of(diagnosis_code_vars), # apply a function to each of these variables
                    .names='_tempflag_{.col}', # create a new variable, with the same name but with 'flag_' as a prefix
                    .fns= ~ 1*(match(.,icdcodes,nomatch=0)>0))) %>% # Check whether the diagnosis code matches a regular expression
      mutate(tempnewvar_=do.call(pmax,select(.,starts_with('_tempflag_')))) %>% # Set flag to the maximum value across all flags
      select(-starts_with('_tempflag_')) # Remove the separate flag variables for each diagnosis code
  }


  # Create a flag variable based on specified regular expression
  if (is.null(icdcodes) & !is.null(regex)) {
    data <- data %>%
      mutate(across(.cols=all_of(diagnosis_code_vars), # apply a function to each of these variables
                    .names='_tempflag_{.col}', # create a new variable, with the same name but with 'flag_' as a prefix
                    .fns= ~str_detect(.,regex))) %>% # Check whether the diagnosis code matches a regular expression
      mutate(tempnewvar_=do.call(pmax,select(.,starts_with('_tempflag_')))) %>% # Set flag to the maximum value across all flags
      select(-starts_with('_tempflag_')) # Remove the separate flag variables for each diagnosis code

  }


  # Rename flag variable using the name provided
  colnames(data)[colnames(data)=='tempnewvar_'] <- quo_name(newvar)


  # Return the dataset with the new flag variable
  return(data)

}





# Define a function that calculates the Charlson Comorbidity Index (CCI) ----
charlson <- function(data=data,
                     newvar=cci,
                     diagnosis_code_vars=paste0('diagnosis_code',c('P',seq(1:50))),
                     keep_flags=TRUE,
                     flag_vars=c('ami','chf','pvd','cva','dementia','pulmonary_disease',
                                 'ctd','peptic_ulcer','liver_disease','diabetes',
                                 'diabetes_complications','paraplegia','renal_disease',
                                 'cancer','metastatic_cancer','severe_liver_disease','hiv'),
                     categories=TRUE,
                     category_var=ccicat,
                     category_breaks=c(0,1,3,Inf),
                     category_labels=NULL,
                     regex_ami='(?i)^(I21|I22|I25\\.2)',
                     regex_chf='(?i)^(I50)',
                     regex_pvd='(?i)^(I71|I79\\.0|I73\\.9|R02|Z95\\.8|Z95\\.9)',
                     regex_cva=paste0(
                       '(?i)^(I60|I61|I62|I63|I64|I65|I66',
                       '|I67\\.0|I67\\.1|I67\\.2',
                       '|I67\\.4|I67\\.5|I67\\.6|I67\\.7|I67\\.8|I67\\.9',
                       '|I68\\.1|I68\\.2|I68\\.8|I69|G46',
                       '|G45\\.0|G45\\.1|G45\\.2|G45\\.4|G45\\.8|G45\\.9)'),
                     regex_dementia='(?i)^(F00|F01|F02|F05\\.1)',
                     regex_pulmonary_disease=paste0(
                       '(?i)^(J40|J41|J42|J43|J44|J45|J46|J47',
                       '|J60|J61|J62|J63|J64|J65|J66|J67)'),
                     regex_ctd=paste0(
                       '(?i)^(M05\\.0|M05\\.1|M05\\.2|M05\\.3',
                       '|M05\\.8|M05\\.9|M06\\.0',
                       '|M32|M34|M33\\.2|M06\\.3|M06\\.9|M35\\.3)'),
                     regex_peptic_ulcer='(?i)^(K25|K26|K27|K28)',
                     regex_liver_disease=paste0(
                       '(?i)^(K74\\.2|K74\\.3|K74\\.4|K74\\.5|K74\\.6',
                       '|K70\\.2|K70\\.3|K71\\.7|K73|K74\\.0)'),
                     regex_diabetes=paste0(
                       '(?i)^(E10\\.9|E11\\.9|E13\\.9|E14\\.9',
                       '|E10\\.1|E11\\.1|E13\\.1|E14\\.1',
                       '|E10\\.5|E11\\.5|E13\\.5|E14\\.5)'),
                     regex_diabetes_complications=paste0(
                       '(?i)^(E10\\.2|E11\\.2|E13\\.2|E14\\.2',
                       '|E10\\.3|E11\\.3|E13\\.3|E14\\.3',
                       '|E10\\.4|E11\\.4|E13\\.4|E14\\.4)'),
                     regex_paraplegia=paste0(
                       '(?i)^(G82\\.0|G82\\.1|G82\\.2',
                       '|G04\\.1|G81)'),
                     regex_renal_disease=paste0(
                       '(?i)^(N05\\.2|N05\\.3|N05\\.4|N05\\.5|N05\\.6',
                       '|N07\\.2|N07\\.3|N07\\.4',
                       '|N01|N03|N18|N19|N25)'),
                     regex_cancer=paste0(
                       '(?i)^(C40|C41|C43|C95|C96',
                       paste0('|C0',seq(0,9),collapse=''),
                       paste0('|C',seq(10,39),collapse=''),
                       paste0('|C',seq(45,49),collapse=''),
                       paste0('|C',seq(50,59),collapse=''),
                       paste0('|C',seq(60,69),collapse=''),
                       paste0('|C',seq(70,76),collapse=''),
                       paste0('|C',seq(80,85),collapse=''),
                       '|C88\\.3|C88\\.7|C88\\.9|C90\\.0|C90\\.1|C94\\.7',
                       '|C91|C92|C93',
                       '|C94\\.0|C94\\.1|C94\\.2|C94\\.3',
                       '|C94.51)'),
                     regex_metastatic_cancer='(?i)^(C77|C78|C79|C80)',
                     regex_severe_liver_disease='(?i)^(K72\\.9|K76\\.6|K76\\.7|K72\\.1)',
                     regex_hiv='(?i)^(B20|B21|B22|B23|B24)',
                     icdcodes_ami=NULL,
                     icdcodes_chf=NULL,
                     icdcodes_pvd=NULL,
                     icdcodes_cva=NULL,
                     icdcodes_dementia=NULL,
                     icdcodes_pulmonary_disease=NULL,
                     icdcodes_ctd=NULL,
                     icdcodes_peptic_ulcer=NULL,
                     icdcodes_liver_disease=NULL,
                     icdcodes_diabetes=NULL,
                     icdcodes_diabetes_complications=NULL,
                     icdcodes_paraplegia=NULL,
                     icdcodes_renal_disease=NULL,
                     icdcodes_cancer=NULL,
                     icdcodes_metastatic_cancer=NULL,
                     icdcodes_severe_liver_disease=NULL,
                     icdcodes_hiv=NULL,
                     weight_ami=1,
                     weight_chf=1,
                     weight_pvd=1,
                     weight_cva=1,
                     weight_dementia=1,
                     weight_pulmonary_disease=1,
                     weight_ctd=1,
                     weight_peptic_ulcer=1,
                     weight_liver_disease=1,
                     weight_diabetes=1,
                     weight_diabetes_complications=2,
                     weight_paraplegia=2,
                     weight_renal_disease=2,
                     weight_cancer=2,
                     weight_metastatic_cancer=6,
                     weight_severe_liver_disease=3,
                     weight_hiv=6) {

  # Inputs:
  # data - the dataset where the Charlson Comorbidity Index should be calculated
  # diagnosis_code_vars - a character vector of the existing variables containing diagnosis codes that should be checked
  # newvar - a new variable which will contain the Charlson Comorbidity Index
  # keep_flags - a logical indicating whether the separate flags for each comorbidity should be kept in the final dataset
  # flag_vars - a character vector of length 17 containing variable names that will be used to create speparate flags for each comorbidity
  # categories - a logical indicating whether the Charlson Comorbidtiy Index should be categorised
  # category_var - a new variable containing the categorised Charlson Comorbidity Index
  # category_breaks - a numerical vector with the break points that should be used to categorise the Charlson Comorbidity Index
  #
  # regex_ami - a character string containing a regular expression to flag ICD codes for Acute myocardial infarction
  # regex_chf - a character string containing a regular expression to flag ICD codes for Congestive heart failure
  # regex_pvd - a character string containing a regular expression to flag ICD codes for Peripheral vascular disease
  # regex_cva - a character string containing a regular expression to flag ICD codes for Cerebral vascular disease
  # regex_dementia - a character string containing a regular expression to flag ICD codes for Dementia
  # regex_pulmonary_disease - a character string containing a regular expression to flag ICD codes for Pulmonary disease
  # regex_ctd - a character string containing a regular expression to flag ICD codes for Connective tissue disorder
  # regex_peptic_ulcer - a character string containing a regular expression to flag ICD codes for Peptic ulcer
  # regex_liver_disease - a character string containing a regular expression to flag ICD codes for Liver disease
  # regex_diabetes - a character string containing a regular expression to flag ICD codes for Diabetes
  # regex_diabetes_complications - a character string containing a regular expression to flag ICD codes for Diabetes complications
  # regex_paraplegia - a character string containing a regular expression to flag ICD codes for Paraplegia
  # regex_renal_disease - a character string containing a regular expression to flag ICD codes for Renal disease
  # regex_cancer - a character string containing a regular expression to flag ICD codes for Cancer
  # regex_metastatic_cancer - a character string containing a regular expression to flag ICD codes for Metastatic cancer
  # regex_severe_liver_disease - a character string containing a regular expression to flag ICD codes for Severe liver disease
  # regex_hiv - a character string containing a regular expression to flag ICD codes for HIV
  #
  # icdcodes_ami - a character vectore containing ICD codes to be matched exactly for Acute myocardial infarction. Overrides regular expression
  # icdcodes_chf - a character vectore containing ICD codes to be matched exactly for Congestive heart failure. Overrides regular expression
  # icdcodes_pvd - a character vectore containing ICD codes to be matched exactly for Peripheral vascular disease. Overrides regular expression
  # icdcodes_cva - a character vectore containing ICD codes to be matched exactly for Cerebral vascular disease. Overrides regular expression
  # icdcodes_dementia - a character vectore containing ICD codes to be matched exactly for Dementia. Overrides regular expression
  # icdcodes_pulmonary_disease - a character vectore containing ICD codes to be matched exactly for Pulmonary disease. Overrides regular expression
  # icdcodes_ctd - a character vectore containing ICD codes to be matched exactly for Connective tissue disorder. Overrides regular expression
  # icdcodes_peptic_ulcer - a character vectore containing ICD codes to be matched exactly for Peptic ulcer. Overrides regular expression
  # icdcodes_liver_disease - a character vectore containing ICD codes to be matched exactly for Liver disease. Overrides regular expression
  # icdcodes_diabetes - a character vectore containing ICD codes to be matched exactly for Diabetes. Overrides regular expression
  # icdcodes_diabetes_complications - a character vectore containing ICD codes to be matched exactly for Diabetes complications. Overrides regular expression
  # icdcodes_paraplegia - a character vectore containing ICD codes to be matched exactly for Paraplegia. Overrides regular expression
  # icdcodes_renal_disease - a character vectore containing ICD codes to be matched exactly for Renal disease. Overrides regular expression
  # icdcodes_cancer - a character vectore containing ICD codes to be matched exactly for Cancer. Overrides regular expression
  # icdcodes_metastatic_cancer - a character vectore containing ICD codes to be matched exactly for Metastatic cancer. Overrides regular expression
  # icdcodes_severe_liver_disease - a character vectore containing ICD codes to be matched exactly for Severe liver disease. Overrides regular expression
  # icdcodes_hiv - a character vectore containing ICD codes to be matched exactly for HIV. Overrides regular expression
  #
  # weight_ami - the weight that should be applied to the comorbidity flag
  # weight_chf - the weight that should be applied to the comorbidity flag
  # weight_pvd - the weight that should be applied to the comorbidity flag
  # weight_cva - the weight that should be applied to the comorbidity flag
  # weight_dementia - the weight that should be applied to the comorbidity flag
  # weight_pulmonary_disease - the weight that should be applied to the comorbidity flag
  # weight_ctd - the weight that should be applied to the comorbidity flag
  # weight_peptic_ulcer - the weight that should be applied to the comorbidity flag
  # weight_liver_disease - the weight that should be applied to the comorbidity flag
  # weight_diabetes - the weight that should be applied to the comorbidity flag
  # weight_diabetes_complications - the weight that should be applied to the comorbidity flag
  # weight_paraplegia - the weight that should be applied to the comorbidity flag
  # weight_renal_disease - the weight that should be applied to the comorbidity flag
  # weight_cancer - the weight that should be applied to the comorbidity flag
  # weight_metastatic_cancer - the weight that should be applied to the comorbidity flag
  # weight_severe_liver_disease - the weight that should be applied to the comorbidity flag
  # weight_hiv - the weight that should be applied to the comorbidity flag


  # Load required libraries
  library('tidyr')
  library('dplyr')


  # Enquote variables
  newvar <- enquo(newvar)
  category_var <- enquo(category_var)


  # If ICD codes are specified, replace regex with NULL
  if (!is.null(icdcodes_ami)) {regex_ami <- NULL}
  if (!is.null(icdcodes_chf)) {regex_chf <- NULL}
  if (!is.null(icdcodes_pvd)) {regex_pvd <- NULL}
  if (!is.null(icdcodes_cva)) {regex_cva <- NULL}
  if (!is.null(icdcodes_dementia)) {regex_dementia <- NULL}
  if (!is.null(icdcodes_pulmonary_disease)) {regex_pulmonary_disease <- NULL}
  if (!is.null(icdcodes_ctd)) {regex_ctd <- NULL}
  if (!is.null(icdcodes_peptic_ulcer)) {regex_peptic_ulcer <- NULL}
  if (!is.null(icdcodes_liver_disease)) {regex_liver_disease <- NULL}
  if (!is.null(icdcodes_diabetes)) {regex_diabetes <- NULL}
  if (!is.null(icdcodes_diabetes_complications)) {regex_diabetes_complications <- NULL}
  if (!is.null(icdcodes_paraplegia)) {regex_paraplegia <- NULL}
  if (!is.null(icdcodes_renal_disease)) {regex_renal_disease <- NULL}
  if (!is.null(icdcodes_cancer)) {regex_cancer <- NULL}
  if (!is.null(icdcodes_metastatic_cancer)) {regex_metastatic_cancer <- NULL}
  if (!is.null(icdcodes_severe_liver_disease)) {regex_severe_liver_disease <- NULL}
  if (!is.null(icdcodes_hiv)) {regex_hiv <- NULL}


  # Create a flag for Acute myocardial infarction
  data <- data %>%
    icdcode_flag(ami_,
                 diagnosis_code_vars=diagnosis_code_vars,
                 regex=regex_ami,
                 icdcodes=icdcodes_ami)

  # Create a flag for Congestive heart failure
  data <- data %>%
    icdcode_flag(chf_,
                 diagnosis_code_vars=diagnosis_code_vars,
                 regex=regex_chf,
                 icdcodes=icdcodes_chf)

  # Create a flag for Peripheral vascular disease
  data <- data %>%
    icdcode_flag(pvd_,
                 diagnosis_code_vars=diagnosis_code_vars,
                 regex=regex_pvd,
                 icdcodes=icdcodes_pvd)

  # Create a flag for Cerebral vascular accident
  data <- data %>%
    icdcode_flag(cva_,
                 diagnosis_code_vars=diagnosis_code_vars,
                 regex=regex_cva,
                 icdcodes=icdcodes_cva)

  # Create a flag for Dementia
  data <- data %>%
    icdcode_flag(dementia_,
                 diagnosis_code_vars=diagnosis_code_vars,
                 regex=regex_dementia,
                 icdcodes=icdcodes_dementia)

  # Create a flag for Pulmonary disease
  data <- data %>%
    icdcode_flag(pulmonary_disease_,
                 diagnosis_code_vars=diagnosis_code_vars,
                 regex=regex_pulmonary_disease,
                 icdcodes=icdcodes_pulmonary_disease)

  # Create a flag for Connective tissue disorder
  data <- data %>%
    icdcode_flag(ctd_,
                 diagnosis_code_vars=diagnosis_code_vars,
                 regex=regex_ctd)

  # Create a flag for Peptic ulcer
  data <- data %>%
    icdcode_flag(peptic_ulcer_,
                 diagnosis_code_vars=diagnosis_code_vars,
                 regex=regex_peptic_ulcer,
                 icdcodes=icdcodes_peptic_ulcer)

  # Create a flag for liver disease
  data <- data %>%
    icdcode_flag(liver_disease_,
                 diagnosis_code_vars=diagnosis_code_vars,
                 regex=regex_liver_disease,
                 icdcodes=icdcodes_liver_disease)

  # Create a flag for Diabetes
  data <- data %>%
    icdcode_flag(diabetes_,
                 diagnosis_code_vars=diagnosis_code_vars,
                 regex=regex_diabetes,
                 icdcodes=icdcodes_diabetes)

  # Create a flag for Diabetes complications
  data <- data %>%
    icdcode_flag(diabetes_complications_,
                 diagnosis_code_vars=diagnosis_code_vars,
                 regex=regex_diabetes_complications,
                 icdcodes=icdcodes_diabetes_complications)

  # Create a flag for Paraplegia
  data <- data %>%
    icdcode_flag(paraplegia_,
                 diagnosis_code_vars=diagnosis_code_vars,
                 regex=regex_paraplegia,
                 icdcodes=icdcodes_paraplegia)

  # Create a flag for Renal disease
  data <- data %>%
    icdcode_flag(renal_disease_,
                 diagnosis_code_vars=diagnosis_code_vars,
                 regex=regex_renal_disease,
                 icdcodes=icdcodes_renal_disease)

  # Create a flag for Cancer
  data <- data %>%
    icdcode_flag(cancer_,
                 diagnosis_code_vars=diagnosis_code_vars,
                 regex=regex_cancer,
                 icdcodes=icdcodes_cancer)

  # Create a flag for Metastatic cancer
  data <- data %>%
    icdcode_flag(metastatic_cancer_,
                 diagnosis_code_vars=diagnosis_code_vars,
                 regex=regex_metastatic_cancer,
                 icdcodes=icdcodes_metastatic_cancer)

  # Create a flag for Severe liver disease
  data <- data %>%
    icdcode_flag(severe_liver_disease_,
                 diagnosis_code_vars=diagnosis_code_vars,
                 regex=regex_severe_liver_disease,
                 icdcodes=icdcodes_severe_liver_disease)

  # Create a flag for HIV
  data <- data %>%
    icdcode_flag(hiv_,
                 diagnosis_code_vars=diagnosis_code_vars,
                 regex=regex_hiv,
                 icdcodes=icdcodes_hiv)


  # Apply weights to each flag
  data <- data %>%
    mutate(
      weighted_ami_=ami_*weight_ami,
      weighted_chf_=chf_*weight_chf,
      weighted_pvd_=pvd_*weight_pvd,
      weighted_cva_=cva_*weight_cva,
      weighted_dementia_=dementia_*weight_dementia,
      weighted_pulmonary_disease_=pulmonary_disease_*weight_pulmonary_disease,
      weighted_ctd_=ctd_*weight_ctd,
      weighted_peptic_ulcer_=peptic_ulcer_*weight_peptic_ulcer,
      weighted_liver_disease_=liver_disease_*weight_liver_disease,
      weighted_diabetes_=diabetes_*weight_diabetes,
      weighted_diabetes_complications_=diabetes_complications_*weight_diabetes_complications,
      weighted_paraplegia_=paraplegia_*weight_paraplegia,
      weighted_renal_disease_=renal_disease_*weight_renal_disease,
      weighted_cancer_=cancer_*weight_cancer,
      weighted_metastatic_cancer_=metastatic_cancer_*weight_metastatic_cancer,
      weighted_severe_liver_disease_=severe_liver_disease_*weight_severe_liver_disease,
      weighted_hiv_=hiv_*weight_hiv)



  # Create the Charlson Comorbidity Index
  data <- data %>%
    mutate(cci_=
      weighted_ami_+
      weighted_chf_+
      weighted_pvd_+
      weighted_cva_+
      weighted_dementia_+
      weighted_pulmonary_disease_+
      weighted_ctd_+
      weighted_peptic_ulcer_+
      weighted_liver_disease_+
      weighted_diabetes_+
      weighted_diabetes_complications_+
      weighted_paraplegia_+
      weighted_renal_disease_+
      weighted_cancer_+
      weighted_metastatic_cancer_+
      weighted_severe_liver_disease_+
      weighted_hiv_) %>%
    select(-weighted_ami_,
           -weighted_chf_,
           -weighted_pvd_,
           -weighted_cva_,
           -weighted_dementia_,
           -weighted_pulmonary_disease_,
           -weighted_ctd_,
           -weighted_peptic_ulcer_,
           -weighted_liver_disease_,
           -weighted_diabetes_,
           -weighted_diabetes_complications_,
           -weighted_paraplegia_,
           -weighted_renal_disease_,
           -weighted_cancer_,
           -weighted_metastatic_cancer_,
           -weighted_severe_liver_disease_,
           -weighted_hiv_)


  # Create categories of Charlson Comorbidity Index
  if (categories==TRUE) {
    data <- data %>%
      mutate(
        ccicat_=cut(cci_,
                    breaks=category_breaks,
                    right=FALSE))
  }

  # Apply category labels to Charlson Comorbidity categories
  if (categories==TRUE & !is.null(category_labels)) {
    data <- data %>%
      mutate(
        ccicat_=factor(ccicat_,
                       labels=category_labels))
  }


  # Remove flag variables if option to keep them is not selected
  if (keep_flags==FALSE) {
    data <- data %>%
      select(-ami_,
             -chf_,
             -pvd_,
             -cva_,
             -dementia_,
             -pulmonary_disease_,
             -ctd_,
             -peptic_ulcer_,
             -liver_disease_,
             -diabetes_,
             -diabetes_complications_,
             -paraplegia_,
             -renal_disease_,
             -cancer_,
             -metastatic_cancer_,
             -severe_liver_disease_,
             -hiv_)
  }


  # Rename variables using the names provided
  colnames(data)[colnames(data)=='cci_'] <- quo_name(newvar)

  if (categories==TRUE) {
    colnames(data)[colnames(data)=='ccicat_'] <- quo_name(category_var)
  }

  if (keep_flags==TRUE) {
    tempflagvars <- c('ami_','chf_','pvd_','cva_','dementia_','pulmonary_disease_',
                      'ctd_','peptic_ulcer_','liver_disease_','diabetes_',
                      'diabetes_complications_','paraplegia_','renal_disease_',
                      'cancer_','metastatic_cancer_','severe_liver_disease_','hiv_')
    colnames(data)[colnames(data) %in% tempflagvars] <- flag_vars
  }


  # Return the dataset with the newly created variables
  return(data)

}
