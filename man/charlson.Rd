\name{charlson}
\alias{charlson}

\title{charlson}

\description{
Calculate the Charlson Comorbidity Index (CCI) based on ICD-10 diagnosis codes. Also includes a function to create a binary flag variable (0 or 1) for whether any speicfied ICD codes are present

}
\usage{
icdcode_flag(data = data,
             newvar = flag,
             diagnosis_code_vars = paste0('diagnosis_code',c('P',seq(1:50))),
             icdcodes=NULL,
             regex=NULL)


charlson(data = data,
         newvar = cci,
         diagnosis_code_vars = paste0('diagnosis_code',c('P',seq(1:50))),
         keep_flags = TRUE,
         flag_vars = c('ami','chf','pvd','cva','dementia','pulmonary_disease',
                     'ctd','peptic_ulcer','liver_disease','diabetes',
                     'diabetes_complications','paraplegia','renal_disease',
                     'cancer','metastatic_cancer','severe_liver_disease','hiv'),
         categories = TRUE,
         category_var = ccicat,
         category_breaks = c(0,1,3,Inf),
         category_labels = NULL,
         regex_ami = '^(I21|I22|I25\\.2)',
         regex_chf = '^(I50)',
         regex_pvd = '^(I71|I79\\.0|I73\\.9|R02|Z95\\.8|Z95\\.9)',
         regex_cva = paste0(
           '^(I60|I61|I62|I63|I64|I65|I66',
           '|I67\\.0|I67\\.1|I67\\.2',
           '|I67\\.4|I67\\.5|I67\\.6|I67\\.7|I67\\.8|I67\\.9',
           '|I68\\.1|I68\\.2|I68\\.8|I69|G46',
           '|G45\\.0|G45\\.1|G45\\.2|G45\\.4|G45\\.8|G45\\.9)'),
         regex_dementia = '^(F00|F01|F02|F05\\.1)',
         regex_pulmonary_disease = paste0(
           '^(J40|J41|J42|J43|J44|J45|J46|J47',
           '|J60|J61|J62|J63|J64|J65|J66|J67)'),
         regex_ctd = paste0(
           '^(M05\\.0|M05\\.1|M05\\.2|M05\\.3',
           '|M05\\.8|M05\\.9|M06\\.0',
           '|M32|M34|M33\\.2|M06\\.3|M06\\.9|M35\\.3)'),
         regex_peptic_ulcer = '^(K25|K26|K27|K28)',
         regex_liver_disease = paste0(
           '^(K74\\.2|K74\\.3|K74\\.4|K74\\.5|K74\\.6',
           '|K70\\.2|K70\\.3|K71\\.7|K73|K74\\.0)'),
         regex_diabetes = paste0(
           '^(E10\\.9|E11\\.9|E13\\.9|E14\\.9',
           '|E10\\.1|E11\\.1|E13\\.1|E14\\.1',
           '|E10\\.5|E11\\.5|E13\\.5|E14\\.5)'),
         regex_diabetes_complications = paste0(
           '^(E10\\.2|E11\\.2|E13\\.2|E14\\.2',
           '|E10\\.3|E11\\.3|E13\\.3|E14\\.3',
           '|E10\\.4|E11\\.4|E13\\.4|E14\\.4)'),
         regex_paraplegia = paste0(
           '^(G82\\.0|G82\\.1|G82\\.2',
           '|G04\\.1|G81)'),
         regex_renal_disease = paste0(
           '^(N05\\.2|N05\\.3|N05\\.4|N05\\.5|N05\\.6',
           '|N07\\.2|N07\\.3|N07\\.4',
           '|N01|N03|N18|N19|N25)'),
         regex_cancer = paste0(
           '^(C40|C41|C43|C95|C96',
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
         regex_metastatic_cancer = '^(C77|C78|C79|C80)',
         regex_severe_liver_disease = '^(K72\\.9|K76\\.6|K76\\.7|K72\\.1)',
         regex_hiv = '^(B20|B21|B22|B23|B24)',
         icdcodes_ami = NULL,
         icdcodes_chf = NULL,
         icdcodes_pvd = NULL,
         icdcodes_cva = NULL,
         icdcodes_dementia = NULL,
         icdcodes_pulmonary_disease = NULL,
         icdcodes_ctd = NULL,
         icdcodes_peptic_ulcer = NULL,
         icdcodes_liver_disease = NULL,
         icdcodes_diabetes = NULL,
         icdcodes_diabetes_complications = NULL,
         icdcodes_paraplegia = NULL,
         icdcodes_renal_disease = NULL,
         icdcodes_cancer = NULL,
         icdcodes_metastatic_cancer = NULL,
         icdcodes_severe_liver_disease = NULL,
         icdcodes_hiv = NULL,
         weight_ami = 1,
         weight_chf = 1,
         weight_pvd = 1,
         weight_cva = 1,
         weight_dementia = 1,
         weight_pulmonary_disease = 1,
         weight_ctd = 1,
         weight_peptic_ulcer = 1,
         weight_liver_disease = 1,
         weight_diabetes = 1,
         weight_diabetes_complications = 2,
         weight_paraplegia = 2,
         weight_renal_disease = 2,
         weight_cancer = 2,
         weight_metastatic_cancer = 6,
         weight_severe_liver_disease = 3,
         weight_hiv = 6)


}

\arguments{
  \item{data}{
A tidy dataset where the new variables should be added
}
  \item{newvar}{
icdcode_flag(): The name of a new variable to be created containing a binary (0 or 1) flag for whether the specified ICD codes

charlson(): The name of a new variable to be created containing the Charlson Comorbidity Index

}
  \item{diagnosis_code_vars}{
A character vector of variable names that contain ICD-10 diagnosis codes
}
  \item{icdcodes}{
A character vector of ICD codes to be matched exactly. Exactly one of icdcodes or regex must be specified
}
  \item{regex}{
A character string containing a regular expression to match ICD codes. Exactly one of icdcodes or regex must be specified
}
  \item{keep_flags}{
A logical indicating whether the 17 individual comorbidity flag variables (0 or 1) should be kept in the dataset
}
  \item{flag_vars}{
A character vector with the names that should be given to the 17 individual comorbidity flag variables (0 or 1)
}
  \item{categories}{
A logical indicating whether a factor variable should be created to place Charlson Comrobdiity Index into groups
}
  \item{category_var}{
The name of the new factor variable for the groups of Charlson Comorbidity Index
}
  \item{category_breaks}{
A numeric vector of length >= 2, containing the break points to be used for determining the groups of Charlson Comorbidity Index

For example: c(0,1,3,Inf) creates 3 groups, (>=0 & <1), (>=1 & <3), (>=3)
}
  \item{category_labels}{
The character vector of the labels to be applied to each group of Charlson Comorbidity Index

For example, with category_breaks = c(0,1,3,Inf) the labels could be '0', '1-2', '3+'
}
  \item{regex_ami}{
A character string containing a regular expression to flag ICD codes for Acute myocardial infarction
}
  \item{regex_chf}{
A character string containing a regular expression to flag ICD codes for Congestive heart failure
}
  \item{regex_pvd}{
A character string containing a regular expression to flag ICD codes for Peripheral vascular disease
}
  \item{regex_cva}{
A character string containing a regular expression to flag ICD codes for Cerebral vascular disease
}
  \item{regex_dementia}{
A character string containing a regular expression to flag ICD codes for Dementia
}
  \item{regex_pulmonary_disease}{
A character string containing a regular expression to flag ICD codes for Pulmonary disease
}
  \item{regex_ctd}{
A character string containing a regular expression to flag ICD codes for Connective tissue disorder
}
  \item{regex_peptic_ulcer}{
A character string containing a regular expression to flag ICD codes for Peptic ulcer
}
  \item{regex_liver_disease}{
A character string containing a regular expression to flag ICD codes for Liver disease
}
  \item{regex_diabetes}{
A character string containing a regular expression to flag ICD codes for Diabetes
}
  \item{regex_diabetes_complications}{
A character string containing a regular expression to flag ICD codes for Diabetes complications
}
  \item{regex_paraplegia}{
A character string containing a regular expression to flag ICD codes for Paraplegia
}
  \item{regex_renal_disease}{
A character string containing a regular expression to flag ICD codes for Renal disease
}
  \item{regex_cancer}{
A character string containing a regular expression to flag ICD codes for Cancer
}
  \item{regex_metastatic_cancer}{
A character string containing a regular expression to flag ICD codes for Metastatic cancer
}
  \item{regex_severe_liver_disease}{
A character string containing a regular expression to flag ICD codes for Severe liver disease
}
  \item{regex_hiv}{
A character string containing a regular expression to flag ICD codes for HIV
}
  \item{icdcodes_ami}{
A character vector containing ICD codes to be matched exactly to flag ICD codes for Acute myocardial infarction. Overrides regular expression
}
  \item{icdcodes_chf}{
A character vector containing ICD codes to be matched exactly to flag ICD codes for Congestive heart failure. Overrides regular expression
}
  \item{icdcodes_pvd}{
A character vector containing ICD codes to be matched exactly to flag ICD codes for Peripheral vascular disease. Overrides regular expression
}
  \item{icdcodes_cva}{
A character vector containing ICD codes to be matched exactly to flag ICD codes for Cerebral vascular disease. Overrides regular expression
}
  \item{icdcodes_dementia}{
A character vector containing ICD codes to be matched exactly to flag ICD codes for Dementia. Overrides regular expression
}
  \item{icdcodes_pulmonary_disease}{
A character vector containing ICD codes to be matched exactly  to flag ICD codes for Pulmonary disease. Overrides regular expression
}
  \item{icdcodes_ctd}{
A character vector containing ICD codes to be matched exactly to flag ICD codes for Connective tissue disorder. Overrides regular expression
}
  \item{icdcodes_peptic_ulcer}{
A character vector containing ICD codes to be matched exactly to flag ICD codes for Peptic ulcer. Overrides regular expression
}
  \item{icdcodes_liver_disease}{
A character vector containing ICD codes to be matched exactly to flag ICD codes for Liver disease. Overrides regular expression
}
  \item{icdcodes_diabetes}{
A character vector containing ICD codes to be matched exactly to flag ICD codes for Diabetes. Overrides regular expression
}
  \item{icdcodes_diabetes_complications}{
A character vector containing ICD codes to be matched exactly to flag ICD codes for Diabetes complications. Overrides regular expression
}
  \item{icdcodes_paraplegia}{
A character vector containing ICD codes to be matched exactly to flag ICD codes for Paraplegia. Overrides regular expression
}
  \item{icdcodes_renal_disease}{
A character vector containing ICD codes to be matched exactly to flag ICD codes for Renal disease. Overrides regular expression
}
  \item{icdcodes_cancer}{
A character vector containing ICD codes to be matched exactly to flag ICD codes for Cancer. Overrides regular expression
}
  \item{icdcodes_metastatic_cancer}{
A character vector containing ICD codes to be matched exactly to flag ICD codes for Metastatic cancer. Overrides regular expression
}
  \item{icdcodes_severe_liver_disease}{
A character vector containing ICD codes to be matched exactly to flag ICD codes for Severe liver disease. Overrides regular expression
}
  \item{icdcodes_hiv}{
A character vector containing ICD codes to be matched exactly to flag ICD codes for HIV. Overrides regular expression
}
  \item{weight_ami}{
The weight that should be applied to the comorbidity flag for Acute myocardial infarction
}
  \item{weight_chf}{
The weight that should be applied to the comorbidity flag for Congestive heart failure
}
  \item{weight_pvd}{
The weight that should be applied to the comorbidity flag for Peripheral vascular disease
}
  \item{weight_cva}{
The weight that should be applied to the comorbidity flag for Cerebral vascular disease
}
  \item{weight_dementia}{
The weight that should be applied to the comorbidity flag for Dementia
}
  \item{weight_pulmonary_disease}{
The weight that should be applied to the comorbidity flag for Pulmonary disease
}
  \item{weight_ctd}{
The weight that should be applied to the comorbidity flag for Connective tissue disorder
}
  \item{weight_peptic_ulcer}{
The weight that should be applied to the comorbidity flag for Peptic ulcer
}
  \item{weight_liver_disease}{
The weight that should be applied to the comorbidity flag for Liver disease
}
  \item{weight_diabetes}{
The weight that should be applied to the comorbidity flag for Diabetes
}
  \item{weight_diabetes_complications}{
The weight that should be applied to the comorbidity flag for Diabetes complications
}
  \item{weight_paraplegia}{
The weight that should be applied to the comorbidity flag for Paraplegia
}
  \item{weight_renal_disease}{
The weight that should be applied to the comorbidity flag for Renal disease
}
  \item{weight_cancer}{
The weight that should be applied to the comorbidity flag for Cancer
}
  \item{weight_metastatic_cancer}{
The weight that should be applied to the comorbidity flag for Metastatic cancer
}
  \item{weight_severe_liver_disease}{
The weight that should be applied to the comorbidity flag for Severe liver disease
}
  \item{weight_hiv}{
The weight that should be applied to the comorbidity flag for HIV
}
}


\examples{
example_icdcodes_data

example_icdcodes_data \%>\% charlson()
}


\references{
Code adapted from 'Module 5 exercise solutions.SAS' from University of Sydney course PUBH5215 Analysis of Linked Health Data

Example dataset created using comorbidity::sample_diag()
Gasparini A, Salmasian H, Williman J 2020
https://cran.r-project.org/web/packages/comorbidity/index.html

Charlson ME, Pompei P, Ales KL, et al. A new method of classifying prognostic comorbidity in longitudinal studies: development and validation. Journal of Chronic Diseases 1987; 40:373-383.
}

\author{
James Hedley

Adpated from code written by Judy Simpson, Patrick Kelly, and Erin Cvejic
}

