## Data

- peptides
  - columns: visit_id, visit_month, patient_id, UniProt, Peptide, PeptideAbundance
  - shape: 981834 x 6
  - null: 0
  - unique peptides: 968
- proteins
  - columns: visit_id, visit_month, patient_id, UniProt, NPX
  - shape: 232741 x 5
  - null: 0
  - unique proteins: 227
- train_clinical_data
  - columns: visit_id, visit_month, patient_id, updrs_1, updrs_2, updrs_3, updrs_4, upd23b_clinical_state_on_medication
  - shape: 2615 x 8
  - null
    - updrs_1: 1
    - updrs_2: 2
    - updrs_3: 25
    - updrs_4: 1038
      - updrs_4 > 0: 616
      - updrs_4 in case (upd23b ~ : Null) : Null(326%), 0(92%)  
    - upd23b_clinical_state_on_medication: 1327
- supplemental_clinical_data(without peptide, protein data)
  - columns: visit_id, visit_month, patient_id, updrs_1, updrs_2, updrs_3, updrs_4, upd23b_clinical_state_on_medication
  - shape: 2223 x 8
  - null:
    - updrs_1: 213
    - updrs_2: 214
    - updrs_3: 5
    - updrs_4: 928
    - upd23b_clinical_state_on_medication: 1101



## Data Visualization

- histogram - medication state on visit_month 
- histogram - updrs values on medication state
- boxplot - updrs values on medication state 
- violinplot - 5 peptides of high coefficient of variation(on patient_id)
- boxplot - 5 peptides of high coefficient of variation(on patient_id) per visit month

- correlation heatmap - updrs, 5 peptides. 
- violinplot, boxplot, correlation heatmap - 5 proteins ~
- histogram - visit month, medication, updrs - cli, suppl
- barplot - updrs sum per visit_month - cli, suppl
- histogram - updrs sum per patient_id - cli, suppl
- lineplot - updrs mean per medication
- line, scatterplot - updrs, medication for 1 patient
- lineplot - medication state per visit_month - cli, suppl



### sub

- pep, pro
  - top 5 coefficient of variation 
  - correlation with updrs
    - pep: ~0.22. only for updrs_2, updrs_3. almost negative relationship 
    - pro: almost same. there is updrs_1 with Q06481
    - corr higher than 0.7 exists between grouped updrs and pep/pro 
- cli, suppl
  - medication state
    - per visit month
    - updrs values
    - updrs mean
  - updrs sum
    - per visit month
    - per patient id
  - correlation with visit_month
    - all corrs are higher than 0.5
    - but when i excluded visit_month == 12*i, corrs were higher. 

## to do

- 7/14
  - olive
    - exercise machine learning example to know how to preprocess data
    - check each updrs part - min, max, range,
      - and determine whether normalization is needed
  - josh
    - plot pca
    - check frequency of each updrs, on visit month
    - understand correlation and each method(pearson first)

## Q

- who's task is it to finds relationship and effects of variables? - human, or machine?
- 