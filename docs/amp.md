## Data analysis

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

## Data information

- updrs
  - https://www.movementdisorders.org/MDS/MDS-Rating-Scales/MDS-Unified-Parkinsons-Disease-Rating-Scale-MDS-UPDRS.htm
  - https://www.movementdisorders.org/MDS-Files1/PDFs/Rating-Scales/MDS-UPDRS_English_FINAL.pdf
  - part 1
    - Non-Motor Aspects of Experiences of Daily Living
    - this portion of scale assesses the non-motor impact of PD on patient's experiences of daily living. There are 13 questions. Part 1A is administered by the rater(six questions) and focuses on complex behaviors. Part 1B is a component of self-administered Patient Questionnaire that covers seven questions on non-motor experiences of daily living. 
    - like cognitive impairment, depressed mood, etc. 
  - part 2
    - Motor Aspects of Experiences of Daily Living
    - like speech, eating tasks, dressing, hygiene, handwriting, etc.
  - Part 3
    - Motor Examination
    - at the top of the form, mark whether the patient is on medication for treating the symptoms of PD and, if on levodopa, the time since the last dose. 
    - also, if the patient is receiving medication for treating the symptoms of PD, mark the patient's clinical state using the following definitions
      - ON is the typical functional state when patients are receiving medication and have a good response. 
      - OFF is the typical functional state when patients have a poor response in spite of taking medications.
  - part 4
    - in this section, the rater uses historical and objective information to assess two motor complications, dyskinesias and motor fluctuations that include OFF-state dystonia. 
- upd23b_clinical_state_on_medication
  - Whether or not the patient was taking medication such as Levodopa during the UPDRS assessment. Expected to mainly affect the scores for Part 3 (motor function). These medications wear off fairly quickly (on the order of one day) so it's common for patients to take the motor function exam twice in a single month, both with and without medication.
- Supplemental_clinical_data
  - Clinical records without any associated CSF samples. This data is intended to provide additional context about the typical progression of Parkinsons.

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
  
- 7/21
  
  - Peptide도 std, mean 다르니까 이거 변환하고 pca 돌려보기
    - 보통 standardization 하면 std가 -4 ~ 4 사이가 나온다고 한다(ando)
    - Data 제공자 측에 preprocessing 된 건지 알아보기(max)
  
  - Updrs n등분해서 pca 지켜보기
  
  - Updrs 파트별로 따로 pca 그려보기

## Q

- who's task is it to finds relationship and effects of variables? - human, or machine?
- 