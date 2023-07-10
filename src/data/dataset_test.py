from make_dataset import ShowData, CustomDataset
from IPython.display import display

path = "/home/jh/PycharmProjects/amp/data"

data = ShowData(
    path,
    "train_peptides.csv",
    "train_proteins.csv",
    "train_clinical_data.csv",
    "supplemental_clinical_data.csv",
)

# print(data.show_df_shape(), data.show_null_count(), sep='\n')
# display(data.display_df_head())
# data.plot_null_count(x='visit_month', y='upd23b_clinical_state_on_medication')

# print(data.show_null_count())
# data.plot_updrs_on_md()
# data.plot_null_count_abt_updrs()
# data.coeff_var_pep()
# data.coeff_var_pep_on_md()
data.corr_top5pep_updrs()
# data.coeff_var_pro()
# data.coeff_var_pro_on_md()
# data.corr_top5pro_updrs()
# data.show_cli_suppl()
# data.show_updrs_sum_ab_visitmonth()
# data.show_updrs_sum_ab_patientid()
# data.show_updrs_mean_ab_md()
# data.plot_patient_updrs_info()
# print(dir(ShowData))
# data.plot_md_ab_visitmonth()
# data.about_data(all_data=True)