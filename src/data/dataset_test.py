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
data.coeff_var_pep()
