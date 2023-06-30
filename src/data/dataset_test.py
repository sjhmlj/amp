from make_dataset import DataPrep, CustomDataset

path = '/home/jh/PycharmProjects/amp/data'

data = DataPrep(path, 'train_peptides.csv', 'train_proteins.csv', 'train_clinical_data.csv', 'supplemental_clinical_data.csv')

print(data.show_df_shape(), data.show_null_count(), sep='\n')