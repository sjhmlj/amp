import os
import pandas as pd
from torch.utils.data import Dataset, DataLoader

print(repr(Dataset))

# for EDA, merging data
class DataPrep(object):
    def __init__(self, path, pep_df, pro_df, cli_df, suppl_df=None):
        self.pep = pd.read_csv(f'{path}/{pep_df}')
        self.pro = pd.read_csv(f'{path}/{pro_df}')
        self.cli = pd.read_csv(f'{path}/{cli_df}')
        self.suppl = pd.read_csv(f'{path}/{suppl_df}')
        self.df_name_list = ['peptides', 'proteins', 'clinical', 'supplemental']
        self.data = [self.pep, self.pro, self.cli, self.suppl]

    def show_df_shape(self):
        df_shape_list = []
        for i in self.data:
            if i is not None:
                df_shape_list.append(i.shape)
        return dict(zip(self.df_name_list, df_shape_list))

    def show_null_count(self):
        null_count_list = []
        for i in self.data:
            if i is not None:
                null_count_list.append((i.isna().sum()))
        return dict(zip(self.df_name_list, null_count_list))

    def




class CustomDataset(Dataset):
    def __init__(self, df):
        if type(df) == pd.core.frame.DataFrame:
            self.df = df
        else:
            raise Exception('pandas.DataFrame is needed')

    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):
        return self.df.iloc[idx, :]
