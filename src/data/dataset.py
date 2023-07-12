import os

import numpy as np
import pandas as pd
from torch.utils.data import Dataset, DataLoader

# print(repr(Dataset))

# for EDA, merging data


class CustomDataset(Dataset):
    def __init__(self, path, peptides, proteins, train_clinical, supplemental_clinical):
        self.pep = pd.read_csv(f'{path}/{peptides}')
        self.pro = pd.read_csv(f'{path}/{proteins}')
        self.cli = pd.read_csv(f'{path}/{train_clinical}')
        self.suppl = pd.read_csv(f'{path}/{supplemental_clinical}')
        self.df = None
    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):
        return self.df.iloc[idx, :]

    def all_df(self):
        return self.pep, self.pro, self.cli, self.suppl
