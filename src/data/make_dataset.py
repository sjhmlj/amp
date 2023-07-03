import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
from plotly.subplots import make_subplots
from torch.utils.data import Dataset, DataLoader
import plotly.graph_objects as go
from IPython.display import display
print(repr(Dataset))

# for EDA, merging data
class ShowData(object):
    def __init__(self, path, pep_df, pro_df, cli_df, suppl_df=None):
        self.pep = pd.read_csv(f'{path}/{pep_df}')
        self.pro = pd.read_csv(f'{path}/{pro_df}')
        self.cli = pd.read_csv(f'{path}/{cli_df}')
        self.suppl = pd.read_csv(f'{path}/{suppl_df}')
        self.df_name_list = ['peptides', 'proteins', 'clinical', 'supplemental']
        self.data = [self.pep, self.pro, self.cli, self.suppl]

    def show_shape(self):
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

    def display_head(self):
        df_head_list = []
        for i in self.data:
            df_head_list.append(i.head())
        return df_head_list

    def plot_null_count(self, df_num=2, x='visit_month', y='upd23b_clinical_state_on_medication', width=800, height=500, title_font_size=20):
        ''' df_num to distinguish which df is needed'''
        df = self.data[df_num].copy()
        df[y] = df[y].fillna('Null')

        fig = px.histogram(df, x=x, color=y, title=f'count of {x}', color_discrete_sequence=px.colors.qualitative.Vivid,
                           width=width, height=height)
        fig.update_layout(template='plotly_dark')
        fig.update_layout(title_font_size=title_font_size)
        fig.show()

    def plot_null_count_abt_updrs(self, df_num=2, y='upd23b_clinical_state_on_medication', width=800, height=500, title_font_size=20):
        df = self.data[df_num].copy()
        df[y] = df[y].fillna('Null')

        fig = make_subplots(rows=4, cols=1)
        for i, updrs_part in enumerate(['1', '2', '3', '4']):
            # fig_ = px.histogram(df, x=f'updrs_{updrs_part}', color=y, title=f'count of updrs part {updrs_part} score ({self.df_name_list[df_num]})', color_discrete_sequence=px.colors.qualitative.Vivid,
            #                width=width, height=height)
            for j in ['On', 'Off', 'Null']:
                df_ = df.query('upd23b_clinical_state_on_medication in @j')
                fig_ = go.Histogram(x=df[f'updrs_{updrs_part}'], y=df_[y])
                fig.add_trace(fig_, row=i+1, col=1)
            fig.update_yaxes(title_text=f'UPDRS Part {updrs_part}', row=i+1)
        fig.update_layout(template='plotly_dark')
        fig.update_layout(title_font_size=title_font_size)
        fig.show()

    def plot_updrs_on_md(self, ):
        df = self.data[2].copy()
        df['upd23b_clinical_state_on_medication'] = df['upd23b_clinical_state_on_medication'].fillna('Null')

        parts = ['1', '2', '3', '4']
        fig = make_subplots(rows=4, cols=2,
                            horizontal_spacing=0.10,
                            vertical_spacing=0.06,
                            column_titles=['Medication: "On"', 'Medication: "Off" or "Null"'])
        for i, medication in enumerate([['On'], ['Off', 'Null']]):
            df_ = df.query('upd23b_clinical_state_on_medication in @medication')
            for j, part in enumerate(parts):
                fig.add_trace(go.Box(x=df_['visit_month'], y=df_[f'updrs_{part}']), row=j+1, col=i+1)
                fig.update_xaxes(title_text='Visit Month', row=j+1, col=i+1)
                fig.update_yaxes(title_text=f'UPDRS Part {part}', row=j+1, col=1)
        fig.update_layout(template='plotly_dark',
                          width=850, height=2000,
                          title_text='<b>UPDRS Score - Visit Month',
                          showlegend=False)
        fig.show()



class DataPrep(object):
    pass

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
