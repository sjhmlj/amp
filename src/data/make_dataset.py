import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
from plotly.subplots import make_subplots
from torch.utils.data import Dataset, DataLoader
import plotly.graph_objects as go
from IPython.display import display

# print(repr(Dataset))

# for EDA, merging data
class ShowData(object):
    def __init__(self, path, pep_df, pro_df, cli_df, suppl_df=None):
        self.pep = pd.read_csv(f"{path}/{pep_df}")
        self.pro = pd.read_csv(f"{path}/{pro_df}")
        self.cli = pd.read_csv(f"{path}/{cli_df}")
        self.suppl = pd.read_csv(f"{path}/{suppl_df}")
        self.df_name_list = ["peptides", "proteins", "clinical", "supplemental"]
        self.data = [self.pep, self.pro, self.cli, self.suppl]
        self.abundance_cv_mean_pep = None
        self.abundance_cv_mean_pro = None

    def about_data(self, df_num=0, all_data=False):
        if all_data == True:
            data = self.data
            for i, item in enumerate(data):
                print('data: ' + self.df_name_list[i])
                display(item.head())
                print(self.show_col(i))
                print('shape: ' + str(self.show_shape(i)))
                print('null count: ', self.show_null_count(i), sep='\n')
                print()
        else:
            data = self.data[df_num]
            print('data: ' + self.df_name_list[df_num])
            display(data.head())
            print(self.show_col(df_num))
            print('shape: ' + str(self.show_shape(df_num)))
            print('null count: ', self.show_null_count(df_num), sep='\n')
    def show_shape(self, df_num=0):
        df = self.data[df_num]
        return df.shape


    def show_null_count(self, df_num=0):
        df = self.data[df_num]
        return df.isna().sum()

    def show_col(self, df_num=0):
        df = self.data[df_num]
        return df.columns

    def plot_null_count(
        self,
        df_num=2,
        x="visit_month",
        y="upd23b_clinical_state_on_medication",
        width=800,
        height=500,
        title_font_size=20,
    ):
        """df_num to distinguish which df is needed"""
        df = self.data[df_num].copy()
        df[y] = df[y].fillna("Null")

        fig = px.histogram(
            df,
            x=x,
            color=y,
            title=f"count of {x} - cli",
            color_discrete_sequence=px.colors.qualitative.Vivid,
            width=width,
            height=height,
        )
        fig.update_layout(template="plotly_dark")
        fig.update_layout(title_font_size=title_font_size)
        fig.show()

    def plot_null_count_abt_updrs(
        self,
        df_num=2,
        y="upd23b_clinical_state_on_medication",
        width=800,
        height=500,
        title_font_size=20,
    ):
        df = self.data[df_num].copy()
        df[y] = df[y].fillna("Null")

        fig = make_subplots(rows=4, cols=1)
        for i, updrs_part in enumerate(["1", "2", "3", "4"]):
            # fig_ = px.histogram(df, x=f'updrs_{updrs_part}', color=y, title=f'count of updrs part {updrs_part} score ({self.df_name_list[df_num]})', color_discrete_sequence=px.colors.qualitative.Vivid,
            #                width=width, height=height)
            for j in ["On", "Off", "Null"]:
                df_ = df.query("upd23b_clinical_state_on_medication in @j")
                fig_ = go.Histogram(x=df[f"updrs_{updrs_part}"], y=df_[y])
                fig.add_trace(fig_, row=i + 1, col=1)
            fig.update_yaxes(title_text=f"UPDRS Part {updrs_part}", row=i + 1)
        fig.update_layout(template="plotly_dark")
        fig.update_layout(title='state_on_medication: On, Off, Null')
        fig.update_layout(title_font_size=title_font_size)
        fig.show()

    def plot_updrs_on_md(
        self,
    ):
        df = self.data[2].copy()
        df["upd23b_clinical_state_on_medication"] = df[
            "upd23b_clinical_state_on_medication"
        ].fillna("Null")

        parts = ["1", "2", "3", "4"]
        fig = make_subplots(
            rows=4,
            cols=2,
            horizontal_spacing=0.10,
            vertical_spacing=0.06,
            column_titles=['Medication: "On"', 'Medication: "Off" or "Null"'],
        )
        for i, medication in enumerate([["On"], ["Off", "Null"]]):
            df_ = df.query("upd23b_clinical_state_on_medication in @medication")
            for j, part in enumerate(parts):
                fig.add_trace(
                    go.Box(x=df_["visit_month"], y=df_[f"updrs_{part}"]),
                    row=j + 1,
                    col=i + 1,
                )
                fig.update_xaxes(title_text="Visit Month", row=j + 1, col=i + 1)
                fig.update_yaxes(title_text=f"UPDRS Part {part}", row=j + 1, col=1)
        fig.update_layout(
            template="plotly_dark",
            width=850,
            height=2000,
            title_text="<b>UPDRS Score - Visit Month",
            showlegend=False,
        )
        fig.show()

    def make_abundance_cv_mean_pep(self):
        if self.abundance_cv_mean_pep is not None:
            return self.abundance_cv_mean_pep

        df = self.pep[["patient_id", "Peptide", "PeptideAbundance"]]
        df_agg = df.groupby(["patient_id", "Peptide"])["PeptideAbundance"].aggregate(
            ["mean", "std"]
        )
        df_agg["CV_PeptideAbundance[%]"] = df_agg["std"] / df_agg["mean"] * 100

        abundance_cv_mean = (
            df_agg.groupby("Peptide")["CV_PeptideAbundance[%]"].mean().reset_index()
        )
        abundance_cv_mean = abundance_cv_mean.sort_values(
            by="CV_PeptideAbundance[%]", ascending=False
        ).reset_index()
        self.abundance_cv_mean_pep = abundance_cv_mean
        return self.abundance_cv_mean_pep
    def coeff_var_pep(self):
        df = self.pep[["patient_id", "Peptide", "PeptideAbundance"]].copy()
        df_agg = df.groupby(["patient_id", "Peptide"])["PeptideAbundance"].aggregate(
            ["mean", "std"]
        )
        df_agg["CV_PeptideAbundance[%]"] = df_agg["std"] / df_agg["mean"] * 100
        abundance_cv_mean = self.make_abundance_cv_mean_pep()
        peptide_cv_top5 = abundance_cv_mean[:5]["Peptide"]

        df_agg_top5 = df_agg.query("Peptide in @peptide_cv_top5").reset_index()
        df_agg_top5["order"] = 0
        for i, peptide in enumerate(peptide_cv_top5):

            df_agg_top5.loc[df_agg_top5["Peptide"] == peptide, 'order'] = i
        df_agg_top5.sort_values(by="order", inplace=True)

        fig = px.violin(
            df_agg_top5,
            y="Peptide",
            x="CV_PeptideAbundance[%]",
            color="Peptide",
            box=True,
            title="<b>Coefficient of Variation (top 5) of PeptideAbundance per patient_id",
            width=900,
            height=800,
        )
        fig.update_layout(
            template="plotly_dark",
            width=800,
            height=500,
            showlegend=False,
            xaxis=dict(
                title="Coefficient of Variation [%] of PeptideAbundance per patient_id"
            ),
            yaxis=dict(title="<b>Peptide"),
        )
        fig.show()
        

    def coeff_var_pep_on_md(self):
        abundance_cv_mean = self.make_abundance_cv_mean_pep()
        df = pd.merge(
            self.pep,
            self.cli[["visit_id", "upd23b_clinical_state_on_medication"]],
            on="visit_id",
        )
        df["upd23b_clinical_state_on_medication"] = df[
            "upd23b_clinical_state_on_medication"
        ].fillna("Null")
        peptide_cv_top5 = abundance_cv_mean[:5]["Peptide"].reset_index(drop=True)
        tmp_df = df[
            [
                "visit_month",
                "patient_id",
                "Peptide",
                "PeptideAbundance",
                "upd23b_clinical_state_on_medication",
            ]
        ]
        top5_pep_df = tmp_df.query("Peptide in @peptide_cv_top5")

        num_patient = 50
        fig = make_subplots(
            rows=5,
            cols=2,
            horizontal_spacing=0.10,
            vertical_spacing=0.06,
            column_titles=["<b>Medication: On", "<b>Medication: Off or Null"],
            shared_yaxes=True,
        )
        for i, medication in enumerate([["On"], ["Off", "Null"]]):
            top5_pep_df_ = top5_pep_df.query(
                "upd23b_clinical_state_on_medication in @medication"
            )
            for j, peptide in enumerate(peptide_cv_top5):
                peptide_df = top5_pep_df_.query("Peptide == @peptide")  ## 이 부분 다르게 했음
                fig.add_trace(
                    go.Box(
                        x=peptide_df["visit_month"], y=peptide_df["PeptideAbundance"]
                    ),
                    row=j + 1,
                    col=i + 1,
                )
                fig.update_xaxes(title_text="Visit Month", row=j + 1, col=i + 1)
                fig.update_yaxes(title_text=f"{peptide} Abundance", row=j + 1, col=1)
        fig.update_layout(
            template="plotly_dark",
            width=800,
            height=2000,
            title_text="<b>Peptide Abundance - Visit Month (Highest top5 CV)",
            showlegend=False,
        )
        fig.show()

    def corr_top5pep_updrs(self):
        num_peptide_candidates = 5
        pep_candidates = self.make_abundance_cv_mean_pep().loc[:num_peptide_candidates-1, 'Peptide']

        pep_candidate_df = self.pep.query('Peptide in @pep_candidates')
        visit_ids = pep_candidate_df.visit_id.unique()
        pep_dict_list = []

        for visit_id in visit_ids:
            pep_df = pep_candidate_df.query(f'visit_id=="{visit_id}"')
            peptides = pep_df['Peptide'].values
            PeptideAbundances = pep_df['PeptideAbundance'].values
            peptide_dict = dict(zip(pep_candidates, [np.nan]*num_peptide_candidates))
            for peptide, PeptideAbundance in zip(peptides, PeptideAbundances):
                peptide_dict[peptide] = PeptideAbundance
            peptide_dict['visit_id'] = visit_id
            pep_dict_list.append(peptide_dict)

        pep_candidate_df2 = pd.DataFrame(pep_dict_list)
        cli_pep_df = pd.merge(self.cli, pep_candidate_df2, on='visit_id')
        cli_pep_df_ = cli_pep_df.drop(['visit_id', 'patient_id', 'visit_month', 'upd23b_clinical_state_on_medication'], axis=1)
        corr = cli_pep_df_.corr()
        fig = px.imshow(corr, width=750, height=750, title='<b>Correlation: UPDR scores and Peptide Abundance')
        fig.show()

    def make_abundance_cv_mean_pro(self):
        if self.abundance_cv_mean_pro is not None:
            return self.abundance_cv_mean_pro
        df = self.pro[["patient_id", "UniProt", "NPX"]]
        df_agg = df.groupby(["patient_id", "UniProt"])["NPX"].aggregate(
            ["mean", "std"]
        )
        df_agg["CV_NPX[%]"] = df_agg["std"] / df_agg["mean"] * 100

        abundance_cv_mean = (
            df_agg.groupby("UniProt")["CV_NPX[%]"].mean().reset_index()
        )
        abundance_cv_mean = abundance_cv_mean.sort_values(
            by="CV_NPX[%]", ascending=False
        ).reset_index()
        self.abundance_cv_mean_pro = abundance_cv_mean
        return self.abundance_cv_mean_pro
    def coeff_var_pro(self):
        df = self.pro[["patient_id", "UniProt", "NPX"]]
        df_agg = df.groupby(["patient_id", "UniProt"])["NPX"].aggregate(
            ["mean", "std"]
        )
        df_agg["CV_NPX[%]"] = df_agg["std"] / df_agg["mean"] * 100
        abundance_cv_mean = self.make_abundance_cv_mean_pro()
        protein_cv_top5 = abundance_cv_mean[:5]["UniProt"]

        df_agg_top5 = df_agg.query("UniProt in @protein_cv_top5").reset_index()
        df_agg_top5["order"] = 0
        for i, protein in enumerate(protein_cv_top5):
            df_agg_top5.loc[df_agg_top5["UniProt"] == protein, 'order'] = i
        df_agg_top5.sort_values(by="order", inplace=True)

        fig = px.violin(
            df_agg_top5,
            y="UniProt",
            x="CV_NPX[%]",
            color="UniProt",
            box=True,
            title="<b>Coefficient of Variation (top 5) of NPX per patient_id",
            width=900,
            height=800,
        )
        fig.update_layout(
            template="plotly_dark",
            width=800,
            height=500,
            showlegend=False,
            xaxis=dict(
                title="Coefficient of Variation [%] of NPX per patient_id"
            ),
            yaxis=dict(title="<b>UniProt"),
        )
        fig.show()

    def coeff_var_pro_on_md(self):
        abundance_cv_mean = self.make_abundance_cv_mean_pro()
        df = pd.merge(
            self.pro,
            self.cli[["visit_id", "upd23b_clinical_state_on_medication"]],
            on="visit_id",
        )
        df["upd23b_clinical_state_on_medication"] = df[
            "upd23b_clinical_state_on_medication"
        ].fillna("Null")
        protein_cv_top5 = abundance_cv_mean[:5]["UniProt"].reset_index(drop=True)
        tmp_df = df[
            [
                "visit_month",
                "patient_id",
                "UniProt",
                "NPX",
                "upd23b_clinical_state_on_medication",
            ]
        ]
        top5_pro_df = tmp_df.query("UniProt in @protein_cv_top5")

        num_patient = 50
        fig = make_subplots(
            rows=5,
            cols=2,
            horizontal_spacing=0.10,
            vertical_spacing=0.06,
            column_titles=["<b>Medication: On", "<b>Medication: Off or Null"],
            shared_yaxes=True,
        )
        for i, medication in enumerate([["On"], ["Off", "Null"]]):
            top5_pro_df_ = top5_pro_df.query(
                "upd23b_clinical_state_on_medication in @medication"
            )
            for j, protein in enumerate(protein_cv_top5):
                protein_df = top5_pro_df_.query("UniProt == @protein")
                fig.add_trace(
                    go.Box(
                        x=protein_df["visit_month"], y=protein_df["NPX"]
                    ),
                    row=j + 1,
                    col=i + 1,
                )
                fig.update_xaxes(title_text="Visit Month", row=j + 1, col=i + 1)
                fig.update_yaxes(title_text=f"{protein} Abundance", row=j + 1, col=1)
        fig.update_layout(
            template="plotly_dark",
            width=800,
            height=2000,
            title_text="<b>NPX - Visit Month (Highest top5 CV)",
            showlegend=False,
        )
        fig.show()

    def corr_top5pro_updrs(self):
        num_protein_candidates = 5
        pro_candidates = self.make_abundance_cv_mean_pro().loc[:num_protein_candidates-1, 'UniProt']

        pro_candidate_df = self.pro.query('UniProt in @pro_candidates')
        visit_ids = pro_candidate_df.visit_id.unique()
        pro_dict_list = []
        for visit_id in visit_ids:
            pro_df = pro_candidate_df.query(f'visit_id=="{visit_id}"')
            proteins = pro_df['UniProt'].values
            NPXs = pro_df['NPX'].values
            protein_dict = dict(zip(pro_candidates, [np.nan]*num_protein_candidates))
            for protein, NPX in zip(proteins, NPXs):
                protein_dict[protein] = NPX
            protein_dict['visit_id'] = visit_id
            pro_dict_list.append(protein_dict)

        pro_candidate_df2 = pd.DataFrame(pro_dict_list)
        cli_pro_df = pd.merge(self.cli, pro_candidate_df2, on='visit_id')
        cli_pro_df_ = cli_pro_df.drop(['visit_id', 'patient_id', 'visit_month', 'upd23b_clinical_state_on_medication'], axis=1)
        corr = cli_pro_df_.corr()
        fig = px.imshow(corr, width=750, height=750, title='<b>Correlation: UPDR scores and Peptide Abundance')
        fig.show()

    def show_cli_suppl(self):

        med = 'upd23b_clinical_state_on_medication'
        cli, suppl = self.cli.copy(), self.suppl.copy()
        cli[med], suppl[med] = cli[med].fillna('Null'), suppl[med].fillna('Null')
        all_df = pd.concat([cli, suppl]).reset_index(drop=True)
        all_df['cli_or_suppl'] = 'cli'
        all_df.loc[len(cli):, 'cli_or_suppl'] = 'suppl'
        display(all_df)
        print(f'# train_clinical_data: {len(self.cli)}')
        print(f'# supplemental_clinical_data: {len(self.suppl)}')

        features = ['visit_month', 'upd23b_clinical_state_on_medication', 'updrs_1', 'updrs_2',
        'updrs_3', 'updrs_4']

        fig = make_subplots(rows=len(features),horizontal_spacing=0.10,
            vertical_spacing=0.06,)
        for i, x in enumerate(features):
            for cli_or_suppl in ['cli', 'suppl']:
                all_df_ = all_df.query('cli_or_suppl == @cli_or_suppl')
                fig.add_trace(
                    go.Histogram(x=all_df_[x], name=f'{cli_or_suppl}'), row=i+1, col=1
                )
            fig.update_yaxes(title_text=f'<b>Count of {x}', row=i + 1)
        fig.update_layout(template='plotly_dark', width=1600, height=2000)
        fig.update_layout(title_font_size=20)

        fig.show()

    def show_updrs_sum_ab_visitmonth(self):
        cli, suppl = self.cli.copy(), self.suppl.copy()
        cli['updrs_sum'] = 0
        suppl['updrs_sum'] = 0
        for i in range(1, 5):
            cli['updrs_sum'] += cli[f'updrs_{i}']
            suppl['updrs_sum'] += suppl[f'updrs_{i}']
        cli_ = cli.groupby('visit_month')['updrs_sum'].mean()
        suppl_ = suppl.groupby('visit_month')['updrs_sum'].mean()

        fig, axs = plt.subplots(nrows=2)

        fig.suptitle('updrs sum about visit month')
        for j, df in enumerate([cli_, suppl_]):
            # sns.histplot(data=pd.DataFrame({'visit_month':df.index, 'sum':df.values}), x='visit_month', y='sum', ax=axs[j], discrete=True)
            sns.barplot(data=pd.DataFrame({'visit_month':df.index, 'sum':df.values}), x='visit_month', y='sum', ax=axs[j],color='skyblue')
        plt.tight_layout()
        plt.show()

    def show_updrs_sum_ab_patientid(self):
        cli, suppl = self.cli.copy(), self.suppl.copy()
        cli, suppl = cli.fillna(0), suppl.fillna(0)
        cli['updrs_sum'] = 0
        suppl['updrs_sum'] = 0

        for i in range(1, 5):
            cli['updrs_sum'] += cli[f'updrs_{i}']
            suppl['updrs_sum'] += suppl[f'updrs_{i}']

        cli_ = cli.groupby('patient_id')['updrs_sum'].mean()
        suppl_ = suppl.groupby('patient_id')['updrs_sum'].mean()

        fig, axs = plt.subplots(nrows=2)
        fig.suptitle('updrs sum about patient id')
        for j, df in enumerate([cli_, suppl_]):
            sns.histplot(data=df,
                         ax=axs[j], discrete=True)
        plt.tight_layout()
        plt.show()

    def show_updrs_mean_ab_md(self):
        cli, suppl = self.cli.copy(), self.suppl.copy()
        _ = 'upd23b_clinical_state_on_medication'
        cli[_] = cli[_].fillna('Null')
        suppl[_] = suppl[_].fillna('Null')

        fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(12, 12))
        for i, df in enumerate([cli, suppl]):
            for j, md in enumerate(['On', 'Off', 'Null']):
                # df_ = df[df[_] == md]
                # df_ = df_.groupby(f'{_}').mean()
                df_ = df.query('upd23b_clinical_state_on_medication in @md').loc[:,['visit_month', 'updrs_1', \
                                'updrs_2', 'updrs_3', 'updrs_4']]
                df_ = df_.groupby('visit_month').mean().reset_index()

                # sns.lineplot(data=df_[['updrs_1','updrs_2','updrs_3', 'updrs_4']], ax=axs[j][i])
                sns.lineplot(data=pd.melt(df_, ['visit_month']), x='visit_month', y='value', hue='variable', ax=axs[j][i])
                axs[j][i].legend(loc='upper right')
                if i == 0:
                    axs[j][i].set_title(f'cli, medication:{md}')
                else:
                    axs[j][i].set_title(f'suppl, medication:{md}')

        fig.suptitle('updrs mean about medication on, off, null')
        plt.tight_layout()
        plt.show()

    def plot_patient_updrs_info(self, idx=None):
        cli = self.cli.copy()
        if idx is None:
            patients = cli.patient_id.unique()
            patient = np.random.choice(patients)
            idx = patient

        patient_df = cli.loc[cli.patient_id == idx].loc[:,
                     ['visit_month', 'updrs_1', 'updrs_2', 'updrs_3', 'updrs_4',
                     'upd23b_clinical_state_on_medication']]
        patient_df['upd23b_clinical_state_on_medication'].replace({'On':1, 'Off':0}, inplace=True)
        patient_df_ = patient_df.drop(columns='upd23b_clinical_state_on_medication')
        patient_df_ = pd.melt(patient_df_, ['visit_month'])
        sns.lineplot(data=patient_df_, x='visit_month', y='value', hue='variable')
        plt.scatter(x=patient_df['visit_month'], y=patient_df['upd23b_clinical_state_on_medication'], label='on medicine')
        plt.title(f'patient {idx} clinical info')
        plt.legend(loc='upper right')
        plt.grid(True)
        plt.show()

    def plot_md_ab_visitmonth(self):
        cli, suppl = self.cli.copy(), self.suppl.copy()
        _ = 'upd23b_clinical_state_on_medication'
        cli[_] = cli[_].fillna('Null')
        suppl[_] = suppl[_].fillna('Null')

        fig, axs = plt.subplots(1, 2, figsize=(10, 6))
        for j, df in enumerate([cli, suppl]):
            for i, md in enumerate(['On', 'Off', 'Null']):
                df_ = df.query(f'{_} in @md')
                df_.loc[:, ['cnt']] = 1
                df_ = df_[['visit_month', 'cnt']].groupby('visit_month').sum().reset_index()
                axs[j].plot(df_['visit_month'], df_['cnt'], label=f'{md}')
            axs[j].legend()
            if j == 0:
                axs[j].set_title('data: cli')
            else:
                axs[j].set_title('data: suppl')

            axs[j].set_xlabel('visit month')
            axs[j].set_ylabel('count')
        fig.suptitle('medication state on visit month')
        plt.show()


class DataPrep(object):
    pass


class CustomDataset(Dataset):
    def __init__(self, df):
        if type(df) == pd.core.frame.DataFrame:
            self.df = df
        else:
            raise Exception("pandas.DataFrame is needed")

    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):
        return self.df.iloc[idx, :]
