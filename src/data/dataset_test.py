import matplotlib.pyplot as plt

import eda
from IPython.display import display




# num = 10
# for i in range(num):
#     eda.plot_patient_updrs_info()


# peptides = eda.find_high_corr_pro(baseline=0.2, method='spearman')
# print(peptides)
# def t1d_(itr):
#     lst = []
#     for i in itr:
#         if len(i) != 1:
#             for j in i:
#                 lst.append(j)
#         else:
#             lst.append(i)
#
# print(len(peptides), min(map(lambda x:x[1], t1d_(peptides.values()))))

# df = eda.corr_top5pep_updrs()
# for i in range(1, 5):
#     display(df[f'updrs_{i}'])

# eda.plot_corr_updrs()
# a, b = eda.find_corr_visit_updrs()
# print(a, b, sep='\n')
#
# a, b = eda.find_corr_groupedvisit_updrs()
# print(a,b, sep='\n')

print(eda.find_high_corr_pep_grouped(baseline=0.7))
# print(eda.find_high_corr_pro_grouped(baseline=0.6))


# eda.plot_md_ab_visitmonth()
# eda.show_updrs_sum_ab_patientid()
# eda.show_updrs_sum_ab_visitmonth()
# eda.show_updrs_mean_ab_md()
# eda.show_cli_suppl()