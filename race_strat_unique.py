import os
os.environ["R_HOME"] = r"C:\Program Files\R\R-4.0.2"
os.environ["PATH"]   = r"C:\Program Files\R\R-4.0.2\bin\x64" + ";" + os.environ["PATH"]

import pandas as pd
import numpy as np
from scipy import stats
import clarite

black_bonf = pd.read_csv("Significant_Results_Bonferroni_0.05_std_black.txt", sep="\t")
mexican_bonf = pd.read_csv("Significant_Results_Bonferroni_0.05_std_mexican.txt", sep="\t")
white_bonf = pd.read_csv("Significant_Results_Bonferroni_0.05_std_white.txt", sep="\t")
org_bonf = pd.read_csv("Significant_Results_Bonferroni_0.05_std.txt", sep="\t")

black_FDR = pd.read_csv("Significant_Results_FDR_0.1_black.txt", sep="\t")
mexican_FDR = pd.read_csv("Significant_Results_FDR_0.1_mexican.txt", sep="\t")
white_FDR = pd.read_csv("Significant_Results_FDR_0.1_white.txt", sep="\t")

merged_inner_wm = pd.merge(left=white_bonf, right=mexican_bonf, left_on=['Variable', 'Outcome'], right_on=['Variable', 'Outcome'])
merged_inner_bm = pd.merge(left=black_bonf, right=mexican_bonf, left_on=['Variable', 'Outcome'], right_on=['Variable', 'Outcome'])
merged_inner_wb = pd.merge(left=white_bonf, right=black_bonf, left_on=['Variable', 'Outcome'], right_on=['Variable', 'Outcome'])
merged_inner_bmw = pd.merge(left=merged_inner_bm, right=white_bonf, left_on=['Variable', 'Outcome'], right_on=['Variable', 'Outcome'])
 
comb_bm_bwm = pd.concat([merged_inner_bmw, merged_inner_bm], ignore_index=True, axis=0)
unique_bm = comb_bm_bwm.drop_duplicates(subset=['Variable', 'Outcome'], keep=False)
unique_bm.to_csv("All_groups_bm_std.txt", sep="\t")

comb_wm_bwm = pd.concat([merged_inner_bmw, merged_inner_wm], ignore_index=True, axis=0)
unique_wm = comb_wm_bwm.drop_duplicates(subset=['Variable', 'Outcome'], keep=False)
unique_wm.to_csv("All_groups_wm_std.txt", sep="\t")

comb_wb_bwm = pd.concat([merged_inner_bmw, merged_inner_wb], ignore_index=True, axis=0)
unique_wb = comb_wb_bwm.drop_duplicates(subset=['Variable', 'Outcome'], keep=False)
unique_wb.to_csv("All_groups_wb_std.txt", sep="\t")

merged_inner_bmw.to_csv("merged_inner_bmw_std.txt", sep="\t")


#add column for race
black_bonf = black_bonf.assign(Group='AA')
white_bonf = white_bonf.assign(Group='EA')
mexican_bonf = mexican_bonf.assign(Group='MA')
org_bonf = org_bonf.assign(Group="Original")

merged_inner_wo = pd.merge(left=white_bonf, right=org_bonf, left_on=['Variable', 'Outcome'], right_on=['Variable', 'Outcome'])

unique_aa = black_bonf.merge(org_bonf.drop_duplicates(), on=['Variable','Outcome'], 
                   how='left', indicator=True)

comb_ma_aa = pd.concat([unique_ma, unique_aa], ignore_index=False, axis=0)
comb_all = pd.concat([comb_ma_aa, unique_ea], ignore_index=False, axis=0)


#rename varaible and outcome columns to match names in the data dictionary
VarDescriptioncsv.loc[VarDescriptioncsv['var'] == 'LBXIRN', 'var_desc'] = 'Iron, Frozen Serum (ug/dL)'
VarDescriptioncsv.loc[VarDescriptioncsv['var'] == 'LBXSIR', 'var_desc'] = 'Iron, Count (ug/dL)'

data_dict = VarDescriptioncsv.set_index("var")["var_desc"].to_dict()
comb_all["var"] = comb_all["Variable"].apply(lambda s: data_dict.get(s, ""))
comb_all["outcome"] = comb_all["Outcome"].apply(lambda s: data_dict.get(s, ""))

comb_all = comb_all.drop(['Variable', 'Outcome'], axis = 1)
comb_all = comb_all.drop(comb_all.columns[12:24], axis = 1)
comb_all.rename(columns = {'var' : 'Variable', 'outcome' : 'Outcome'}, inplace = True)
first_column = comb_all.pop('Variable')
second_column = comb_all.pop('Outcome')
comb_all.insert(0, 'Variable', first_column)
comb_all.insert(1, 'Outcome', second_column)
comb_all.replace(to_replace="Vitamin_A (ug/dL)", value="Retinol (Vitamin A) (ug/dL)")
comb_all["Variable"].replace({"g-Tocopherol (ug/dL)":"gamma-tocopherol (Vitamin E) (ug/dl)"}, inplace=True)
comb_all.to_csv("comb_all.txt", na_rep="NA", sep="\t")
