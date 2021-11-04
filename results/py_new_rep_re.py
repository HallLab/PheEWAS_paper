import os
os.environ["R_HOME"] = r"C:\Program Files\R\R-4.0.2"
os.environ["PATH"]   = r"C:\Program Files\R\R-4.0.2\bin\x64" + ";" + os.environ["PATH"]

import pandas as pd
import numpy as np
from scipy import stats
import clarite





# Define data files
data_main_table_over18 = "MainTable_keepvar_over18.tsv"
data_main_table = "MainTable.csv"
data_var_categories = "VarCat_nopf.txt"

##result tables
data_dis = "Discovery_Results_python.txt"
data_rep = "Replication_Results_python.txt"

nhanes_dis = clarite.load.from_tsv("Discovery_Results_python.txt", index_col=("index"))
nhanes_rep = clarite.load.from_tsv("Replication_Results_python.txt")

# Load Data
nhanes = clarite.load.from_tsv(data_main_table_over18, index_col="ID")


# Load Weight Values (and SDMVPSU/SDMVSTRA)
###########################################
survey_design_discovery = pd.read_csv("weights_discovery.txt", sep="\t")\
                            .rename(columns={'SEQN':'ID'})\
                            .set_index("ID")\
                            .drop(columns="SDDSRVYR")
survey_design_replication = pd.read_csv("weights_replication.txt", sep="\t")\
                            .rename(columns={'SEQN':'ID'})\
                            .set_index("ID")\
                            .drop(columns="SDDSRVYR")

# Get weight mapping data
#########################
var_weights = pd.read_csv("VarWeights_new.csv")
weights_discovery = dict()
weights_replication = dict()
for idx, row in var_weights.iterrows():
    weights_discovery[row['variable_name']] = row['discovery']
    weights_replication[row['variable_name']] = row['replication']

####Get column names of phenotypes
phenotypes = pd.read_csv("all_dis_data.txt", sep="\t")
pheno = list(phenotypes.drop(['ID', 'LBXSCRINV'], axis=1))

####Bring over QC'ed data
discovery_catbin = clarite.load.from_csv("discovery_catbin_final_exp.txt", sep="\t")
discovery_cont_phenotype = clarite.load.from_csv("discovery_cont_phenotype_exp.txt", sep="\t")
discovery_cont_phenotype = discovery_cont_phenotype.drop("LBXSCRINV", axis=1)
replication_catbin = clarite.load.from_csv("replication_catbin_final_exp.txt", sep="\t")
replication_cont_phenotype = clarite.load.from_csv("replication_cont_phenotype_exp.txt", sep="\t")
replication_cont_phenotype = replication_cont_phenotype.drop("LBXSCRINV", axis=1)

####Import and Categorize the files
categorical_vars = pd.read_csv("categorical.txt", sep="\t").drop('ID', axis=1).columns
binary_vars = pd.read_csv("binary.txt", sep="\t").drop(['ID', 'current_loud_noise', 'ordinary_salt', 'no_salt'], axis=1).columns
discovery_catbin = clarite.modify.make_binary(discovery_catbin, only=binary_vars)
discovery_catbin = clarite.modify.make_categorical(discovery_catbin, only=categorical_vars)
replication_catbin = clarite.modify.make_binary(replication_catbin, only=binary_vars)
replication_catbin = clarite.modify.make_categorical(replication_catbin, only=categorical_vars)



####Assign covariates & drop white column
cov = ["other_hispanic", "other_eth", "mexican", "black", "female", "SES_LEVEL", "RIDAGEYR", "SDDSRVYR", "BMXBMI"]
drop_white_dis = discovery_catbin.drop(['white'], axis=1)
drop_white_rep = replication_catbin.drop([ 'white'], axis=1)

####Assign covariates & drop ethnicity columns
covariates = ["female", "SES_LEVEL", "RIDAGEYR", "SDDSRVYR", "BMXBMI"]
no_eth_catbin_black = discovery_catbin.drop(['other_hispanic', 'other_eth', 'mexican', 'white'], axis=1)
no_eth_catbin_black2 = replication_catbin.drop(['other_hispanic', 'other_eth', 'mexican', 'white'], axis=1)
no_eth_catbin_mexican = discovery_catbin.drop(['other_hispanic', 'other_eth', 'black', 'white'], axis=1)
no_eth_catbin_mexican2 = replication_catbin.drop(['other_hispanic', 'other_eth', 'black', 'white'], axis=1)
no_eth_catbin_white = discovery_catbin.drop(['other_hispanic', 'other_eth', 'mexican', 'black'], axis=1)
no_eth_catbin_white2 = replication_catbin.drop(['other_hispanic', 'other_eth', 'mexican', 'black'], axis=1)

####Merge columns for original analysis
discovery_all = clarite.modify.merge_variables(drop_white_dis, discovery_cont_phenotype, how="left")
replication_all = clarite.modify.merge_variables(drop_white_rep, replication_cont_phenotype, how="left")

####Set up for interactions 
discovery_all_exposures = clarite.modify.merge_variables(discovery_catbin, discovery_cont_phenotype, how="left")
discovery_all_exposures = discovery_all_exposures.drop(['white', 'black', 'mexican', "female", "SES_LEVEL", "RIDAGEYR", "SDDSRVYR", "BMXBMI"], axis=1)
discovery_exp_no_pheno = discovery_all_exposures.drop(pheno, axis=1)
discovery_exp_list = list(discovery_exp_no_pheno)

####Merge variables (both have 9,063 rows)
discovery_all_black = clarite.modify.merge_variables(no_eth_catbin_black, discovery_cont_phenotype, how="left")
replication_all_black = clarite.modify.merge_variables(no_eth_catbin_black2, replication_cont_phenotype, how="left")
discovery_all_mexican = clarite.modify.merge_variables(no_eth_catbin_mexican, discovery_cont_phenotype, how="left")
replication_all_mexican = clarite.modify.merge_variables(no_eth_catbin_mexican2, replication_cont_phenotype, how="left")
discovery_all_white = clarite.modify.merge_variables(no_eth_catbin_white, discovery_cont_phenotype, how="left")
replication_all_white = clarite.modify.merge_variables(no_eth_catbin_white2, replication_cont_phenotype, how="left")


# Get survey design for discovery and replication datasets
survey_design_discovery_match = survey_design_discovery.loc[discovery_all.index]
survey_design_replication_match = survey_design_replication.loc[replication_all.index]

# Get survey info only for the data
# Note: Ideally this doesn't happen (only subset is used to drop rows from the data)
survey_design_discovery_aa = survey_design_discovery.loc[discovery_all_black.index]
survey_design_replication_aa = survey_design_replication.loc[replication_all_black.index]
survey_design_discovery_mx = survey_design_discovery.loc[discovery_all_mexican.index]
survey_design_replication_mx = survey_design_replication.loc[replication_all_mexican.index]
survey_design_discovery_wh = survey_design_discovery.loc[discovery_all_white.index]
survey_design_replication_wh = survey_design_replication.loc[replication_all_white.index]

# Make survey design for original analysis
design_dis = clarite.survey.SurveyDesignSpec(survey_design_discovery_match,
                                         weights=weights_discovery, cluster="SDMVPSU", strata="SDMVSTRA", drop_unweighted= True,
                                         fpc=None, nest=True)

design_rep = clarite.survey.SurveyDesignSpec(survey_design_replication_match,
                                         weights=weights_replication, cluster="SDMVPSU", strata="SDMVSTRA", drop_unweighted= True,
                                         fpc=None, nest=True)

# Make survey design
design_aa = clarite.survey.SurveyDesignSpec(survey_design_discovery_aa,
                                         weights=weights_discovery, cluster="SDMVPSU", strata="SDMVSTRA", single_cluster= 'adjust', drop_unweighted= True,
                                         fpc=None, nest=True)

design2_aa = clarite.survey.SurveyDesignSpec(survey_design_replication_aa,
                                         weights=weights_replication, cluster="SDMVPSU", strata="SDMVSTRA", single_cluster= 'adjust', drop_unweighted= True,
                                         fpc=None, nest=True)

design_mx = clarite.survey.SurveyDesignSpec(survey_design_discovery_mx,
                                         weights=weights_discovery, cluster="SDMVPSU", strata="SDMVSTRA", single_cluster= 'adjust', drop_unweighted= True,
                                         fpc=None, nest=True)

design2_mx = clarite.survey.SurveyDesignSpec(survey_design_replication_mx,
                                         weights=weights_replication, cluster="SDMVPSU", strata="SDMVSTRA", single_cluster= 'adjust', drop_unweighted= True,
                                         fpc=None, nest=True)

design_wh = clarite.survey.SurveyDesignSpec(survey_design_discovery_wh,
                                         weights=weights_discovery, cluster="SDMVPSU", strata="SDMVSTRA", single_cluster= 'adjust', drop_unweighted= True,
                                         fpc=None, nest=True)

design2_wh = clarite.survey.SurveyDesignSpec(survey_design_replication_wh,
                                         weights=weights_replication, cluster="SDMVPSU", strata="SDMVSTRA", single_cluster= 'adjust', drop_unweighted= True,
                                         fpc=None, nest=True)

# Subset
design_aa.subset(discovery_all_black['black'] == 1)
design2_aa.subset(replication_all_black['black'] == 1)

design_mx.subset(discovery_all_mexican['mexican'] == 1)
design2_mx.subset(replication_all_mexican['mexican'] == 1)

design_wh.subset(discovery_all_white['white'] == 1)
design2_wh.subset(replication_all_white['white'] == 1)

# Drop 'black' since it was used to subset and is now a constant
discovery_all_black = discovery_all_black.drop(columns=['black'])
replication_all_black = replication_all_black.drop(columns=['black'])

discovery_all_mexican = discovery_all_mexican.drop(columns=['mexican'])
replication_all_mexican = replication_all_mexican.drop(columns=['mexican'])

discovery_all_white = discovery_all_white.drop(columns=['white'])
replication_all_white = replication_all_white.drop(columns=['white'])

# Run
results = []
for current_pheno in pheno:
    print(current_pheno)
    drop_phenos = [p for p in pheno if p != current_pheno]
    ewas_result_dis = clarite.analyze.ewas(outcome= current_pheno,
                                       covariates= cov,
                                       min_n=200,
                                       data=discovery_all.drop(columns=drop_phenos),
                                       survey_design_spec=design_dis)
    results.append(ewas_result_dis)
result = pd.concat(results)
clarite.analyze.add_corrected_pvalues(result)


result.to_csv("Discovery_Results_python.txt", sep="\t")

# Get a dictionary of phenotype : list of significant variables
significant_results = result[result['pvalue_fdr']<0.1].index.values
sr_dict = dict()
for var, phenotype in significant_results:
    if phenotype in sr_dict:
        sr_dict[phenotype].append(var)
    else:
        sr_dict[phenotype] = [var]

# Run EWAS
replication_results = []
for phenotype, variable_list in sr_dict.items():
    print(phenotype)
    data = replication_all[[phenotype] + cov + variable_list]
    ewas_result_rep = clarite.analyze.ewas(outcome= phenotype,
                                                 covariates= cov, 
                                                 min_n=200,
                                                 data=data,
                                                 survey_design_spec=design_rep)
    replication_results.append(ewas_result_rep)
result_rep = pd.concat(replication_results)
clarite.analyze.add_corrected_pvalues(result_rep)

result_rep.to_csv("Replication_Results_python.txt", sep="\t")


# Save combined reports
    # Combine results
ewas_keep_cols = ['pvalue', 'pvalue_bonferroni', 'pvalue_fdr', 'N', 'Beta']
combined = pd.merge(result[['Variable_type'] + ewas_keep_cols],
                        result_rep[ewas_keep_cols],
                        left_index=True, right_index=True, suffixes=("_disc", "_repl"))
    
# FDR < 0.1 in both
fdr_significant = combined.loc[(combined['pvalue_fdr_disc'] <= 0.1) & (combined['pvalue_fdr_repl'] <= 0.1),]
fdr_significant = fdr_significant.assign(m=fdr_significant[['pvalue_fdr_disc', 'pvalue_fdr_repl']].mean(axis=1))\
                                     .sort_values('m').drop('m', axis=1)
fdr_significant.to_csv("Significant_Results_FDR_0.1_original.txt", sep="\t")
print(f"{len(fdr_significant)} variables had FDR < 0.1 in both discovery and replication")
    
    # Bonferroni < 0.05 in both
bonf_significant05 = combined.loc[(combined['pvalue_bonferroni_disc'] <= 0.05) & (combined['pvalue_bonferroni_repl'] <= 0.05),]
bonf_significant05 = bonf_significant05.assign(m=fdr_significant[['pvalue_bonferroni_disc', 'pvalue_bonferroni_repl']].mean(axis=1))\
                                           .sort_values('m').drop('m', axis=1)
bonf_significant05.to_csv("Significant_Results_Bonferroni_0.05_original.txt", sep="\t")
print(f"{len(bonf_significant05)} variables had Bonferroni < 0.05 in both discovery and replication")
    
# Bonferroni < 0.01 in both
bonf_significant01 = combined.loc[(combined['pvalue_bonferroni_disc'] <= 0.01) & (combined['pvalue_bonferroni_repl'] <= 0.01),]
bonf_significant01 = bonf_significant01.assign(m=fdr_significant[['pvalue_bonferroni_disc', 'pvalue_bonferroni_repl']].mean(axis=1))\
                                           .sort_values('m').drop('m', axis=1)
bonf_significant01.to_csv("Significant_Results_Bonferroni_0.01_original.txt", sep="\t")
print(f"{len(bonf_significant01)} variables had Bonferroni < 0.01 in both discovery and replication")

# Manhattan Plots
data_categories = pd.read_csv(data_var_categories, sep="\t").set_index('Variable')
data_categories.columns = ['category']
data_categories = data_categories['category'].to_dict()

clarite.plot.manhattan({'Discovery': result, 'Replication': result_rep}, num_labeled= 0, fdr= 0.1,
                       categories=data_categories, title="Weighted PheEWAS Results", filename='ewas_plot_org.png',
                       figsize=(14, 10))
    
                       


# Interactions
# Create the "race" variable, starting with all samples set to "other"
discovery_all_exposures["race"] = "other"
# Update for each variable
discovery_all_exposures["race"].loc[discovery_all_exposures["white"]==1] = "white"
discovery_all_exposures["race"].loc[discovery_all_exposures["black"]==1] = "black"
discovery_all_exposures["race"].loc[discovery_all_exposures["mexican"]==1] = "mexican"
# Make it a categorical and drop the original columns
discovery_all_exposures["race"] = discovery_all_exposures["race"].astype('category')
discovery_all_exposures = discovery_all_exposures.drop(columns=["white", "black", "mexican"])

#get model for interaction
interactions = [tuple(line.strip().split("\t")) for line in open("NHANES_interaction_race_model.txt").readlines()]

results = []
for x in pheno:
    print(x)
    interaction_dis = clarite.analyze.interaction_test(data= discovery_all_exposures,
                                       covariates= covariates,
                                       min_n=200,
                                       outcome_variable= x,
                                       interactions= interactions,
                                       report_betas=True)
    results.append(interaction_dis)
result = pd.concat(results)
clarite.analyze.add_corrected_pvalues(result, pvalue='LRT_pvalue', groupby=["Term1", "Term2"])

result.to_csv("Discovery_Results_interaction.txt", sep="\t")
discovery_all_exposures.to_csv("discovery_all_exposures.txt", sep="\t")

