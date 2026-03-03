import os
import pandas as pd

columns = [
           "disease",	"dataset",	"gwas_n",	"study"
          ]

disease = [
           {'disease': 'IBD', 'dataset': 'ebi-a-GCST004131', 'gwas_n': 59957, "model":"predixcan_mashr_eqtl"},
           {'disease': 'LDL', 'dataset': 'ukb-d-30780_irnt', 'gwas_n': 343621, "model":"predixcan_mashr_eqtl"},
           {'disease': 'SBP', 'dataset': 'ukb-a-360', 'gwas_n': 317754, "model":"predixcan_mashr_eqtl"}
           ]

disease = [
           {'disease': 'Asthma', 'dataset': 'Asthma_Bothsex_inv_var_meta_GBMI_052021_nbbkgt1_GBMI', 'gwas_n': 1800785},
           {'disease': 'Asthma', 'dataset': 'GCST006911_Shrine_N_2018', 'gwas_n': 30810},
           {'disease': 'Asthma', 'dataset': 'GCST010042_Han_Y_2020',    'gwas_n': 303859},
           {'disease': 'COPD',   'dataset': 'GCST90016593_Kim_W_2020',  'gwas_n': 14590}
           ]


df = pd.DataFrame(disease)
holder = []
counter = 0
for study in os.listdir(expression_path):
    if study.endswith(".db"):
        counter += 1
        tmp = df.copy()
        study = study.split('.')[0]
        tmp["study"] = study
        holder.append(tmp)
print(counter)
meta = pd.concat(holder).sort_values(by=["disease", "study"])
meta['model'] = 'no_model'
# meta['model'] = "PrediXcan_enet"

# expression_path = "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/data/Prediction_models/PredixcanV8/mashr_eqtl"
expression_path = "/home/jottensmeier/BioAdHoc/R24_Data/DICE-LUNG-combined/PREDICTION_MODELS/DICE-LUNG-combined/Models/Transcriptome/PrediXcan"

# out = f"/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/1_workflow/config/meta_gtex_{counter}.tsv"
# out = "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/1_workflow/config/meta_gtex_3_DICE_GWAS.tsv"

# out = "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/1_workflow/config/meta_DICE-LUNG-comb-tissue_enet_DICE_GWAS.tsv"
out = "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/1_workflow/config/meta_DICE-LUNG-comb-tissue_ss_DICE_GWAS.tsv"
meta.to_csv(out, sep="\t", header=True, index=False)