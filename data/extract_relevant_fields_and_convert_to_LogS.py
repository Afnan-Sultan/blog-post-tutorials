import pandas as pd
import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import Descriptors

# Setup output directories
output_dir_full = "aqueous_solubility/aggregated"
output_dir_sources = "aqueous_solubility/by_source"
os.makedirs(output_dir_full, exist_ok=True)
os.makedirs(output_dir_sources, exist_ok=True)

min_num_molecules = 1000


# ---------- Utility Functions ----------
def get_molwt(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Descriptors.MolWt(mol) if mol else np.nan


def convert_log_ug_per_ml_to_logS(log_ug_per_ml, smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return np.nan
        mw = Descriptors.MolWt(mol)
        if pd.isna(mw) or mw <= 0:
            return np.nan
        ug_per_ml = 10 ** float(log_ug_per_ml)
        mol_l = (ug_per_ml * 1e-6) / mw
        return np.log10(mol_l) if mol_l > 0 else np.nan
    except:
        return np.nan


def group_and_save_by_condition(df, group_cols, smiles_col, logs_col, dataset_name, output_dir, min_count=1000):
    counts = df.groupby(group_cols, observed=True)[smiles_col].nunique().sort_values(ascending=False)
    counts_filtered = counts[counts > min_count]

    print(f"🔢 {dataset_name} contains {counts.index.droplevel(list(range(len(group_cols) - 1))).nunique()} unique source group(s).")
    print(f"📋 Groups with >{min_count} unique solutes:")
    if counts_filtered.empty:
        print("⚠️ No groups found.\n")
        return

    for keys, count in counts_filtered.items():
        keys = keys if isinstance(keys, tuple) else (keys,)
        key_str = " | ".join([f"{k:.2f}°C" if isinstance(k, float) else str(k) for k in keys])
        print(f"  {key_str}: {count} molecules")

        mask = pd.Series(True, index=df.index)
        for col, val in zip(group_cols, keys):
            mask &= df[col] == val

        df_filtered = df.loc[mask, [smiles_col, logs_col]].rename(
            columns={smiles_col: "SMILES", logs_col: "LogS"}
        ).dropna()

        filename_parts = [dataset_name] + [str(round(k)) if isinstance(k, float) else str(k)[:50] for k in keys]
        filename = "_".join(filename_parts).replace("/", "_") + ".csv"
        df_filtered.to_csv(os.path.join(output_dir, filename), index=False)
    print('\n')

# ---------- 1. AqSolDB ----------
df_aqsoldb = pd.read_csv('raw_datasets/AqSolDB.csv')
df_aqsoldb['Source'] = df_aqsoldb['ID'].astype(str).str.extract(r'^([A-Z])')
df_aqsoldb = df_aqsoldb.dropna(subset=['Source'])
df_aqsoldb['LogS'] = df_aqsoldb['Solubility']
df_aqsoldb_filtered = df_aqsoldb[['SMILES', 'LogS', 'Source']]
df_aqsoldb_filtered.to_csv(f'{output_dir_full}/AqSolDB.csv', index=False)

group_and_save_by_condition(
    df=df_aqsoldb_filtered,
    group_cols=['Source'],
    smiles_col='SMILES',
    logs_col='LogS',
    dataset_name="AqSolDB",
    output_dir=output_dir_sources,
    min_count=min_num_molecules
)

# ---------- 2. ESol ----------
df_esol = pd.read_csv('raw_datasets/ESol.csv')
df_esol = df_esol[['SMILES', 'measured log(solubility:mol/L)']].copy()
df_esol.columns = ['SMILES', 'LogS']

# Calculate MW and assign MW-based source groups
df_esol['MW'] = df_esol['SMILES'].apply(get_molwt)
df_esol['MW_Group'] = pd.cut(
    df_esol['MW'],
    bins=[-np.inf, 200, 300, np.inf],
    labels=["<200", "200-300", ">400"]
)

df_esol.to_csv(f'{output_dir_full}/ESol.csv', index=False)
df_esol.to_csv(f'{output_dir_sources}/ESol.csv', index=False)

# Group and save by MW-based pseudo-source
group_and_save_by_condition(
    df=df_esol,
    group_cols=['MW_Group'],
    smiles_col='SMILES',
    logs_col='LogS',
    dataset_name="ESol",
    output_dir=output_dir_sources,
    min_count=min_num_molecules
)

# ---------- 3. AZ Solubility ----------
df_az = pd.read_csv('raw_datasets/solubility_az.csv', sep=';')
df_az_filtered = df_az[['Smiles', 'Standard Value']].copy()
df_az_filtered['Standard Value'] = df_az_filtered['Standard Value'].astype(float)
df_az_filtered['LogS'] = df_az_filtered['Standard Value'].apply(lambda x: np.log10(x * 1e-9))
df_az_filtered = df_az_filtered.rename(columns={'Smiles': 'SMILES'})[['SMILES', 'LogS']]
df_az_filtered.to_csv(f'{output_dir_sources}/solubility_az.csv', index=False)

# ---------- 4. Biogen (Fang et al.) ----------
df_fang = pd.read_csv('raw_datasets/solubility_biogen.csv')
df_fang = df_fang.dropna(subset=['SMILES', 'LOG SOLUBILITY PH 6.8 (ug/mL)'])
df_fang['LogS'] = df_fang.apply(
    lambda row: convert_log_ug_per_ml_to_logS(row['LOG SOLUBILITY PH 6.8 (ug/mL)'], row['SMILES']), axis=1)
df_fang_filtered = df_fang[['SMILES', 'LogS']]
df_fang_filtered.to_csv(f'{output_dir_sources}/solubility_biogen.csv', index=False)

# ---------- 5. BigSolDB ----------
df_bigsol = pd.read_csv('raw_datasets/BigSolDBv2.0.csv')
df_bigsol = df_bigsol[df_bigsol['Solvent'].str.lower() == 'water']
df_bigsol['Temperature_C'] = df_bigsol['Temperature_K'] - 273.15

group_and_save_by_condition(
    df=df_bigsol,
    group_cols=['Temperature_C', 'Source'],
    smiles_col='SMILES_Solute',
    logs_col='LogS(mol/L)',
    dataset_name="BigSolDB",
    output_dir=output_dir_sources,
    min_count=min_num_molecules
)

df_bigsol_out = df_bigsol[['SMILES_Solute', 'LogS(mol/L)']].rename(
    columns={'SMILES_Solute': 'SMILES', 'LogS(mol/L)': 'LogS'})
df_bigsol_out.to_csv(f"{output_dir_full}/BigSolDB.csv", index=False)

# ---------- 6. SOMAS ----------
df_somas = pd.read_csv('raw_datasets/SOMAS.csv')
df_somas = df_somas.dropna(subset=['SMILES', 'Experimental Solubility in Water', 'Temperature'])

# Compute MW for conversion from mg/L to mol/L
df_somas['MW'] = df_somas['SMILES'].apply(get_molwt)

# Convert temperature to Celsius
df_somas['Temperature_C'] = df_somas['Temperature'] - 273.15

# Convert mg/L → mol/L → LogS
def convert_mg_per_l_to_logs(mg_per_l, mw):
    try:
        if pd.isna(mg_per_l) or pd.isna(mw) or mw <= 0:
            return np.nan
        mol_l = (mg_per_l * 1e-3) / mw  # mg → g → mol
        return np.log10(mol_l) if mol_l > 0 else np.nan
    except:
        return np.nan

df_somas['LogS'] = df_somas.apply(
    lambda row: convert_mg_per_l_to_logs(row['Experimental Solubility in Water'], row['MW']), axis=1
)

# Group and save
group_and_save_by_condition(
    df=df_somas,
    group_cols=['Temperature_C', 'Experiment Reference'],
    smiles_col='SMILES',
    logs_col='LogS',
    dataset_name="SOMAS",
    output_dir=output_dir_sources,
    min_count=min_num_molecules
)

# Save full file
df_somas[['SMILES', 'LogS']].dropna().to_csv(f"{output_dir_full}/SOMAS.csv", index=False)
