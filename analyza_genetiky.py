import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

script_dir = os.path.dirname(os.path.abspath(__file__))
df = pd.read_csv(os.path.join(script_dir, 'processed_final_dataset.csv'), low_memory=False)

gen_df = df[df['Result_type'] == 'genetic'].pivot_table(
    index=['Patient_ID', 'Pohlavie', 'Rok_Nar'], 
    columns='Lab_Nazov', 
    values='HFE_status', 
    aggfunc='first'
).reset_index()

rename_dict = {'HFE C187G(H63D)': 'H63D', 'HFE A193T(S65C)': 'S65C', 'HFE G845A(C282Y)': 'C282Y'}
gen_df = gen_df.rename(columns=rename_dict)

# Vytvorenie Genotype Pattern
mutation_cols = ['H63D', 'S65C', 'C282Y']
gen_df['Genotype_Pattern'] = gen_df[mutation_cols].apply(
    lambda row: ' | '.join([f"{col}:{row[col]}" for col in mutation_cols]), axis=1
)

# 2. Analýza vzorov
patterns = gen_df['Genotype_Pattern'].value_counts()
compound_het = gen_df[(gen_df['H63D'] == 'positive') & (gen_df['C282Y'] == 'positive')]
rare_threshold = len(gen_df) * 0.01
rare_cases = gen_df[gen_df['Genotype_Pattern'].map(patterns) < rare_threshold]

print(f"--- GENETICKÁ ANALÝZA ---")
print(f"Celkový počet pacientov s genetikou: {len(gen_df)}")
print(f"Počet compound heterozygótov (H63D+C282Y): {len(compound_het)}")
print(f"Počet vzácnych kombinácií (<1%): {len(rare_cases)}")

plt.figure(figsize=(14, 10))

# Graf 1: Distribúcia vzorov
plt.subplot(2, 1, 1)
patterns.plot(kind='barh', color='skyblue')
plt.title('Distribúcia genetických vzorov (HFE kombinácie)')
plt.xlabel('Počet pacientov')

# Graf 2: Veková distribúcia nosičov mutácie C282Y
plt.subplot(2, 1, 2)
gen_df['Vek'] = 2024 - gen_df['Rok_Nar']
sns.histplot(data=gen_df, x='Vek', hue='C282Y', multiple="stack", palette='magma')
plt.title('Veková štruktúra pacientov podľa mutácie C282Y')

plt.tight_layout()
plt.savefig(os.path.join(script_dir, 'graf_geneticka_distribucia.png'), dpi=300)
print(f"Obrázok uložený: graf_geneticka_distribucia.png")
plt.show()

rare_cases.to_csv(os.path.join(script_dir, 'vzacne_geneticke_pripady.csv'), index=False)