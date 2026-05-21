import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
df = pd.read_csv(os.path.join(script_dir, 'processed_final_dataset.csv'), low_memory=False)
df_summary = pd.read_csv(os.path.join(script_dir, 'patient_profiles.csv'), low_memory=False)


pivot_num = df.pivot_table(index='Patient_ID', columns='Lab_Nazov', values='Vysl_numeric', aggfunc='mean')
pivot_gen = df[df['Result_type'] == 'genetic'].pivot_table(index='Patient_ID', columns='Lab_Nazov', values='HFE_status', aggfunc='first')

final_df = pd.merge(pivot_num, pivot_gen, on='Patient_ID', how='inner').reset_index()
final_df = pd.merge(final_df, df_summary[['Patient_ID', 'Pohlavie', 'Rok_Nar']], on='Patient_ID', how='inner')

iron_col, crp_col, alt_col, c282y_col = 'Železo_S', 'CRP_S', 'ALT_S', 'HFE G845A(C282Y)'
needed_cols = [iron_col, crp_col, alt_col, 'Feritin_S', 'Kreatinín_S', 'Kys. močová_S']
final_df[needed_cols] = final_df[needed_cols].apply(pd.to_numeric, errors='coerce')

# Vypocet veku
final_df['Vek'] = 2024 - final_df['Rok_Nar']

# Detekcia anomálií (Z-Score > 3)
def get_outliers(df, col):
    return df[df[col] > (df[col].mean() + 3 * df[col].std())]

print("--- ANALÝZA KLINICKÝCH ANOMÁLIÍ ---")

# Štatistika kumulácie železa (Vek/Pohlavie) pre nosičov mutácie
positives = final_df[final_df[c282y_col] == 'positive']
m_avg = positives[positives['Pohlavie'] == 'M'][iron_col].mean()
f_under_50 = positives[(positives['Pohlavie'] == 'F') & (positives['Vek'] < 50)][iron_col].mean()
f_over_50 = positives[(positives['Pohlavie'] == 'F') & (positives['Vek'] >= 50)][iron_col].mean()

print(f"Priemerné železo u mužov (C282Y+): {m_avg:.2f}")
print(f"Priemerné železo u žien pod 50r (C282Y+): {f_under_50:.2f}")
print(f"Priemerné železo u žien nad 50r (C282Y+): {f_over_50:.2f}")

# Riziková hemochromatóza
risk_group = final_df[(final_df[iron_col] > 27) & (final_df[c282y_col] == 'positive')]
print(f"Vysoké riziko (Mutácia + Fe > 27): {len(risk_group)} pacientov")

# Zápalové skreslenie
inflammation = final_df[(final_df[iron_col] > 25) & (final_df[crp_col] > 10)]
print(f"Zápalové skreslenie (Vysoké Fe kvôli CRP): {len(inflammation)} pacientov")

# 3. Vizualizácie
fig, axes = plt.subplots(3, 2, figsize=(16, 18))

# Korelačná matica 
cols_corr = [c for c in needed_cols if c in final_df.columns]
sns.heatmap(final_df[cols_corr].corr(), annot=True, cmap='coolwarm', ax=axes[0,0])
axes[0,0].set_title('Korelácie biochemických parametrov')

# Boxplot - Vplyv mutácie na Železo
sns.boxplot(x=c282y_col, y=iron_col, data=final_df, palette='Set2', ax=axes[0,1])
axes[0,1].set_title('Hladina Železa podľa C282Y statusu')

# Scatter plot Železo vs CRP (hľadanie anomálií)
sns.scatterplot(data=final_df, x=crp_col, y=iron_col, hue=c282y_col, alpha=0.6, ax=axes[1,0])
axes[1,0].axvline(10, color='r', linestyle='--')
axes[1,0].axhline(27, color='g', linestyle='--')
axes[1,0].set_title('Identifikácia zápalových vs genetických anomálií')

# ALT (Pečeň) vs Železo
sns.regplot(data=final_df, x=iron_col, y=alt_col, scatter_kws={'alpha':0.3}, line_kws={'color':'red'}, ax=axes[1,1])
axes[1,1].set_title('Vzťah medzi Železom a poškodením pečene (ALT)')

# GRAF kumulacie ŽELEZA (Vek vs Pohlavie)
sns.regplot(data=positives[positives['Pohlavie']=='M'], x='Vek', y=iron_col, label='Muži', color='blue', scatter_kws={'alpha':0.3}, ax=axes[2,0])
sns.regplot(data=positives[positives['Pohlavie']=='F'], x='Vek', y=iron_col, label='Ženy', color='red', scatter_kws={'alpha':0.3}, ax=axes[2,0])
axes[2,0].axvline(50, color='black', linestyle=':', label='Menopauza')
axes[2,0].set_title('Kumulácia železa u nosičov mutácie (Vek a Pohlavie)')
axes[2,0].legend()

sns.histplot(data=final_df, x='Vek', hue='Pohlavie', multiple="stack", ax=axes[2,1], palette='pastel')
axes[2,1].set_title('Veková distribúcia pacientov v datasete')

plt.tight_layout()
plt.savefig(os.path.join(script_dir, 'graf_klinicke_anomalie.png'), dpi=300)
print(f"Obrázok uložený: graf_klinicke_anomalie.png")
plt.show()

final_df.to_csv(os.path.join(script_dir, 'kompletna_analyza_pacientov.csv'), index=False, encoding='utf-8-sig')