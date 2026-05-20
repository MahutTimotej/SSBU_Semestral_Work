import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

INPUT_FILE = Path("output/processed_final_dataset.csv")
OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(exist_ok=True)

# Minimalny pocet pacientov v skupine, aby sa robil test.
MIN_GROUP_SIZE = 5

# Tri HFE markery, ktore su v datasete.
HFE_MARKERS = {
    "HFE_H63D": "HFE C187G(H63D)",
    "HFE_S65C": "HFE A193T(S65C)",
    "HFE_C282Y": "HFE G845A(C282Y)",
}


# ---------------------------------------------------------------------
# 1) NACITANIE A KLASIFIKACIA GENOTYPOV
# ---------------------------------------------------------------------

def classify_genotype(text):
    """
    Z textoveho vysledku urci genotyp:
      - 'wt'    = wild-type (negativny)
      - 'het'   = heterozygot (standardna + mutantna alela)
      - 'hom'   = homozygot (len mutantne alely)
      - NaN     = nevyhodnotitelne
    """
    if pd.isna(text):
        return np.nan

    t = str(text).lower()

    if "negat" in t or "štandardné alely" in t or "standardne alely" in t:
        return "wt"
    if "štandardná a mutantná" in t or "standardna a mutantna" in t:
        return "het"
    if "len mutantné alely" in t or "len mutantne alely" in t:
        return "hom"

    return np.nan


def build_patient_genotypes(df):
    """
    Pre kazdeho pacienta zisti genotyp pre kazdy zo 3 HFE markerov.
    Ak ma pacient viacero HFE vysetreni, zoberie sa to najzavaznejsie
    (hom > het > wt) - v praxi by mali byt rovnake, ale poistka.
    """
    hfe = df[df["Is_HFE"]].copy()
    hfe["Genotype"] = hfe["Vysl_Value"].apply(classify_genotype)

    severity = {"wt": 0, "het": 1, "hom": 2}
    hfe["Severity"] = hfe["Genotype"].map(severity)

    rows = []
    for pid, g in hfe.groupby("Patient_ID"):
        row = {"Patient_ID": pid}
        for short, full in HFE_MARKERS.items():
            sub = g[g["Lab_Nazov"] == full]
            if sub.empty or sub["Severity"].isna().all():
                row[short] = np.nan
            else:
                worst = sub.loc[sub["Severity"].idxmax()]
                row[short] = worst["Genotype"]
        rows.append(row)

    geno = pd.DataFrame(rows)

    # Spolocna kategoria - ma pacient akukolvek HFE mutaciu?
    def overall(r):
        vals = [r["HFE_H63D"], r["HFE_S65C"], r["HFE_C282Y"]]
        if all(pd.isna(v) for v in vals):
            return np.nan
        if "hom" in vals:
            return "any_hom"
        if "het" in vals:
            return "any_het"
        return "all_wt"

    geno["HFE_overall"] = geno.apply(overall, axis=1)
    return geno


# ---------------------------------------------------------------------
# 2) PRIPRAVA LAB DAT (long -> wide po vysetreni)
# ---------------------------------------------------------------------

def prepare_lab_data(df, patient_geno):
    """
    Z long formatu vytvori wide tabulku, kde 1 riadok = 1 vysetrenie pacienta,
    stlpce = jednotlive laboratorne parametre (numericke hodnoty).
    Ku kazdemu vysetreniu sa pripoji genotyp pacienta.
    """
    lab = df[(~df["Is_HFE"]) & (df["Vysl_numeric"].notna())].copy()

    # Niektore parametre maju malo merani - vyhodime tie s < 30 hodnotami,
    # aby sme nerobili statistiku na neprezitelnych vzorkach.
    counts = lab["Lab_Nazov"].value_counts()
    keep = counts[counts >= 30].index
    lab = lab[lab["Lab_Nazov"].isin(keep)]

    # Pivot: jeden riadok = jedno vysetrenie (Cislo_Z + Patient_ID).
    # Ak je v ramci jedneho vysetrenia ten isty parameter viackrat,
    # zoberieme medianu, aby outlier-y nestrhli vysledok.
    wide = lab.pivot_table(
        index=["Patient_ID", "Cislo_Z"],
        columns="Lab_Nazov",
        values="Vysl_numeric",
        aggfunc="median",
    ).reset_index()

    wide = wide.merge(patient_geno, on="Patient_ID", how="left")
    return wide, list(keep)


# ---------------------------------------------------------------------
# 3) DESKRIPTIVNE STATISTIKY PODLA GENOTYPU
# ---------------------------------------------------------------------

def descriptive_stats(wide, lab_params, marker):
    """
    Pre kazdy lab parameter spocita medianu, IQR a pocet hodnot
    v jednotlivych skupinach genotypu.
    """
    rows = []
    for param in lab_params:
        if param not in wide.columns:
            continue
        for geno_val, sub in wide.groupby(marker):
            values = sub[param].dropna()
            if len(values) < MIN_GROUP_SIZE:
                continue
            rows.append({
                "Marker": marker,
                "Genotyp": geno_val,
                "Lab_parameter": param,
                "N": len(values),
                "Median": round(values.median(), 3),
                "Q1": round(values.quantile(0.25), 3),
                "Q3": round(values.quantile(0.75), 3),
                "Priemer": round(values.mean(), 3),
                "SD": round(values.std(), 3),
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------
# 4) STATISTICKE TESTY ASOCIACII
# ---------------------------------------------------------------------

def benjamini_hochberg(pvals):
    """
    FDR korekcia Benjamini-Hochberg.
    Ked robime desiatky testov naraz, surove p-hodnoty preceuju asociacie.
    """
    p = np.array(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / (np.arange(n) + 1)
    # Monotonia - kazda adjustovana p musi byt aspon taka ako nasledujuca.
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n)
    out[order] = adj
    return out


def test_associations(wide, lab_params, marker):
    """
    Pre kazdy lab parameter porovna jeho rozdelenie medzi skupinami genotypu.
    Pouziva Mann-Whitney (2 skupiny) alebo Kruskal-Wallis (3 skupiny),
    lebo lab hodnoty zriedka maju normalne rozdelenie.
    """
    results = []
    for param in lab_params:
        if param not in wide.columns:
            continue
        groups, labels = [], []
        for geno_val, sub in wide.groupby(marker):
            vals = sub[param].dropna().values
            if len(vals) >= MIN_GROUP_SIZE:
                groups.append(vals)
                labels.append(str(geno_val))

        if len(groups) < 2:
            continue

        try:
            if len(groups) == 2:
                stat, p = stats.mannwhitneyu(groups[0], groups[1], alternative="two-sided")
                test = "Mann-Whitney U"
            else:
                stat, p = stats.kruskal(*groups)
                test = "Kruskal-Wallis"
        except ValueError:
            continue

        medians = {f"median_{lab}": round(float(np.median(g)), 3)
                   for lab, g in zip(labels, groups)}
        sizes = {f"n_{lab}": int(len(g)) for lab, g in zip(labels, groups)}

        results.append({
            "Marker": marker,
            "Lab_parameter": param,
            "Test": test,
            "Statistika": round(float(stat), 3),
            "p_hodnota": float(p),
            **sizes,
            **medians,
        })

    if not results:
        return pd.DataFrame()

    df_res = pd.DataFrame(results)
    df_res["p_FDR"] = benjamini_hochberg(df_res["p_hodnota"].values)
    df_res["Vyznamne_FDR"] = df_res["p_FDR"] < 0.05
    df_res = df_res.sort_values("p_FDR").reset_index(drop=True)
    return df_res


# ---------------------------------------------------------------------
# 5) KORELACIE MEDZI LAB PARAMETRAMI (uzitocne pre body 4 a 5)
# ---------------------------------------------------------------------

def lab_correlations(wide, lab_params):
    """
    Spearman korelacie - robustne k outlier-om.
    Vyplivne hornu trojuholnikovu cast matice ako dlhu tabulku.
    """
    cols = [c for c in lab_params if c in wide.columns]
    rows = []
    for i, a in enumerate(cols):
        for b in cols[i + 1:]:
            pair = wide[[a, b]].dropna()
            if len(pair) < 30:
                continue
            rho, p = stats.spearmanr(pair[a], pair[b])
            if pd.isna(rho):
                continue
            rows.append({
                "Param_A": a,
                "Param_B": b,
                "N": len(pair),
                "Spearman_rho": round(float(rho), 3),
                "p_hodnota": float(p),
            })
    df_res = pd.DataFrame(rows)
    if df_res.empty:
        return df_res
    df_res["p_FDR"] = benjamini_hochberg(df_res["p_hodnota"].values)
    df_res["abs_rho"] = df_res["Spearman_rho"].abs()
    df_res = df_res.sort_values("abs_rho", ascending=False).drop(columns="abs_rho")
    return df_res.reset_index(drop=True)


# ---------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------

def main():
    print("Nacitavam predspracovany dataset...")
    df = pd.read_csv(INPUT_FILE, encoding="utf-8-sig", low_memory=False)
    print(f"  {len(df):,} riadkov, {df['Patient_ID'].nunique()} pacientov")

    print("\nKlasifikujem HFE genotypy pacientov...")
    geno = build_patient_genotypes(df)
    print("  Rozdelenie pre kazdy marker:")
    for m in HFE_MARKERS:
        print(f"    {m}:")
        print(geno[m].value_counts(dropna=False).to_string().replace("\n", "\n      "))

    geno.to_csv(OUTPUT_DIR / "bod3_genotypy_pacientov.csv",
                index=False, encoding="utf-8-sig")

    print("\nPripravujem lab data (long -> wide)...")
    wide, lab_params = prepare_lab_data(df, geno)
    print(f"  {len(wide):,} vysetreni, {len(lab_params)} laboratornych parametrov")

    print("\nPocitam deskriptivne statistiky podla genotypu...")
    desc_all = []
    for short in list(HFE_MARKERS) + ["HFE_overall"]:
        desc_all.append(descriptive_stats(wide, lab_params, short))
    desc = pd.concat(desc_all, ignore_index=True)
    desc.to_csv(OUTPUT_DIR / "bod3_deskriptivne_statistiky.csv",
                index=False, encoding="utf-8-sig")
    print(f"  {len(desc)} riadkov ulozenych")

    print("\nTestujem asociacie genotyp vs. laboratorne parametre...")
    asoc_all = []
    for short in list(HFE_MARKERS) + ["HFE_overall"]:
        res = test_associations(wide, lab_params, short)
        if not res.empty:
            asoc_all.append(res)
    asoc = pd.concat(asoc_all, ignore_index=True) if asoc_all else pd.DataFrame()
    asoc = asoc.sort_values("p_FDR").reset_index(drop=True)
    asoc.to_csv(OUTPUT_DIR / "bod3_asociacie.csv",
                index=False, encoding="utf-8-sig")
    print(f"  {len(asoc)} testov, {asoc['Vyznamne_FDR'].sum()} vyznamnych po FDR korekcii")

    if len(asoc) > 0:
        print("\n  TOP 10 najsilnejsich asociacii (najnizsie p_FDR):")
        cols = ["Marker", "Lab_parameter", "Test", "p_hodnota", "p_FDR"]
        print(asoc[cols].head(10).to_string(index=False))

    print("\nPocitam korelacie medzi laboratornymi parametrami (Spearman)...")
    corr = lab_correlations(wide, lab_params)
    corr.to_csv(OUTPUT_DIR / "bod3_korelacie_lab.csv",
                index=False, encoding="utf-8-sig")
    print(f"  {len(corr)} parov ulozenych")

    if len(corr) > 0:
        print("\n  TOP 10 najsilnejsich korelacii (|rho|):")
        print(corr.head(10).to_string(index=False))

    print("\nHotovo. Vystupy:")
    for name in ["bod3_genotypy_pacientov.csv",
                 "bod3_deskriptivne_statistiky.csv",
                 "bod3_asociacie.csv",
                 "bod3_korelacie_lab.csv"]:
        print(f"  output/{name}")


if __name__ == "__main__":
    main()