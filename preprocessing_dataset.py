import pandas as pd
import numpy as np
import re
from pathlib import Path


INPUT_FILE = Path("data/SSBU26_dataset_fix.xls.xlsx")
OUTPUT_DIR = Path("output")

OUTPUT_DIR.mkdir(exist_ok=True)


def clean_columns(df):
    """
    Odstrani medzery z nazvov stlpcov.
    """
    df.columns = [str(col).strip() for col in df.columns]
    return df


def parse_date(value):
    """
    Prevedie datum z Excelu alebo textu na normalny datum.
    V datasete su datumy miesane, preto to riesime takto.
    """
    if pd.isna(value):
        return pd.NaT

    if isinstance(value, (int, float)) and not isinstance(value, bool):
        try:
            return pd.to_datetime("1899-12-30") + pd.to_timedelta(int(value), unit="D")
        except Exception:
            return pd.NaT

    return pd.to_datetime(value, errors="coerce")


def clean_text_value(value):
    """
    Upravi textovy vysledok, aby tam neboli zbytocne medzery.
    """
    if pd.isna(value):
        return np.nan

    text = str(value).strip()
    text = re.sub(r"\s+", " ", text)
    return text


def extract_numeric_value(value):
    """
    Pokusi sa vytiahnut cislo z vysledku.
    Napriklad:
    5,4 -> 5.4
    < 2,90 -> 2.90
    hemolyza -> NaN
    """
    if pd.isna(value):
        return np.nan

    text = str(value).strip().lower()
    text = text.replace(",", ".")

    invalid_values = [
        "",
        "-",
        "negat�vne",
        "negativne",
        "nepo��tan�",
        "nepocitane",
        "hemol�za",
        "hemolyza"
    ]

    if text in invalid_values:
        return np.nan

    match = re.search(r"[-+]?\d*\.?\d+", text)

    if match:
        try:
            return float(match.group())
        except Exception:
            return np.nan

    return np.nan


def get_value_flag(value):
    """
    Oznaci hodnoty, ktore boli zadane ako < alebo >.
    """
    if pd.isna(value):
        return ""

    text = str(value).strip()

    if text.startswith("<"):
        return "<"

    if text.startswith(">"):
        return ">"

    return ""


def is_hfe_result(lab_name):
    """
    Zisti, ci ide o HFE geneticke vysetrenie.
    """
    if pd.isna(lab_name):
        return False

    return str(lab_name).strip().startswith("HFE")


def classify_hfe_result(value):
    """
    Zjednodusi geneticky vysledok na kategoriu.
    """
    if pd.isna(value):
        return ""

    text = str(value).lower()

    if "negat" in text or "�tandardn� alely" in text or "standardne alely" in text:
        return "negative"

    if "heterozygot" in text:
        return "heterozygote"

    if "homozygot" in text:
        return "homozygote"

    if "mutant" in text:
        return "positive"

    return "unknown"


def create_patient_id(date_value):
    """
    Patient_ID vytvorime z fiktivneho datumu narodenia.
    """
    if pd.isna(date_value):
        return np.nan

    return date_value.strftime("%Y-%m-%d")


def get_patient_hfe_status(group):
    """
    Vytvori celkovy HFE status pacienta.
    """
    hfe_statuses = set(group["HFE_status"].dropna().astype(str))

    if "homozygote" in hfe_statuses:
        return "homozygote"

    if "heterozygote" in hfe_statuses:
        return "heterozygote"

    if "positive" in hfe_statuses:
        return "positive"

    if "negative" in hfe_statuses:
        return "negative"

    return "not_tested"


def main():
    print("Nacitavam Excel subor...")

    df = pd.read_excel(INPUT_FILE)
    df = clean_columns(df)

    print("Povodny pocet riadkov:", len(df))
    print("Stlpce v datasete:")
    print(list(df.columns))

    required_columns = [
        "Cislo_Z",
        "Datum_O",
        "Cas_O",
        "Datum_Nar",
        "Rok_Nar",
        "Pohlavie",
        "Lab_Nazov",
        "Vysl_Value"
    ]

    missing_columns = [col for col in required_columns if col not in df.columns]

    if missing_columns:
        print("Chybaju tieto stlpce:")
        print(missing_columns)
        return

    print("Cistim datumy...")

    df["Datum_O_clean"] = df["Datum_O"].apply(parse_date)
    df["Datum_Nar_clean"] = df["Datum_Nar"].apply(parse_date)

    print("Cistim vysledky...")

    df["Vysl_text_clean"] = df["Vysl_Value"].apply(clean_text_value)
    df["Vysl_numeric"] = df["Vysl_Value"].apply(extract_numeric_value)
    df["Vysl_flag"] = df["Vysl_Value"].apply(get_value_flag)

    print("Vytvaram Patient_ID...")

    df["Patient_ID"] = df["Datum_Nar_clean"].apply(create_patient_id)

    print("Oznacujem HFE geneticke vysledky...")

    df["Is_HFE"] = df["Lab_Nazov"].apply(is_hfe_result)
    df["HFE_status"] = ""

    df.loc[df["Is_HFE"], "HFE_status"] = df.loc[df["Is_HFE"], "Vysl_Value"].apply(classify_hfe_result)

    print("Oznacujem typ vysledku...")

    df["Result_type"] = "text"
    df.loc[df["Vysl_numeric"].notna(), "Result_type"] = "numeric"
    df.loc[df["Is_HFE"], "Result_type"] = "genetic"
    df.loc[df["Vysl_Value"].isna(), "Result_type"] = "missing"

    print("Odstranujem riadky bez pacienta alebo nazvu vysetrenia...")

    before_drop = len(df)

    df = df[df["Patient_ID"].notna()].copy()
    df = df[df["Lab_Nazov"].notna()].copy()

    after_drop = len(df)

    print("Odstranenych riadkov:", before_drop - after_drop)

    print("Vytvaram profil pacienta...")

    profiles = []

    for patient_id, group in df.groupby("Patient_ID"):
        sex_values = group["Pohlavie"].dropna().unique()

        if len(sex_values) == 1:
            sex = sex_values[0]
        elif len(sex_values) > 1:
            sex = "mixed"
        else:
            sex = ""

        birth_date = group["Datum_Nar_clean"].dropna().iloc[0]

        profiles.append({
            "Patient_ID": patient_id,
            "Datum_Nar": birth_date,
            "Rok_Nar": group["Rok_Nar"].dropna().iloc[0] if group["Rok_Nar"].notna().any() else np.nan,
            "Pohlavie": sex,
            "Pocet_riadkov": len(group),
            "Pocet_vysetreni": group["Cislo_Z"].nunique(),
            "Pocet_lab_parametrov": group["Lab_Nazov"].nunique(),
            "Prvy_odber": group["Datum_O_clean"].min(),
            "Posledny_odber": group["Datum_O_clean"].max(),
            "Ma_HFE_vysledok": bool(group["Is_HFE"].any()),
            "HFE_status": get_patient_hfe_status(group)
        })

    patient_profiles = pd.DataFrame(profiles)

    print("Ukladam vysledne subory...")

    df.to_csv(
        OUTPUT_DIR / "processed_final_dataset.csv",
        index=False,
        encoding="utf-8-sig"
    )

    patient_profiles.to_csv(
        OUTPUT_DIR / "patient_profiles.csv",
        index=False,
        encoding="utf-8-sig"
    )

    print()
    print("Hotovo.")
    print("Vystupy su ulozene v priecinku output.")
    print()
    print("processed_final_dataset.csv:", len(df), "riadkov")
    print("patient_profiles.csv:", len(patient_profiles), "pacientov")
    print()
    print("Zakladna kontrola:")
    print("Pocet pacientov:", df["Patient_ID"].nunique())
    print("Pocet lab parametrov:", df["Lab_Nazov"].nunique())
    print("Pocet HFE riadkov:", df["Is_HFE"].sum())
    print("Pocet chybajucich vysledkov:", df["Vysl_Value"].isna().sum())


if __name__ == "__main__":
    main()