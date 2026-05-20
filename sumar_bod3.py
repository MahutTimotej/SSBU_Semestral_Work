"""
Sumar najdolezitejsich vysledkov bodu 3.
Vytvori prehladny textovy report pre dokumentaciu.
 
Pouzitie:
    python sumar_bod3.py
"""
 
import pandas as pd
from pathlib import Path
 
OUTPUT_DIR = Path("output")
 
 
def main():
    asoc = pd.read_csv(OUTPUT_DIR / "bod3_asociacie.csv", encoding="utf-8-sig")
    desc = pd.read_csv(OUTPUT_DIR / "bod3_deskriptivne_statistiky.csv", encoding="utf-8-sig")
    geno = pd.read_csv(OUTPUT_DIR / "bod3_genotypy_pacientov.csv", encoding="utf-8-sig")
 
    lines = []
    lines.append("=" * 70)
    lines.append("BOD 3 - SUMAR NAJDOLEZITEJSICH VYSLEDKOV")
    lines.append("=" * 70)
    lines.append("")
 
    # 1) Rozdelenie genotypov
    lines.append("1) ROZDELENIE GENOTYPOV V POPULACII")
    lines.append("-" * 70)
    lines.append(f"Celkovy pocet pacientov s HFE testom: {len(geno)}")
    lines.append("")
    for col, full in [("HFE_H63D", "H63D (C187G)"),
                      ("HFE_S65C", "S65C (A193T)"),
                      ("HFE_C282Y", "C282Y (G845A)")]:
        counts = geno[col].value_counts(dropna=False)
        total = counts.sum()
        lines.append(f"  {full}:")
        for g in ["wt", "het", "hom"]:
            n = counts.get(g, 0)
            pct = 100 * n / total if total else 0
            lines.append(f"    {g:5s} {n:4d} ({pct:5.1f}%)")
        lines.append("")
 
    # 2) Vyznamne asociacie
    lines.append("2) VYZNAMNE ASOCIACIE (FDR < 0.05)")
    lines.append("-" * 70)
    sig = asoc[asoc["Vyznamne_FDR"]].copy()
    lines.append(f"Pocet vyznamnych asociacii: {len(sig)} z {len(asoc)} testovanych")
    lines.append("")
    lines.append("TOP 15 najsilnejsich:")
    lines.append(f"{'Marker':<14}{'Lab parameter':<45}{'p_FDR':>12}")
    for _, r in sig.head(15).iterrows():
        lines.append(f"{r['Marker']:<14}{r['Lab_parameter']:<45}{r['p_FDR']:>12.2e}")
    lines.append("")
 
    # 3) Klucove biologicke nalezy
    lines.append("3) KLUCOVE BIOLOGICKE NALEZY")
    lines.append("-" * 70)
    lines.append("")
 
    key_params = [
        ("Železo_S", "Zelezo v serume (umol/l)"),
        ("Feritin_S", "Feritin (ug/l)"),
        ("Transferín_S", "Transferin (g/l)"),
        ("ALT_S", "ALT (ukat/l)"),
        ("AST_S", "AST (ukat/l)"),
        ("GMT_S", "GMT (ukat/l)"),
        ("Bili-celkový_S", "Bilirubin celkovy (umol/l)"),
        ("Krvný obraz: Hemoglobín_HGB", "Hemoglobin (g/l)"),
    ]
 
    for param, label in key_params:
        sub = desc[(desc["Marker"] == "HFE_C282Y") & (desc["Lab_parameter"] == param)]
        if sub.empty:
            continue
        lines.append(f"  {label}")
        for _, r in sub.iterrows():
            lines.append(f"    {r['Genotyp']:5s} N={r['N']:4d}  "
                         f"median={r['Median']:8.2f}  "
                         f"IQR=[{r['Q1']:.2f}, {r['Q3']:.2f}]")
        # p-hodnota z asociacii
        p_row = asoc[(asoc["Marker"] == "HFE_C282Y") &
                     (asoc["Lab_parameter"] == param)]
        if not p_row.empty:
            p = p_row.iloc[0]["p_FDR"]
            lines.append(f"    -> p_FDR = {p:.3e}")
        lines.append("")
 

    text = "\n".join(lines)
    out = OUTPUT_DIR / "bod3_sumar.txt"
    out.write_text(text, encoding="utf-8")
    print(text)
    print(f"\nUlozene: {out}")
 
 
if __name__ == "__main__":
    main()
 