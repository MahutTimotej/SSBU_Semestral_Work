# SSBU analýza dát

Tento projekt je prototyp interaktívnej webovej aplikácie postavenej na frameworku Shiny pre Python. Slúži na vizualizáciu a analýzu prepojených genetických výsledkov (HFE mutácie) a laboratórnych vyšetrení (železo, CRP, pečeňové testy).

Aplikácia je rozdelená na dve hlavné časti, ktoré sú dostupné cez horné záložky.

## Hlavné funkcie

### 1. Profil pacienta a anomálie
Táto záložka je určená na analýzu konkrétneho pacienta vybraného zo zoznamu.
* **Vyhľadávanie podľa ID:** Rýchly výber pacienta.
* **Základné údaje pacienta:** Zobrazenie pohlavia, dopočítaného veku, počtu odberov a nameraných hodnôt.
* **Genetický HFE status:** Agregovaný výsledok z dostupných genetických vyšetrení pacienta.
* **Detekcia klinických anomálií:** 
  * *Kritické riziko hemochromatózy* (ak je pacient HFE pozitívny/homozygot a hladina železa presiahne 27 umol/l).
  * *Zápalové skreslenie* (ak je železo nad 25 umol/l a zároveň CRP nad 10 mg/l, čo indikuje, že zvýšené železo môže byť sekundárnym prejavom zápalu).
* **História vyšetrení:** Kompletná tabuľka nameraných hodnôt s možnosťou filtrovania a zoradenia.
* **Graf porovnania s populáciou:** Boxplot distribúcie železa v populácii podľa HFE statusu, v ktorom sú červenou farbou vyznačené reálne hodnoty vybraného pacienta.

### 2. Skupinové vizualizácie a vzťahy
Táto časť slúži na analýzu celého datasetu a hľadanie trendov v populácii.
* **Jednorozmerná analýza:** Histogramy a boxploty vybraných biochemických parametrov s možnosťou filtrovania podľa pohlavia a genetického statusu.
* **Dvojrozmerná analýza:** Scatter plot na porovnanie vzťahu medzi dvoma vybranými parametrami (napr. ALT vs AST).
* **Korelačná matica:** Heatmapa zobrazujúca Spearmanovu koreláciu medzi kľúčovými biochemickými parametrami (železo, feritín, transferín, pečeňové testy).

## Požiadavky na spustenie

Pred spustením je potrebné mať nainštalovaný Python a nasledovné knižnice:

```bash
pip install shiny pandas matplotlib seaborn
```

## Štruktúra projektu

Projekt má nasledovné usporiadanie priečinkov a súborov:

```text
SEMESTRALKA/
 ├── data/
 │    └── SSBU26_dataset_fix.xls.xlsx       # Zdrojový Excel súbor s dátami
 ├── output/                                # Generované výstupy z analýz a preprocessingu
 │    ├── bod3_asociacie.csv
 │    ├── bod3_deskriptivne_statistiky.csv
 │    ├── bod3_genotypy_pacientov.csv
 │    ├── bod3_korelacie_lab.csv
 │    ├── bod3_sumar.txt
 │    ├── patient_profiles.csv
 │    └── processed_final_dataset.csv
 ├── analyza_anomalii.py                    # Analýza klinických anomálií a ukladania železa
 ├── analyza_bod3.py                        # Štatistické testy asociácií a korelácií
 ├── analyza_genetiky.py                    # Analýza a vizualizácia HFE genotypových vzorov
 ├── app.py                                 # Hlavný kód interaktívnej Shiny aplikácie
 ├── documentation.docx                     # Sprievodná dokumentácia k projektu
 ├── preprocessing_dataset.py               # Preprocessing a čistenie vstupných dát
 ├── README.md                              # Tento súbor s dokumentáciou
 └── sumar_bod3.py                          # Skript na generovanie textového sumáru k bodu 3
```

Aplikácia načítava vyčistené dáta priamo zo zložky `output/`, ktorú generuje skript `preprocessing_dataset.py`.

## Ako aplikáciu spustiť

Aplikáciu spustíš z príkazového riadku v koreňovom adresári projektu:

```bash
python -m shiny run app.py --reload
```

Po spustení sa v termináli zobrazí adresa (štandardne `http://127.0.0.1:8000`), ktorú stačí otvoriť v prehliadači. Parameter `--reload` zabezpečí, že pri úprave kódu v `app.py` sa aplikácia sama reštartuje a zmeny sa prejavia hneď po obnovení stránky.