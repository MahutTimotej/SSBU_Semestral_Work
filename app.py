from shiny import App, ui, render, reactive
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / "output"

try:
    df_lab = pd.read_csv(OUTPUT_DIR / "processed_final_dataset.csv", low_memory=False)
    df_profiles = pd.read_csv(OUTPUT_DIR / "patient_profiles.csv", low_memory=False)
    df_lab['Vysl_numeric'] = pd.to_numeric(df_lab['Vysl_numeric'], errors='coerce')
    patient_ids = df_profiles['Patient_ID'].dropna().unique().tolist()
except FileNotFoundError:
    print(f"CHYBA: Nenašiel sa priečinok '{OUTPUT_DIR}' alebo potrebné CSV súbory.")
    df_lab = pd.DataFrame()
    df_profiles = pd.DataFrame()
    patient_ids = []
except Exception as e:
    print(f"Neočakávaná chyba pri načítaní dát: {e}")
    df_lab = pd.DataFrame()
    df_profiles = pd.DataFrame()
    patient_ids = []

app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.h3("Vyhľadávanie"),
        ui.input_selectize(
            "patient_id", 
            "Zadaj ID pacienta:", 
            choices=patient_ids,
            multiple=False
        ),
        ui.hr()
    ),
    
    ui.h2("Profil pacienta"),
    
    ui.output_ui("alerts"),
    
    ui.navset_card_underline(
        # Záložka 1: Profil a Dáta
        ui.nav_panel(
            "Profil a Laboratórne výsledky",
            ui.layout_columns(
                ui.card(
                    ui.card_header("Základné informácie"),
                    ui.output_ui("patient_info")
                ),
                ui.card(
                    ui.card_header("HFE Genetický status"),
                    ui.output_ui("genetic_info")
                )
            ),
            ui.card(
                ui.card_header("História vyšetrení"),
                ui.output_data_frame("lab_table")
            )
        ),
        
        ui.nav_panel(
            "Vizualizácia v populácii",
            ui.card(
                ui.card_header("Hladina železa pacienta v porovnaní s ostatnými (podľa HFE statusu)"),
                ui.output_plot("iron_plot")
            )
        )
    )
)

def server(input, output, session):

    @reactive.Calc
    def get_patient_profile():
        pid = input.patient_id()
        if not pid or df_profiles.empty:
            return None
        # Vráti prvý riadok patriaci pacientovi (vždy je len jeden)
        match = df_profiles[df_profiles['Patient_ID'] == pid]
        return match.iloc[0] if not match.empty else None

    @reactive.Calc
    def get_patient_labs():
        pid = input.patient_id()
        if not pid or df_lab.empty:
            return pd.DataFrame()
        return df_lab[df_lab['Patient_ID'] == pid]

    @output
    @render.ui
    def patient_info():
        prof = get_patient_profile()
        if prof is None:
            return ui.p("Vyberte pacienta pre zobrazenie detailov.")
        
        vek = 2026 - prof['Rok_Nar'] if pd.notna(prof['Rok_Nar']) else "Neznámy"
        
        return ui.HTML(f"""
            <ul style="font-size: 1.1em; line-height: 1.6;">
                <li><b>Pohlavie:</b> {prof.get('Pohlavie', 'Neznáme')}</li>
                <li><b>Rok narodenia:</b> {prof.get('Rok_Nar', 'Neznámy')} (Vek: ~{vek} r.)</li>
                <li><b>Celkový počet vyšetrení:</b> {prof.get('Pocet_vysetreni', 0)}</li>
                <li><b>Počet nameraných parametrov:</b> {prof.get('Pocet_riadkov', 0)} (počet riadkov)</li>
                <li><b>Prvý odber:</b> {prof.get('Prvy_odber', 'Chýba')}</li>
                <li><b>Posledný odber:</b> {prof.get('Posledny_odber', 'Chýba')}</li>
            </ul>
        """)

    # Vykreslenie genetiky
    @output
    @render.ui
    def genetic_info():
        prof = get_patient_profile()
        if prof is None:
            return ui.p("Vyberte pacienta.")
        
        status = str(prof.get('HFE_status', 'not_tested'))
        color = "#d9534f" if status in ['homozygote', 'positive', 'heterozygote'] else "#5cb85c"
        if status == 'not_tested': color = "#777"
        
        return ui.HTML(f"""
            <h3 style='color:{color}; text-transform: uppercase;'>{status}</h3>
        """)

    @output
    @render.ui
    def alerts():
        labs = get_patient_labs()
        prof = get_patient_profile()
        if labs.empty or prof is None:
            return ui.HTML("")
            
        alerts_html = ""
        
        fe_vals = labs[labs['Lab_Nazov'].str.contains('Železo', na=False, case=False)]['Vysl_numeric']
        crp_vals = labs[labs['Lab_Nazov'].str.contains('CRP', na=False, case=False)]['Vysl_numeric']
        
        max_fe = fe_vals.max() if not fe_vals.empty else 0
        max_crp = crp_vals.max() if not crp_vals.empty else 0
        is_hfe_pos = prof.get('HFE_status', '') in ['homozygote', 'positive']
        
        if max_fe > 27 and is_hfe_pos:
            alerts_html += f"""
            <div style="background-color:#ffe6e6; padding:15px; border-left: 6px solid #d9534f; border-radius: 4px; margin-bottom:15px;">
                <h4 style="margin-top:0; color:#d9534f;">Kritické upozornenie: Genetická anomália (Hemochromatóza)</h4>
                Pacient má potvrdenú mutáciu a kriticky vysokú hladinu železa (<b>{max_fe} μmol/l</b>). Riziko poškodenia orgánov.
            </div>
            """
            
        if max_fe > 25 and max_crp > 10:
            alerts_html += f"""
            <div style="background-color:#fff8e6; padding:15px; border-left: 6px solid #f0ad4e; border-radius: 4px; margin-bottom:15px;">
                <h4 style="margin-top:0; color:#f0ad4e;">Upozornenie: Zápalové skreslenie hodnôt</h4>
                Vysoká hladina železa (<b>{max_fe} μmol/l</b>) môže byť sekundárnym prejavom zápalu, nakoľko pacient má zvýšené CRP (<b>{max_crp} mg/l</b>). 
            </div>
            """
            
        return ui.HTML(alerts_html)

    @output
    @render.data_frame
    def lab_table():
        labs = get_patient_labs()
        if labs.empty:
            return pd.DataFrame()
        # Vrátime len relevantné stĺpce a premenujeme ich pre lepšiu čitateľnosť
        show_df = labs[['Datum_O_clean', 'Lab_Nazov', 'Vysl_numeric', 'Vysl_text_clean', 'Result_type']].copy()
        show_df.columns = ['Dátum odberu', 'Názov vyšetrenia', 'Numerický výsledok', 'Textový výsledok', 'Typ']
        show_df = show_df.sort_values(by='Dátum odberu', ascending=False)
        return render.DataTable(show_df, filters=True)

    @output
    @render.plot
    def iron_plot():
        labs = get_patient_labs()
        pid = input.patient_id()
        if df_lab.empty or df_profiles.empty:
            return None

        fe_pop = df_lab[df_lab['Lab_Nazov'] == 'Železo_S'].dropna(subset=['Vysl_numeric']).copy()
        
        if 'HFE_status' in fe_pop.columns:
            fe_pop = fe_pop.drop(columns=['HFE_status'])
        fe_pop = fe_pop.merge(df_profiles[['Patient_ID', 'HFE_status']], on='Patient_ID', how='inner')
        
        fe_pop = fe_pop[fe_pop['HFE_status'].isin(['homozygote', 'heterozygote', 'positive', 'negative'])]

        sns.set_theme(style="whitegrid")
        fig, ax = plt.subplots(figsize=(10, 6))
        
        sns.boxplot(data=fe_pop, x='HFE_status', y='Vysl_numeric', ax=ax, color='lightblue')
        
        ax.set_title("Hodnoty Železa: Populácia vs. Vybraný pacient", fontsize=14, pad=15)
        ax.set_ylabel("Železo v sére (μmol/l)", fontsize=11)
        ax.set_xlabel("Genetický HFE status pacienta", fontsize=11)
        ax.axhline(27, color='#d9534f', linestyle='--', linewidth=2, label='Kritická hranica (27 μmol/l)')
        
        patient_fe = labs[labs['Lab_Nazov'] == 'Železo_S']['Vysl_numeric'].dropna()
        if not patient_fe.empty:
            prof = get_patient_profile()
            p_status = prof.get('HFE_status', 'not_tested') if prof is not None else 'not_tested'
            
            if p_status in fe_pop['HFE_status'].unique():
                ax.scatter([p_status]*len(patient_fe), patient_fe, 
                           color='#333333', edgecolor='white', s=120, zorder=10, 
                           label=f'Merania pacienta ({pid})')
            
        ax.legend(loc='upper right')
        plt.tight_layout()
        return fig

app = App(app_ui, server)
