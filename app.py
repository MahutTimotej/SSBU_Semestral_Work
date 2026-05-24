from shiny import App, ui, render, reactive
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

#load ds
BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / "output"

try:
    df = pd.read_csv(OUTPUT_DIR / "processed_final_dataset.csv", low_memory=False)
    profiles = pd.read_csv(OUTPUT_DIR / "patient_profiles.csv", low_memory=False)
    
    df['Vysl_numeric'] = pd.to_numeric(df['Vysl_numeric'], errors='coerce')
    
    profiles_hfe = profiles[["Patient_ID", "HFE_status"]].copy()
    profiles_hfe = profiles_hfe.rename(columns={"HFE_status": "HFE_status_patient"})
    if 'HFE_status_patient' not in df.columns:
        df = df.merge(profiles_hfe, on="Patient_ID", how="left")
    
    patient_ids = profiles['Patient_ID'].dropna().unique().tolist()
    lab_test_count = df["Lab_Nazov"].dropna().nunique() if "Lab_Nazov" in df.columns else 0
    genetic_test_count = 0
    genetic_name_cols = [c for c in profiles.columns if "genet" in c.lower() and "status" not in c.lower()]
    if genetic_name_cols:
        genetic_test_count = profiles[genetic_name_cols[0]].dropna().nunique()
    elif "HFE_status" in profiles.columns:
        genetic_test_count = 1
except FileNotFoundError:
    print(f"CHYBA: Nenašiel sa priečinok '{OUTPUT_DIR}' alebo potrebné CSV súbory.")
    df = pd.DataFrame()
    profiles = pd.DataFrame()
    patient_ids = []
    lab_test_count = 0
    genetic_test_count = 0


selected_parameters = [
    "Kreatinín_S", "Urea_S", "AST_S", "ALT_S", "CRP_S", "GMT_S", "ALP_S",
    "Kálium-S", "Nátrium-S", "Chloridy-S", "Glukóza_ S", "Albumín_S",
    "Bielkoviny_S", "Bili-celkový_S", "Bili-priamy_S", "Železo_S", "Feritin_S",
    "Transferín_S", "Saturácia Trf"
]

available_parameters = [p for p in selected_parameters if p in df["Lab_Nazov"].unique()]

heatmap_parameters = [
    p for p in ["Kreatinín_S", "Urea_S", "AST_S", "ALT_S", "CRP_S", "Glukóza_ S", 
                "Železo_S", "Feritin_S", "Transferín_S", "Saturácia Trf"]
    if p in df["Lab_Nazov"].unique()
]

PARAM_UNITS = {
    "Kreatinín_S": "µmol/l",
    "Urea_S": "mmol/l",
    "AST_S": "U/l",
    "ALT_S": "U/l",
    "CRP_S": "mg/l",
    "GMT_S": "U/l",
    "ALP_S": "U/l",
    "Kálium-S": "mmol/l",
    "Nátrium-S": "mmol/l",
    "Chloridy-S": "mmol/l",
    "Glukóza_ S": "mmol/l",
    "Albumín_S": "g/l",
    "Bielkoviny_S": "g/l",
    "Bili-celkový_S": "µmol/l",
    "Bili-priamy_S": "µmol/l",
    "Železo_S": "µmol/l",
    "Feritin_S": "µg/l",
    "Transferín_S": "g/l",
    "Saturácia Trf": "%",
}

def format_param_label(param):
    unit = PARAM_UNITS.get(param)
    return f"{param} ({unit})" if unit else param

# filter data
def filter_data(parameter, gender_filter, hfe_filter):
    data = df[(df["Lab_Nazov"] == parameter) & (df["Vysl_numeric"].notna())].copy()
    if gender_filter != "Všetci": data = data[data["Pohlavie"] == gender_filter]
    if hfe_filter != "Všetci": data = data[data["HFE_status_patient"] == hfe_filter]
    return data

def filter_long_data(parameters, gender_filter, hfe_filter):
    data = df[(df["Lab_Nazov"].isin(parameters)) & (df["Vysl_numeric"].notna())].copy()
    if gender_filter != "Všetci": data = data[data["Pohlavie"] == gender_filter]
    if hfe_filter != "Všetci": data = data[data["HFE_status_patient"] == hfe_filter]
    return data

app_ui = ui.page_navbar(
    #patient detailt
    ui.nav_panel(
        "Profil Pacienta a Anomálie",
        ui.layout_sidebar(
            ui.sidebar(
                ui.h4("Vyhľadávanie"),
                ui.input_selectize("patient_id", "Zadaj ID pacienta:", choices=patient_ids, multiple=False),
                ui.hr(),
                width=300
            ),
            
            ui.h3("Biomedicínsky profil pacienta"),
            ui.output_ui("alerts"),
            
            ui.layout_columns(
                ui.card(ui.card_header("Základné informácie"), ui.output_ui("patient_info")),
                ui.card(ui.card_header("HFE Genetický status"), ui.output_ui("genetic_info"))
            ),
            ui.card(
                ui.card_header("História vyšetrení pacienta"),
                ui.output_data_frame("lab_table")
            ),
            ui.card(
                ui.card_header("Porovnanie pacienta s populáciou"),
                ui.input_select(
                    "patient_param",
                    "Laboratórny parameter (porovnanie s populáciou):",
                    choices=available_parameters,
                    selected="Železo_S" if "Železo_S" in available_parameters else available_parameters[0]
                ),
                ui.output_plot("iron_plot")
            )
        )
    ),

    #visualisation tab
    ui.nav_panel(
        "Skupinové Vizualizácie a Vzťahy",
        ui.layout_sidebar(
            ui.sidebar(
                ui.h4("Filtre grafov"),
                ui.input_select("parameter", "Laboratórny parameter:", choices=available_parameters, selected=available_parameters[0]),
                ui.input_select("plot_type", "Typ grafu:", choices=["Histogram", "Boxplot", "Boxplot podľa HFE statusu"], selected="Histogram"),
                ui.input_select("gender_filter", "Pohlavie:", choices=["Všetci", "M", "F", "mixed"], selected="Všetci"),
                ui.input_select("hfe_filter", "HFE status:", choices=["Všetci", "positive", "negative", "not_tested"], selected="Všetci"),
                ui.hr(),
                ui.input_select("x_parameter", "Scatter plot - os X:", choices=available_parameters, selected="Urea_S" if "Urea_S" in available_parameters else available_parameters[0]),
                ui.input_select("y_parameter", "Scatter plot - os Y:", choices=available_parameters, selected="Kreatinín_S" if "Kreatinín_S" in available_parameters else available_parameters[0]),
                width=320
            ),
            
            ui.card(
                ui.h5("Základné informácie o datasete"),
                ui.p(f"Celkový počet riadkov v datasete: {len(df):,}"),
                ui.p(f"Celkový počet pacientov: {profiles['Patient_ID'].nunique():,}"),
                ui.p(f"Počet rôznych laboratórnych vyšetrení: {lab_test_count:,}"),
                style="background-color: #f8f9fa; border-left: 5px solid #2f80ed;"
            ),
            
            ui.card(
                ui.card_header("1. Rozdelenie hodnôt vybraného parametra"),
                ui.output_plot("main_plot"),
                ui.output_text_verbatim("summary_text")
            ),
            
            ui.layout_columns(
                ui.card(
                    ui.card_header("2. Vzťah medzi dvoma parametrami"),
                    ui.output_plot("scatter_plot")
                ),
                ui.card(
                    ui.card_header("3. Korelačná heatmapa"),
                    ui.output_plot("correlation_heatmap")
                )
            )
        )
    ),
    
    title="SSBU Analýza Dát",
    id="main_nav",
    bg="#a3cbff"
)


def server(input, output, session):
    #patient logic
    @reactive.Calc
    def get_patient_profile():
        pid = input.patient_id()
        if not pid or profiles.empty: return None
        match = profiles[profiles['Patient_ID'] == pid]
        return match.iloc[0] if not match.empty else None

    @reactive.Calc
    def get_patient_labs():
        pid = input.patient_id()
        if not pid or df.empty: return pd.DataFrame()
        return df[df['Patient_ID'] == pid]

    @output
    @render.ui
    def patient_info():
        prof = get_patient_profile()
        if prof is None: return ui.p("Vyberte pacienta.")
        vek = 2026 - prof['Rok_Nar'] if pd.notna(prof['Rok_Nar']) else "Neznámy"
        return ui.HTML(f"""
            <ul style="font-size: 1.1em; line-height: 1.6;">
                <li><b>Pohlavie:</b> {prof.get('Pohlavie', 'Neznáme')}</li>
                <li><b>Rok narodenia:</b> {prof.get('Rok_Nar', 'Neznámy')} (Vek: ~{vek} r.)</li>
                <li><b>Počet návštev (odberov):</b> {prof.get('Pocet_vysetreni', 0)}</li>
                <li><b>Počet nameraných parametrov:</b> {prof.get('Pocet_riadkov', 0)}</li>
                <li><b>Počet laboratórnych parametrov:</b> {prof.get('Pocet_lab_parametrov', 0)}</li>
                <li><b>Prvý odber:</b> {prof.get('Prvy_odber', 'Chýba')}</li>
                <li><b>Posledný odber:</b> {prof.get('Posledny_odber', 'Chýba')}</li>
            </ul>
        """)

    @output
    @render.ui
    def genetic_info():
        prof = get_patient_profile()
        if prof is None: return ui.p("Vyberte pacienta.")
        status = str(prof.get('HFE_status', 'not_tested'))
        color = "#d9534f" if status in ['homozygote', 'positive', 'heterozygote'] else "#5cb85c"
        if status == 'not_tested': color = "#777"
        return ui.HTML(f"<h3 style='color:{color}; text-transform: uppercase;'>{status}</h3>")

    @output
    @render.ui
    def alerts():
        labs = get_patient_labs()
        prof = get_patient_profile()
        if labs.empty or prof is None: return ui.HTML("")
        alerts_html = ""
        fe_vals = labs[labs['Lab_Nazov'].str.contains('Železo', na=False, case=False)]['Vysl_numeric']
        crp_vals = labs[labs['Lab_Nazov'].str.contains('CRP', na=False, case=False)]['Vysl_numeric']
        max_fe = fe_vals.max() if not fe_vals.empty else 0
        max_crp = crp_vals.max() if not crp_vals.empty else 0
        is_hfe_pos = prof.get('HFE_status', '') in ['homozygote', 'positive']
        
        if max_fe > 27 and is_hfe_pos:
            alerts_html += f'<div style="background-color:#ffe6e6; padding:15px; border-left: 6px solid #d9534f; margin-bottom:15px;"><b>⚠️ Kritické upozornenie:</b> Potvrdená mutácia a vysoké železo ({max_fe} μmol/l). Riziko hemochromatózy.</div>'
        if max_fe > 25 and max_crp > 10:
            alerts_html += f'<div style="background-color:#fff8e6; padding:15px; border-left: 6px solid #f0ad4e; margin-bottom:15px;"><b>🔔 Upozornenie:</b> Vysoké železo ({max_fe} μmol/l) môže byť skreslené zápalom (CRP {max_crp} mg/l).</div>'
        return ui.HTML(alerts_html)

    @output
    @render.data_frame
    def lab_table():
        labs = get_patient_labs()
        if labs.empty: return pd.DataFrame()
        show_df = labs[['Datum_O_clean', 'Lab_Nazov', 'Vysl_numeric', 'Vysl_text_clean', 'Result_type']].copy()
        show_df.columns = ['Dátum odberu', 'Názov vyšetrenia', 'Numerický výsledok', 'Textový výsledok', 'Typ']
        return render.DataTable(show_df.sort_values(by='Dátum odberu', ascending=False), filters=True)

    @output
    @render.plot
    def iron_plot():
        labs = get_patient_labs()
        pid = input.patient_id()
        if df.empty: return None
        selected_param = input.patient_param()
        pop_data = df[df['Lab_Nazov'] == selected_param].dropna(subset=['Vysl_numeric']).copy()
        pop_data = pop_data[pop_data['HFE_status_patient'].isin(['homozygote', 'heterozygote', 'positive', 'negative'])]
        
        sns.set_theme(style="whitegrid")
        fig, ax = plt.subplots(figsize=(10, 5))
        sns.boxplot(data=pop_data, x='HFE_status_patient', y='Vysl_numeric', ax=ax, color='lightblue')
        ax.set_title(f"Hodnoty {selected_param}: Populácia vs. Vybraný pacient", pad=15)
        ax.set_ylabel(format_param_label(selected_param))
        if selected_param == "Železo_S":
            ax.axhline(27, color='#d9534f', linestyle='--', label='Kritická hranica (27)')
        
        patient_vals = labs[labs['Lab_Nazov'] == selected_param]['Vysl_numeric'].dropna()
        prof = get_patient_profile()
        p_status = prof.get('HFE_status', 'not_tested') if prof is not None else 'not_tested'
        if not patient_vals.empty and p_status in pop_data['HFE_status_patient'].unique():
            ax.scatter([p_status]*len(patient_vals), patient_vals, color='#333333', edgecolor='white', s=120, zorder=10, label=f'Pacient ({pid})')
        ax.legend()
        plt.tight_layout()
        return fig

    #visualisation logic
    @reactive.Calc
    def selected_data():
        return filter_data(input.parameter(), input.gender_filter(), input.hfe_filter())

    @output
    @render.plot
    def main_plot():
        data = selected_data()
        fig, ax = plt.subplots(figsize=(8, 5))
        if data.empty:
            ax.text(0.5, 0.5, "Pre zvolené filtre nie sú dostupné dáta.", ha="center", va="center")
            ax.set_axis_off()
            return fig

        plot_type = input.plot_type()
        if plot_type == "Histogram":
            ax.hist(data["Vysl_numeric"], bins=30, color="#2f80ed", edgecolor="white")
            ax.set_title(f"Histogram pre {input.parameter()}")
            ax.set_xlabel(format_param_label(input.parameter()))
            ax.set_ylabel("Počet meraní")
        elif plot_type == "Boxplot":
            # Fix pre Seaborn
            sns.boxplot(y=data["Vysl_numeric"], ax=ax, color='lightblue')
            ax.set_title(f"Boxplot pre {input.parameter()}")
            ax.set_ylabel(format_param_label(input.parameter()))
        elif plot_type == "Boxplot podľa HFE statusu":
            hfe_data = filter_data(input.parameter(), input.gender_filter(), "Všetci")
            hfe_data = hfe_data[hfe_data["HFE_status_patient"].isin(["positive", "negative"])]
            if hfe_data.empty:
                ax.text(0.5, 0.5, "Nie sú dostupné HFE údaje.", ha="center", va="center")
                ax.set_axis_off()
                return fig
            # Fix pre Seaborn
            sns.boxplot(data=hfe_data, x="HFE_status_patient", y="Vysl_numeric", ax=ax, color='lightblue')
            ax.set_title(f"{input.parameter()} podľa HFE statusu")
            ax.set_ylabel(format_param_label(input.parameter()))
        
        fig.tight_layout()
        return fig

    @output
    @render.text
    def summary_text():
        data = selected_data()
        if data.empty: return "Pre zvolené filtre nie sú dostupné dáta."
        return (f"Štatistika pre {input.parameter()}:\n"
                f"Počet záznamov: {len(data)} | Počet pacientov: {data['Patient_ID'].nunique()}\n"
                f"Priemer: {data['Vysl_numeric'].mean():.2f} | Medián: {data['Vysl_numeric'].median():.2f}\n"
                f"Min: {data['Vysl_numeric'].min():.2f} | Max: {data['Vysl_numeric'].max():.2f}")

    @output
    @render.plot
    def scatter_plot():
        x_param, y_param = input.x_parameter(), input.y_parameter()
        data = filter_long_data([x_param, y_param], input.gender_filter(), input.hfe_filter())
        pivot_data = data.pivot_table(index="Patient_ID", columns="Lab_Nazov", values="Vysl_numeric", aggfunc="mean")
        
        fig, ax = plt.subplots(figsize=(7, 5))
        if x_param not in pivot_data.columns or y_param not in pivot_data.columns:
            ax.text(0.5, 0.5, "Nedostatok dát pre túto dvojicu.", ha="center", va="center")
            ax.set_axis_off()
            return fig

        ax.scatter(pivot_data[x_param], pivot_data[y_param], alpha=0.6, color="#2f80ed")
        ax.set_title(f"Vzťah: {x_param} vs {y_param}")
        ax.set_xlabel(format_param_label(x_param))
        ax.set_ylabel(format_param_label(y_param))
        fig.tight_layout()
        return fig

    @output
    @render.plot
    def correlation_heatmap():
        data = filter_long_data(heatmap_parameters, input.gender_filter(), input.hfe_filter())
        pivot_data = data.pivot_table(index="Patient_ID", columns="Lab_Nazov", values="Vysl_numeric", aggfunc="mean")
        
        fig, ax = plt.subplots(figsize=(10, 7))
        if pivot_data.shape[1] < 2:
            ax.text(0.5, 0.5, "Nedostatok parametrov na koreláciu.", ha="center", va="center")
            ax.set_axis_off()
            return fig

        sns.heatmap(pivot_data.corr(), annot=True, cmap="coolwarm", fmt=".2f", ax=ax)
        ax.set_title("Korelačná matica")
        fig.tight_layout()
        return fig

app = App(app_ui, server)
