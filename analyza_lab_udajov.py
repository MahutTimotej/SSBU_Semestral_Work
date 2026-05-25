import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import mannwhitneyu


class GeneticLaboratoryAnalysis:
    def __init__(self, output_dir="output"):
        self.output_dir = output_dir
        self.dataset_path = os.path.join(output_dir, "processed_final_dataset.csv")
        self.results_dir = os.path.join(output_dir, "genetic_lab_analysis")
        self.plots_dir = os.path.join(self.results_dir, "plots")

        self.selected_parameters = [
            "ALT_S",
            "AST_S",
            "GMT_S",
            "CRP_S",
            "Železo_S",
            "Feritin_S",
            "Transferín_S",
            "Saturácia Trf"
        ]

    def load_data(self):
        self.df = pd.read_csv(self.dataset_path, low_memory=False)
        print(f"Dáta načítané: {len(self.df)} riadkov")

    def prepare_data(self):
        # 1. genetické riadky, z nich zistíme HFE status pacienta
        genetic = self.df[self.df["HFE_status"].notna()][["Patient_ID", "HFE_status"]].copy()

        patient_hfe = (
            genetic
            .drop_duplicates()
            .groupby("Patient_ID")["HFE_status"]
            .agg(lambda x: "positive" if "positive" in list(x) else "negative")
            .reset_index()
        )

        # 2. laboratórne numerické riadky
        lab = self.df[self.df["Result_type"] == "numeric"].copy()
        lab["Vysl_numeric"] = pd.to_numeric(lab["Vysl_numeric"], errors="coerce")
        lab = lab.dropna(subset=["Patient_ID", "Lab_Nazov", "Vysl_numeric"])

        # 3. priradíme HFE status k laboratórnym dátam podľa pacienta
        lab = lab.merge(patient_hfe, on="Patient_ID", how="inner", suffixes=("", "_genetic"))

        if "HFE_status_genetic" in lab.columns:
            lab["HFE_status"] = lab["HFE_status_genetic"]

        # 4. necháme iba vybrané parametre
        lab = lab[lab["Lab_Nazov"].isin(self.selected_parameters)]

        self.lab = lab

        print(f"Laboratórne záznamy po priradení HFE: {len(self.lab)}")
        print(f"Počet pacientov: {self.lab['Patient_ID'].nunique()}")
        print("\nVybrané parametre:")
        print(self.lab["Lab_Nazov"].value_counts())

    def aggregate_to_patient_level(self):
        self.patient_data = (
            self.lab
            .groupby(["Patient_ID", "HFE_status", "Lab_Nazov"], as_index=False)
            .agg(
                median_value=("Vysl_numeric", "median"),
                count_values=("Vysl_numeric", "count")
            )
        )

        print(f"\nPacientské záznamy po agregácii: {len(self.patient_data)}")

    def mann_whitney_tests(self):
        results = []

        for parameter in self.selected_parameters:
            subset = self.patient_data[self.patient_data["Lab_Nazov"] == parameter]

            negative = subset[subset["HFE_status"] == "negative"]["median_value"]
            positive = subset[subset["HFE_status"] == "positive"]["median_value"]

            if len(negative) < 5 or len(positive) < 5:
                continue

            stat, p_value = mannwhitneyu(negative, positive, alternative="two-sided")

            neg_median = np.median(negative)
            pos_median = np.median(positive)

            if pos_median > neg_median:
                trend = "vyššie hodnoty pri HFE positive"
            elif pos_median < neg_median:
                trend = "nižšie hodnoty pri HFE positive"
            else:
                trend = "bez rozdielu v mediáne"

            results.append({
                "parameter": parameter,
                "n_negative": len(negative),
                "n_positive": len(positive),
                "negative_median": neg_median,
                "positive_median": pos_median,
                "negative_q1": np.percentile(negative, 25),
                "negative_q3": np.percentile(negative, 75),
                "positive_q1": np.percentile(positive, 25),
                "positive_q3": np.percentile(positive, 75),
                "mann_whitney_u": stat,
                "p_value": p_value,
                "significant": p_value < 0.05,
                "trend": trend
            })

        self.results = pd.DataFrame(results)

        print("\nVýsledky Mann-Whitney U testu:")
        print(self.results[[
            "parameter",
            "n_negative",
            "n_positive",
            "negative_median",
            "positive_median",
            "p_value",
            "significant",
            "trend"
        ]])

    def save_results(self):
        os.makedirs(self.results_dir, exist_ok=True)

        results_path = os.path.join(self.results_dir, "mann_whitney_results.csv")
        self.results.to_csv(results_path, index=False, encoding="utf-8-sig")

        summary_path = os.path.join(self.results_dir, "summary.txt")

        with open(summary_path, "w", encoding="utf-8") as file:
            file.write("Analýza genetických a laboratórnych údajov\n")
            file.write("=" * 50 + "\n\n")
            file.write("Porovnávali sa laboratórne hodnoty medzi HFE negative a HFE positive pacientmi.\n")
            file.write("Pre každého pacienta bol použitý medián hodnoty daného parametra.\n")
            file.write("Keďže dáta nemali normálne rozdelenie, bol použitý Mann-Whitney U test.\n\n")

            for _, row in self.results.iterrows():
                file.write(f"Parameter: {row['parameter']}\n")
                file.write(f"HFE negative medián: {row['negative_median']:.3f}\n")
                file.write(f"HFE positive medián: {row['positive_median']:.3f}\n")
                file.write(f"p-hodnota: {row['p_value']:.5f}\n")
                file.write(f"Štatisticky významné: {row['significant']}\n")
                file.write(f"Trend: {row['trend']}\n\n")

        print(f"\nVýsledky uložené do: {self.results_dir}")

    def create_boxplots(self):
        os.makedirs(self.plots_dir, exist_ok=True)

        for parameter in self.results["parameter"]:
            subset = self.patient_data[self.patient_data["Lab_Nazov"] == parameter]

            negative = subset[subset["HFE_status"] == "negative"]["median_value"]
            positive = subset[subset["HFE_status"] == "positive"]["median_value"]

            p_value = self.results[self.results["parameter"] == parameter]["p_value"].iloc[0]

            plt.figure(figsize=(7, 5))
            plt.boxplot(
                [negative, positive],
                tick_labels=["HFE negative", "HFE positive"],
                showfliers=True
            )

            plt.title(f"{parameter}\np = {p_value:.5f}")
            plt.ylabel("Pacientsky medián hodnoty")
            plt.xlabel("HFE status")
            plt.grid(axis="y", alpha=0.3)
            plt.tight_layout()

            safe_name = parameter.replace(" ", "").replace("/", "").replace("\\", "_")
            plt.savefig(os.path.join(self.plots_dir, f"boxplot_{safe_name}.png"), dpi=300)
            plt.close()

        print(f"Boxploty uložené do: {self.plots_dir}")

    def run(self):
        print("Spúšťam analýzu genetických a laboratórnych údajov...")
        print("-" * 60)

        self.load_data()
        self.prepare_data()
        self.aggregate_to_patient_level()
        self.mann_whitney_tests()
        self.save_results()
        self.create_boxplots()

        print("-" * 60)
        print("Analýza bola dokončená.")


if __name__ == "__main__":
    analysis = GeneticLaboratoryAnalysis(output_dir="output")
    analysis.run()