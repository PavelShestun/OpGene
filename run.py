import os
from src.codon_optimizer import CodonOptimizationApp

if __name__ == "__main__":
    # Установите ваш email для Entrez
    if "ENTREZ_EMAIL" not in os.environ:
        os.environ["ENTREZ_EMAIL"] = "your_email@example.com"

    # Можно выбрать режим: "optimization" (максимизация CAI) или "harmonization" (соответствие распределению)
    mode = "optimization"
    
    app = CodonOptimizationApp()
    app.run(show_plots=False, mode=mode)
