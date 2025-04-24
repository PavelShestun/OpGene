import os
from src.codon_optimizer import CodonOptimizationApp

if __name__ == "__main__":
    # Установите ваш email для Entrez
    os.environ["ENTREZ_EMAIL"] = "your_email@example.com"
    
    app = CodonOptimizationApp()
    app.run(show_plots=False)
