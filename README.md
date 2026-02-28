# 🧬 OpGene Elite: SOTA Codon Optimization Suite

![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Streamlit](https://img.shields.io/badge/UI-Streamlit-FF4B4B)
![BioPython](https://img.shields.io/badge/Powered%20by-BioPython-orange)

**OpGene** is an industrial-grade, State-of-the-Art (SOTA) codon optimization and gene design suite. It utilizes a Memetic Genetic Algorithm to optimize amino acid sequences into DNA that is highly expressed in target organisms while ensuring the sequence is actually manufacturable by DNA synthesis vendors (like Twist Bioscience or IDT) and safe to use.

Whether you are expressing simple peptides in *E. coli* or complex human proteins in *CHO cells*, OpGene handles the biological complexity for you.

---

## ✨ Key Features

### 🔬 Biological Accuracy (Expression Optimization)
* **Dual CAI Modes:** Choose between **"Maximize"** (greedy optimization for simple proteins) and **"Harmonize"** (matches the host's natural codon frequency to prevent misfolding of complex proteins).
* **Ramp Hypothesis Integration:** Automatically uses slower, rarer codons for the first ~15 amino acids to prevent ribosome traffic jams during translation initiation.
* **Codon Pair Bias (CPS):** Evaluates dinucleotide frequencies to maintain optimal ribosome speed.
* **5' mRNA Stability:** Uses `ViennaRNA` to prevent tight RNA hairpins at the start codon that would block ribosome binding.

### 🏭 Manufacturability (Synthesis Readiness)
* **Local GC Content Peaks:** Scans sliding windows to prevent localized GC-rich (>80%) or AT-rich (<20%) regions that stall polymerases.
* **Direct Repeat Avoidance:** Strictly penalizes sequence repeats (≥12bp) to ensure successful physical DNA synthesis.
* **Motif Avoidance:** Automatically removes homopolymers (e.g., `AAAAAA`), internal RBS sites (in bacteria), and high CpG islands (in mammals).

### 🛡️ Biosecurity Screening
* Built-in screening for hazardous sequences (e.g., Antibiotic Resistance Markers, Toxins).
* Supports uploading custom threat signatures via JSON.

### 💻 User Interfaces
* **Web UI (Streamlit):** Interactive dashboard with Plotly visualizations, PDF reports, and GenBank export.
* **CLI:** Command-line interface for CI/CD pipelines and headless servers.
* **Batch Processing:** Optimize hundreds of sequences at once via multi-FASTA upload.

---

## ⚙️ Installation

### Prerequisites
* Python 3.8+
* [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) (Required for RNA folding metrics)

### 1. System Dependencies (Linux/Ubuntu)
```bash
sudo apt-get update
sudo apt-get install -y viennarna
```

### 2. Python Environment Setup
We recommend using Conda to easily manage the ViennaRNA dependency and Python packages.

```bash
git clone https://github.com/PavelShestun/OpGene.git
cd OpGene

# Create and activate conda environment
conda env create -f environment.yml
conda activate codon_optimization

# Install additional UI dependencies (if using the Web App)
pip install streamlit plotly pandas fpdf
```

---

## 🚀 Usage

### Option A: Web Interface (Recommended)
OpGene comes with a rich, interactive web application.
```bash
streamlit run app.py
```
* **Single Mode:** Paste your amino acid sequence, tweak weights, and watch the evolutionary algorithm run. Export results as FASTA, GenBank (`.gb`), or a comprehensive PDF report.
* **Batch Mode:** Upload a `.fasta` file containing multiple protein sequences to optimize them all sequentially and download a summary CSV.

### Option B: Command Line Interface (CLI)
For quick headless execution or script integration:
```bash
export ENTREZ_EMAIL="your_email@example.com"

python src/run.py \
  --organism "Escherichia coli K-12" \
  --taxid 83333 \
  --sequence "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRY" \
  --out optimized.fasta
```

---

## 🧠 Under the Hood (Architecture)

OpGene uses an Object-Oriented architecture based on the **Strategy Pattern** for objective evaluation:

* `OrganismFactory`: Automatically fetches and caches genomic data, Codon Usage Tables, and Codon Pair Scores (CPS) directly from NCBI Entrez using the target organism's Taxonomy ID.
* `GeneticOptimizer`: The core evolutionary engine. It starts with a mixed population (greedy + random) and applies selection, synonymous mutations, and elitism over multiple generations.
* **Objectives Framework** (`src/opgene/objectives/`):
  * `CodonAdaptationObjective`
  * `CodonPairObjective`
  * `GcContentObjective`
  * `MotifAvoidanceObjective`
  * `RnaFoldingObjective`
  * `RepeatAvoidanceObjective`
  * `BiosecurityObjective`

---

## 📁 Project Structure

```text
OpGene/
├── app.py                      # Streamlit Web Dashboard
├── src/
│   ├── run.py                  # CLI Entrypoint
│   └── opgene/                 # Core Library
│       ├── algorithms/         # Genetic and Memetic optimizers
│       ├── objectives/         # Fitness functions (GC, CAI, CPB, etc.)
│       ├── utils/              # PDF reporting, GenBank export
│       ├── factory.py          # Organism profile creation
│       ├── data_loaders.py     # NCBI Entrez integrations
│       └── models.py           # Dataclasses & Enums
├── environment.yml             # Conda environment specs
└── export_codebase.py          # Utility script for LLM codebase dumping
```

---

## 🤝 Contributing
Contributions are welcome! If you want to add new fitness objectives (e.g., Deep Learning based translation efficiency predictors), please fork the repository and submit a Pull Request.

## 📄 License
This project is licensed under the MIT License. See the `LICENSE` file for details.
