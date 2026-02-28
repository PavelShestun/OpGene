import argparse
import logging
from src.opgene.factory import OrganismFactory
from src.opgene.objectives.codon_usage import CodonAdaptationObjective
from src.opgene.objectives.structure import GcContentObjective, RnaFoldingObjective
from src.opgene.objectives.motifs import MotifAvoidanceObjective
from src.opgene.algorithms.genetic import GeneticOptimizer

# Настройка логгера
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("OpGene")

def main():
    parser = argparse.ArgumentParser(description="OpGene: SOTA Codon Optimizer")
    parser.add_argument("--organism", type=str, required=True, help="Target organism name")
    parser.add_argument("--taxid", type=str, required=True, help="NCBI Taxonomy ID")
    parser.add_argument("--sequence", type=str, required=True, help="Amino acid sequence")
    args = parser.parse_args()

    # 1. Создание профиля организма
    logger.info(f"Loading profile for {args.organism}...")
    profile = OrganismFactory.create(args.organism, args.taxid)

    # 2. Настройка целей (Objectives)
    objectives = [
        CodonAdaptationObjective(
            profile.codon_usage, 
            profile.codon_table_id, 
            weight=1.5,
            use_ramp=profile.use_ramp_hypothesis
        ),
        GcContentObjective(
            target_range=profile.ideal_gc_range,
            weight=1.0
        ),
        MotifAvoidanceObjective(
            forbidden=profile.forbidden_motifs,
            suppress_cpg=profile.suppress_cpg,
            weight=5.0 # Критически важно
        ),
        RnaFoldingObjective(weight=2.0) # Требует ViennaRNA
    ]

    # 3. Инициализация оптимизатора
    optimizer = GeneticOptimizer(objectives, pop_size=50, generations=30)

    # 4. Запуск
    logger.info("Starting optimization...")
    dna_seq, metrics = optimizer.run(args.sequence, profile)

    # 5. Результат
    print("\n--- Result ---")
    print(f"Optimized DNA: {dna_seq}")
    print(f"Final Fitness: {metrics['fitness']:.4f}")

if __name__ == "__main__":
    main()
