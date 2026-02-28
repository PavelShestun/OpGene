import argparse
import logging
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Импорты модулей OpGene
from src.opgene.factory import OrganismFactory
from src.opgene.objectives.codon_usage import CodonAdaptationObjective
from src.opgene.objectives.structure import GcContentObjective, RnaFoldingObjective, RepeatAvoidanceObjective
from src.opgene.objectives.motifs import MotifAvoidanceObjective
from src.opgene.algorithms.genetic import GeneticOptimizer

# Настройка логирования
logging.basicConfig(
    level=logging.INFO, 
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("OpGene-CLI")

def main():
    parser = argparse.ArgumentParser(description="OpGene: SOTA Codon Optimizer (CLI)")
    parser.add_argument("--organism", default="Escherichia coli K-12", help="Target organism name")
    parser.add_argument("--taxid", default="83333", help="NCBI Taxonomy ID")
    parser.add_argument("--sequence", required=True, help="Amino acid sequence to optimize")
    parser.add_argument("--out", default="optimized.fasta", help="Output FASTA file")
    parser.add_argument("--email", default=os.getenv("ENTREZ_EMAIL", "test@example.com"), help="Email for NCBI Entrez")
    
    args = parser.parse_args()

    # 1. Создание профиля организма через Фабрику
    logger.info(f"Loading biological profile for {args.organism} (TaxID: {args.taxid})...")
    factory = OrganismFactory(args.email)
    profile = factory.create(args.organism, args.taxid)

    # 2. Настройка стратегий оптимизации (Objectives)
    # Мы используем веса по умолчанию, оптимизированные для промышленного синтеза
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
            weight=5.0  # Высокий вес для критических запретов
        ),
        RnaFoldingObjective(
            weight=1.5
        ),
        RepeatAvoidanceObjective(
            min_repeat_len=12, 
            weight=3.0  # Критично для синтеза ДНК
        )
    ]

    # 3. Инициализация оптимизатора
    # По умолчанию используем 50 особей и 30 поколений для баланса скорость/качество
    optimizer = GeneticOptimizer(objectives, profile, pop_size=50, generations=30)

    # 4. Запуск оптимизации
    aa_clean = "".join(args.sequence.split()).upper()
    logger.info(f"Starting evolution for sequence length {len(aa_clean)} aa...")
    
    dna, results = optimizer.run(aa_seq=aa_clean)

    # 5. Вывод результатов и сохранение
    print("\n" + "="*50)
    print("OPTIMIZATION COMPLETE")
    print("="*50)
    print(f"Final DNA: {dna}")
    print(f"Total Fitness: {results['fitness']:.4f}")
    print("\nDetailed Metrics:")
    for name, msg in results['metrics'].items():
        print(f"  - {name}: {msg}")

    # Сохранение в FASTA
    record = SeqRecord(
        Seq(dna), 
        id="OpGene_Optimized", 
        description=f"Optimized for {args.organism} CAI:{results['metrics'].get('CodonAdaptationObjective')}"
    )
    SeqIO.write(record, args.out, "fasta")
    logger.info(f"Saved optimized sequence to {args.out}")

if __name__ == "__main__":
    main()
