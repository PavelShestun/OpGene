import random
import math
import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
from enum import Enum
from Bio.Data import CodonTable

# --- Настройка логирования ---
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger("OpGene")

# --- Данные (E. coli K-12 Codon Usage для теста) ---
DEFAULT_CODON_USAGE = {
    'TTT': 0.58, 'TTC': 0.42, 'TTA': 0.14, 'TTG': 0.13, 'CTT': 0.12, 'CTC': 0.10,
    'CTA': 0.04, 'CTG': 0.47, 'ATT': 0.49, 'ATC': 0.39, 'ATA': 0.11, 'ATG': 1.00,
    'GTT': 0.28, 'GTC': 0.20, 'GTA': 0.17, 'GTG': 0.35, 'TCT': 0.17, 'TCC': 0.15,
    'TCA': 0.14, 'TCG': 0.14, 'CCT': 0.18, 'CCC': 0.13, 'CCA': 0.20, 'CCG': 0.49,
    'ACT': 0.19, 'ACC': 0.40, 'ACA': 0.17, 'ACG': 0.25, 'GCT': 0.18, 'GCC': 0.26,
    'GCA': 0.23, 'GCG': 0.33, 'TAT': 0.59, 'TAC': 0.41, 'CAT': 0.57, 'CAC': 0.43,
    'CAA': 0.34, 'CAG': 0.66, 'AAT': 0.49, 'AAC': 0.51, 'AAA': 0.74, 'AAG': 0.26,
    'GAT': 0.63, 'GAC': 0.37, 'GAA': 0.68, 'GAG': 0.32, 'TGT': 0.46, 'TGC': 0.54,
    'TGG': 1.00, 'CGT': 0.36, 'CGC': 0.36, 'CGA': 0.07, 'CGG': 0.11, 'AGT': 0.16,
    'AGC': 0.25, 'AGA': 0.07, 'AGG': 0.04, 'GGT': 0.35, 'GGC': 0.37, 'GGA': 0.13,
    'GGG': 0.15
}

# --- 1. Models (Сущности) ---

class OrganismType(Enum):
    PROKARYOTE = "prokaryote"
    EUKARYOTE = "eukaryote"
    MAMMALIAN = "mammalian"

@dataclass
class OrganismProfile:
    name: str
    org_type: OrganismType
    codon_usage: Dict[str, float]
    codon_table_id: int = 11
    # SOTA параметры
    use_ramp_hypothesis: bool = True
    suppress_cpg: bool = False
    ideal_gc_range: tuple = (0.45, 0.55)
    forbidden_motifs: List[str] = field(default_factory=list)

# --- 2. Objectives (Стратегии оценки) ---

@dataclass
class EvaluationResult:
    score: float  # 0.0 to 1.0
    message: str

class OptimizationObjective(ABC):
    """Интерфейс стратегии"""
    def __init__(self, weight: float = 1.0):
        self.weight = weight

    @abstractmethod
    def evaluate(self, sequence: str) -> EvaluationResult:
        pass
    
    @property
    def name(self) -> str:
        return self.__class__.__name__

class CodonAdaptationObjective(OptimizationObjective):
    """CAI с учетом Ramp Hypothesis (медленный старт трансляции)"""
    def __init__(self, codon_usage: dict, table_id: int, weight: float = 1.0, use_ramp: bool = True):
        super().__init__(weight)
        self.use_ramp = use_ramp
        self.weights = self._calculate_weights(codon_usage, table_id)
        self.ramp_length = 15  # Первые 15 кодонов (~45 нуклеотидов)

    def _calculate_weights(self, usage, table_id):
        table = CodonTable.unambiguous_dna_by_id.get(table_id, CodonTable.standard_dna_table)
        weights = {}
        processed_aa = set()
        for codon, aa in table.forward_table.items():
            if aa in processed_aa: continue
            synonyms = [c for c, a in table.forward_table.items() if a == aa]
            freqs = [usage.get(c, 0.01) for c in synonyms]
            max_freq = max(freqs) if freqs else 0.01
            for c in synonyms:
                weights[c] = usage.get(c, 0.01) / max_freq
            processed_aa.add(aa)
        return weights

    def evaluate(self, sequence: str) -> EvaluationResult:
        if not sequence: return EvaluationResult(0.0, "Empty")
        
        # Разделяем последовательность на "Рампу" и "Тело"
        ramp_len_nt = self.ramp_length * 3
        if self.use_ramp and len(sequence) > ramp_len_nt:
            ramp_seq = sequence[:ramp_len_nt]
            body_seq = sequence[ramp_len_nt:]
            
            cai_body = self._calculate_cai(body_seq)
            cai_ramp = self._calculate_cai(ramp_seq)
            
            # SOTA: Рампа должна быть "мягкой" (CAI < 0.7), но не ужасной (CAI > 0.3)
            # Если рампа слишком быстрая (CAI > 0.8), рибосомы могут столкнуться
            penalty = 0.0
            if cai_ramp > 0.8: 
                penalty = 0.3 # Штраф за слишком быструю инициацию
                msg = f"Body CAI: {cai_body:.2f} | Ramp too fast ({cai_ramp:.2f})"
            else:
                msg = f"Body CAI: {cai_body:.2f} | Ramp OK ({cai_ramp:.2f})"
            
            final_score = max(0, cai_body - penalty)
            return EvaluationResult(final_score, msg)
        else:
            cai = self._calculate_cai(sequence)
            return EvaluationResult(cai, f"CAI: {cai:.2f}")

    def _calculate_cai(self, seq):
        log_w = 0
        count = 0
        for i in range(0, len(seq), 3):
            c = seq[i:i+3]
            if c in self.weights:
                log_w += math.log(max(self.weights[c], 0.001))
                count += 1
        return math.exp(log_w/count) if count else 0

class GcContentObjective(OptimizationObjective):
    """GC-состав: проверка общего процента и локальных пиков"""
    def __init__(self, target_range: Tuple[float, float], weight: float = 1.0):
        super().__init__(weight)
        self.min_gc, self.max_gc = target_range
        self.window_size = 40

    def evaluate(self, sequence: str) -> EvaluationResult:
        total_gc = sum(1 for b in sequence if b in "GC") / len(sequence)
        
        # Локальная проверка (SOTA)
        max_deviation = 0.0
        for i in range(len(sequence) - self.window_size):
            sub = sequence[i:i+self.window_size]
            local_gc = sum(1 for b in sub if b in "GC") / len(sub)
            if local_gc > 0.8 or local_gc < 0.2: # Экстремальные значения
                max_deviation = 1.0
        
        score = 1.0
        if not (self.min_gc <= total_gc <= self.max_gc):
            score -= abs(total_gc - (self.min_gc + self.max_gc)/2) * 2
        
        if max_deviation > 0:
            score -= 0.3 # Штраф за локальные пики
            return EvaluationResult(max(0, score), f"GC: {total_gc:.2f} (Local peaks found!)")
        
        return EvaluationResult(max(0, score), f"GC: {total_gc:.2f}")

class MotifAvoidanceObjective(OptimizationObjective):
    """Избегание запрещенных мотивов и гомополимеров"""
    def __init__(self, motifs: List[str], weight: float = 1.0):
        super().__init__(weight)
        self.motifs = motifs
        self.homopolymers = [b*5 for b in "ACGT"] # AAAAA, CCCCC...

    def evaluate(self, sequence: str) -> EvaluationResult:
        penalty = 0.0
        found = []
        for m in self.motifs + self.homopolymers:
            if m in sequence:
                penalty += 0.5
                found.append(m)
        
        score = max(0.0, 1.0 - penalty)
        msg = f"Clean" if score == 1.0 else f"Found: {', '.join(found[:3])}"
        return EvaluationResult(score, msg)

# --- 3. Optimizer (Алгоритм) ---

class GeneticOptimizer:
    def __init__(self, objectives: List[OptimizationObjective], profile: OrganismProfile):
        self.objectives = objectives
        self.profile = profile
        self.back_table = self._build_back_table()

    def _build_back_table(self):
        table = CodonTable.unambiguous_dna_by_id.get(self.profile.codon_table_id, CodonTable.standard_dna_table)
        bt = {}
        for codon, aa in table.forward_table.items():
            if aa != '*': bt.setdefault(aa, []).append(codon)
        return bt

    def calculate_fitness(self, dna: str) -> float:
        total = 0.0
        for obj in self.objectives:
            res = obj.evaluate(dna)
            total += res.score * obj.weight
        return total

    def run(self, aa_seq: str, pop_size=20, generations=15) -> Tuple[str, Dict]:
        # 1. Init Population
        population = []
        # Add Greedy
        population.append("".join([self.back_table[aa][0] for aa in aa_seq]))
        # Add Random
        for _ in range(pop_size - 1):
            dna = "".join([random.choice(self.back_table[aa]) for aa in aa_seq])
            population.append(dna)

        best_dna = population[0]
        best_fit = self.calculate_fitness(best_dna)

        # 2. Evolution Loop
        for gen in range(generations):
            scores = [(dna, self.calculate_fitness(dna)) for dna in population]
            scores.sort(key=lambda x: x[1], reverse=True)
            
            if scores[0][1] > best_fit:
                best_dna, best_fit = scores[0]
            
            # Simple Selection & Mutation (Prototyping)
            new_pop = [scores[0][0], scores[1][0]] # Elitism
            while len(new_pop) < pop_size:
                parent = random.choice(scores[:5])[0] # Top 5
                # Mutation
                child_list = list(parent)
                idx = random.randint(0, len(aa_seq)-1) * 3
                aa = aa_seq[idx//3]
                child_list[idx:idx+3] = list(random.choice(self.back_table[aa]))
                new_pop.append("".join(child_list))
            
            population = new_pop
            if gen % 5 == 0:
                logger.info(f"Gen {gen}: Best Fitness = {best_fit:.4f}")

        # 3. Final Report
        metrics = {}
        for obj in self.objectives:
            res = obj.evaluate(best_dna)
            metrics[obj.name] = {"score": res.score, "msg": res.message}
        
        return best_dna, {"fitness": best_fit, "metrics": metrics}

# --- 4. Main (Facade) ---

def main():
    logger.info("--- Starting OpGene OOP Test ---")
    
    # 1. Создаем профиль (как бы из Фабрики)
    # Пример: Bacillus subtilis (Prokaryote)
    profile = OrganismProfile(
        name="Bacillus_subtilis_Mock",
        org_type=OrganismType.PROKARYOTE,
        codon_usage=DEFAULT_CODON_USAGE, # В реальности загружаем из Entrez
        use_ramp_hypothesis=True,
        forbidden_motifs=["GGGGGG", "CCCCCC"]
    )
    
    # 2. Настраиваем цели (Strategy)
    objectives = [
        CodonAdaptationObjective(profile.codon_usage, profile.codon_table_id, weight=1.5, use_ramp=True),
        GcContentObjective(target_range=(0.4, 0.6), weight=1.0),
        MotifAvoidanceObjective(motifs=profile.forbidden_motifs, weight=2.0)
    ]
    
    # 3. Тестовая последовательность (GFP фрагмент)
    aa_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRY"
    
    logger.info(f"Organism: {profile.name}")
    logger.info(f"Sequence Length: {len(aa_sequence)} aa")
    
    # 4. Запуск
    optimizer = GeneticOptimizer(objectives, profile)
    best_dna, results = optimizer.run(aa_sequence, pop_size=50, generations=20)
    
    # 5. Вывод
    print("\n" + "="*30)
    print("OPTIMIZATION RESULTS")
    print("="*30)
    print(f"Final DNA: {best_dna}")
    print(f"Total Fitness: {results['fitness']:.4f}")
    print("\nDetailed Metrics:")
    for name, data in results['metrics'].items():
        print(f"  [{name}]: {data['score']:.2f} -> {data['msg']}")
    print("="*30)

if __name__ == "__main__":
    main()
