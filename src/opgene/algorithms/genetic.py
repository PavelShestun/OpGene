# src/opgene/algorithms/genetic.py
import random
from typing import List, Tuple, Dict
from Bio.Data import CodonTable
from ..objectives.base import OptimizationObjective
from ..models import OrganismProfile
from .base import BaseOptimizer

class GeneticOptimizer(BaseOptimizer):
    def __init__(self, objectives: List[OptimizationObjective], profile: OrganismProfile, pop_size=30, generations=20):
        super().__init__(objectives)
        self.profile = profile
        self.pop_size = pop_size
        self.generations = generations
        self.back_table = self._build_bt()

    def _build_bt(self):
        t = CodonTable.unambiguous_dna_by_id.get(self.profile.codon_table_id, CodonTable.standard_dna_table)
        bt = {}
        for c, aa in t.forward_table.items():
            if aa != '*': bt.setdefault(aa, []).append(c)
        return bt

    def calculate_fitness(self, dna: str) -> float:
        total = 0.0
        for obj in self.objectives:
            res = obj.evaluate(dna)
            # SOTA: Если критическая метрика провалена, обнуляем фитнес
            if getattr(res, 'is_critical', False) and res.score < 0.5:
                return 0.0
            total += res.score * obj.weight
        return total

    def run(self, aa_seq: str) -> Tuple[str, Dict]: # Убрал profile из аргументов run, он есть в init
        # 1. Init
        population = []
        
        # Особь 1: Жадная (Самые частые кодоны)
        population.append("".join([self.back_table[aa][0] for aa in aa_seq]))
        
        # Остальные: Случайные
        for _ in range(self.pop_size - 1):
            dna = "".join([random.choice(self.back_table[aa]) for aa in aa_seq])
            population.append(dna)

        best_dna = population[0]
        best_fit = self.calculate_fitness(best_dna)

        # 2. Loop
        for gen in range(self.generations):
            # Evaluate
            scores = [(d, self.calculate_fitness(d)) for d in population]
            scores.sort(key=lambda x: x[1], reverse=True)
            
            current_best, current_fit = scores[0]
            if current_fit > best_fit:
                best_dna, best_fit = current_best, current_fit

            # Next Gen (Elitism)
            new_pop = [scores[0][0], scores[1][0]]
            
            # Tournament / Breeding
            while len(new_pop) < self.pop_size:
                # Tournament selection
                candidates = random.sample(scores, 3)
                parent = max(candidates, key=lambda x: x[1])[0]
                
                # Mutation (SOTA: Synonymous swap)
                child_list = list(parent)
                # Мутируем 1-3 кодона
                for _ in range(random.randint(1, 3)):
                    idx = random.randint(0, len(aa_seq)-1) * 3
                    aa = aa_seq[idx//3]
                    opts = self.back_table[aa]
                    if len(opts) > 1:
                        child_list[idx:idx+3] = list(random.choice(opts))
                
                new_pop.append("".join(child_list))
            
            population = new_pop

        # 3. Report
        metrics = {obj.name: obj.evaluate(best_dna).message for obj in self.objectives}
        return best_dna, {"fitness": best_fit, "metrics": metrics}
