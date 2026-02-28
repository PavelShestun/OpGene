import math
from Bio.Data import CodonTable
from .base import OptimizationObjective, EvaluationResult

class CodonAdaptationObjective(OptimizationObjective):
    def __init__(self, codon_usage: dict, table_id: int, weight: float = 1.0, 
                 use_ramp: bool = True, mode: str = "maximize"):
        super().__init__(weight)
        self.mode = mode 
        self.codon_usage = codon_usage
        self.table_id = table_id
        self.use_ramp = use_ramp
        self.weights = self._calculate_relative_adaptiveness()
        self.ramp_length = 15 # Первые 15 кодонов

    def _calculate_relative_adaptiveness(self) -> dict:
        # Логика расчета W (relative adaptiveness)
        table = CodonTable.unambiguous_dna_by_id.get(self.table_id, CodonTable.standard_dna_table)
        weights = {}
        for aa, codons in self._get_synonymous_codons(table).items():
            freqs = [self.codon_usage.get(c, 0.01) for c in codons]
            max_freq = max(freqs) if freqs else 1.0
            for c in codons:
                weights[c] = self.codon_usage.get(c, 0.01) / max_freq
        return weights

    def _get_synonymous_codons(self, table) -> dict:
        syn = {}
        for codon, aa in table.forward_table.items():
            syn.setdefault(aa, []).append(codon)
        return syn

    def evaluate(self, sequence: str) -> EvaluationResult:
        if not sequence: return EvaluationResult(0.0, "Empty")

        # 1. Режим Максимизации (Классический CAI + Ramp)
        if self.mode == "maximize":
            ramp_nt = self.ramp_length * 3
            if self.use_ramp and len(sequence) > ramp_nt:
                cai_body = self._calculate_cai(sequence[ramp_nt:])
                cai_ramp = self._calculate_cai(sequence[:ramp_nt])
                penalty = 0.3 if cai_ramp > 0.8 else 0.0
                score = max(0, cai_body - penalty)
                return EvaluationResult(score, f"CAI: {cai_body:.2f} (Ramp OK)")
            else:
                cai = self._calculate_cai(sequence)
                return EvaluationResult(cai, f"CAI: {cai:.2f}")

        # 2. Режим Гармонизации (SOTA: Имитация ритма хоста)
        else:
            # Считаем частоты кодонов в текущей последовательности
            codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
            total = len(codons)
            our_usage = {c: codons.count(c)/total for c in set(codons)}
            
            # Считаем MSE (Mean Squared Error) между нами и профилем хоста
            error = 0.0
            for codon, host_freq in self.codon_usage.items():
                our_freq = our_usage.get(codon, 0.0)
                error += (host_freq - our_freq) ** 2
            
            # Чем меньше ошибка, тем выше балл (инвертируем MSE)
            # Ошибка 0.1 обычно уже считается большой для частот
            score = max(0.0, 1.0 - (error * 15)) 
            return EvaluationResult(score, f"Harmonized (MSE: {error:.4f})")

    def _eval_harmonize(self, sequence: str) -> EvaluationResult:
        # Считаем частоты кодонов в нашей текущей ДНК
        counts = {}
        total = 0
        for i in range(0, len(sequence), 3):
            c = sequence[i:i+3]
            counts[c] = counts.get(c, 0) + 1
            total += 1
        
        # Сравнение распределения (MSE между нашей ДНК и Хостом)
        error = 0
        for codon, host_freq in self.codon_usage.items():
            our_freq = counts.get(codon, 0) / total
            error += (host_freq - our_freq)**2
            
        score = max(0.0, 1.0 - error * 10)
        return EvaluationResult(score, f"Harmonized (MSE: {error:.4f})")

    def _calculate_cai(self, seq: str) -> float:
        if not seq: return 0.0
        log_w = 0.0
        count = 0
        for i in range(0, len(seq), 3):
            codon = seq[i:i+3]
            if len(codon) == 3 and codon in self.weights:
                log_w += math.log(max(self.weights[codon], 0.001))
                count += 1
        return math.exp(log_w / count) if count else 0.0
