# src/opgene/objectives/codon_pair.py
import math
from .base import OptimizationObjective, EvaluationResult

class CodonPairObjective(OptimizationObjective):
    """
    SOTA: Codon Pair Bias (CPB).
    Отрицательный CPS (Codon Pair Score) замедляет трансляцию.
    """
    def __init__(self, cps_table: dict, weight: float = 1.0):
        super().__init__(weight)
        self.cps_table = cps_table

    def evaluate(self, sequence: str) -> EvaluationResult:
        if not self.cps_table:
            return EvaluationResult(1.0, "No CPS data available")

        # Извлекаем все пары кодонов (дикодоны)
        scores = []
        for i in range(0, len(sequence) - 6, 3):
            pair = sequence[i:i+6]
            # Получаем CPS из таблицы. Если пары нет — берем 0 (нейтрально)
            cps = self.cps_table.get(pair, 0.0)
            scores.append(cps)

        if not scores:
            return EvaluationResult(0.5, "Sequence too short for CPB")

        # Средний CPS по всей длине
        avg_cps = sum(scores) / len(scores)

        # Переводим логарифмический CPS в оценку 0..1 через сигмоиду
        # В индустрии CPS > 0 — хорошо, CPS < -0.1 — уже риск замедления
        score = 1.0 / (1.0 + math.exp(-avg_cps * 10))
        
        return EvaluationResult(score, f"Avg CPS: {avg_cps:.3f}")
