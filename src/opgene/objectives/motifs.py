# src/opgene/objectives/motifs.py
from typing import List
from .base import OptimizationObjective, EvaluationResult

class MotifAvoidanceObjective(OptimizationObjective):
    """
    SOTA: Избегание запрещенных мотивов.
    Включает: сайты рестрикции, гомополимеры, CpG (для млекопитающих) 
    и внутренние сайты посадки рибосом (для бактерий).
    """
    def __init__(self, forbidden: List[str], suppress_cpg: bool = False, 
                 avoid_internal_rbs: bool = True, weight: float = 1.0):
        super().__init__(weight)
        self.forbidden = [m.upper() for m in forbidden]
        self.suppress_cpg = suppress_cpg
        self.avoid_internal_rbs = avoid_internal_rbs
        
        # Стандартные гомополимеры (вызывают "проскальзывание" при синтезе и ПЦР)
        self.homopolymers = ["AAAAAA", "TTTTTT", "GGGGGG", "CCCCCC"]
        
        # SOTA: Мотивы Шайна-Дальгарно (RBS). 
        # Их появление внутри гена может привести к заторам рибосом в бактериях.
        self.rbs_motifs = ["AGGAGG", "GGAGGA", "GAGGAG", "AGGAG", "GGAGG"]

    def evaluate(self, sequence: str) -> EvaluationResult:
        penalty = 0.0
        found = []
        sequence = sequence.upper()

        # 1. Проверка пользовательских мотивов и гомополимеров
        check_list = self.forbidden + self.homopolymers
        for m in check_list:
            if m in sequence:
                count = sequence.count(m)
                penalty += 0.5 * count
                if m not in found:
                    found.append(m)

        # 2. SOTA: Поиск внутренних RBS (пропускаем самое начало гена)
        if self.avoid_internal_rbs:
            # Игнорируем первые 20 нуклеотидов (где старт трансляции)
            body = sequence[20:]
            for rbs in self.rbs_motifs:
                if rbs in body:
                    count = body.count(rbs)
                    penalty += 0.4 * count
                    if "Internal-RBS" not in found:
                        found.append(f"Internal-RBS({rbs})")

        # 3. CpG Подавление (Mammalian SOTA)
        if self.suppress_cpg:
            cpg_count = sequence.count("CG")
            # Эмпирический порог: если динуклеотидов CG слишком много
            if (cpg_count / len(sequence)) > 0.05:
                penalty += 0.2
                found.append("High-CpG")

        score = max(0.0, 1.0 - penalty)
        msg = "Clean" if not found else f"Hits: {', '.join(found[:3])}"
        
        # Критичность: если найдены мотивы, это может блокировать синтез
        return EvaluationResult(score, msg, is_critical=(penalty > 0))
