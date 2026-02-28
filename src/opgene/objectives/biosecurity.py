from .base import OptimizationObjective, EvaluationResult

class BiosecurityObjective(OptimizationObjective):
    """
    SOTA: Проверка на биобезопасность (Screening).
    Проверяет отсутствие сигнатур генов устойчивости к антибиотикам 
    и токсинов, которые могут вызвать блокировку заказа вендором.
    """
    def __init__(self, external_db: dict = None, weight: float = 10.0):
        super().__init__(weight)
        # Базовые угрозы + загруженные пользователем
        self.threat_signatures = {
            "AmpR_marker": "CGCGGAACCCCTATTTGTTTATTTTTCTAA",
            "KanR_marker": "GGTGAGATGGTGCTTTCGATTTAGTGC",
        }
        if external_db:
            self.threat_signatures.update(external_db)

    def evaluate(self, sequence: str) -> EvaluationResult:
        found_threats = []
        for name, sig in self.threat_signatures.items():
            if sig.upper() in sequence.upper():
                found_threats.append(name)
        
        if found_threats:
            return EvaluationResult(0.0, f"THREATS: {', '.join(found_threats)}", is_critical=True)
        return EvaluationResult(1.0, "Bio-Safety: Clear")
