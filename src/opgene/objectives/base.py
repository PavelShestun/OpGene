from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Tuple, Optional

@dataclass
class EvaluationResult:
    score: float  # 0.0 to 1.0
    message: str
    is_critical: bool = False # Если True и score низкий -> штраф x100

class OptimizationObjective(ABC):
    """Базовый класс для любой спецификации оптимизации."""
    
    def __init__(self, weight: float = 1.0):
        self.weight = weight

    @abstractmethod
    def evaluate(self, sequence: str) -> EvaluationResult:
        pass
    
    @property
    def name(self) -> str:
        return self.__class__.__name__
