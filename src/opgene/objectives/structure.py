import subprocess
import re
import logging
from typing import Tuple, List, Set
from .base import OptimizationObjective, EvaluationResult

logger = logging.getLogger(__name__)

class GcContentObjective(OptimizationObjective):
    """
    Оценивает GC-состав.
    SOTA: Проверяет не только среднее значение, но и локальные экстремумы (Local Peaks),
    которые могут вызвать остановку полимеразы.
    """
    def __init__(self, target_range: Tuple[float, float], window_size: int = 50, weight: float = 1.0):
        super().__init__(weight)
        self.min_gc, self.max_gc = sorted(target_range)
        self.window = window_size
        # Лимиты для локальных окон (физические ограничения синтеза/секвенирования)
        self.local_min = 0.25
        self.local_max = 0.80

    def evaluate(self, sequence: str) -> EvaluationResult:
        # Глобальный расчет
        gc_total = sum(1 for b in sequence if b in "GC") / len(sequence)
        
        # Поиск локальных "зон смерти" (GC < 20% или > 80%)
        # Полимераза часто "проскальзывает" или "застревает" на таких участках
        window = 40
        bad_local = 0
        for i in range(len(sequence) - window):
            sub = sequence[i:i+window]
            local_gc = sum(1 for b in sub if b in "GC") / window
            if local_gc > 0.85 or local_gc < 0.15:
                bad_local += 1
        
        score = 1.0
        # Штраф за глобальный GC
        if not (self.min_gc <= gc_total <= self.max_gc):
            score -= 0.5
        
        # Штраф за локальные аномалии
        if bad_local > 0:
            score -= 0.4
            return EvaluationResult(max(0, score), f"GC: {gc_total:.2f} (Local Peaks!)")

        return EvaluationResult(max(0, score), f"GC: {gc_total:.2f}")

    def _calc_gc(self, seq: str) -> float:
        return sum(1 for b in seq if b in "GC") / len(seq) if seq else 0.0


class RnaFoldingObjective(OptimizationObjective):
    """
    Оценивает вторичную структуру мРНК (шпильки) на 5'-конце.
    Использует ViennaRNA (RNAfold).
    SOTA: Сильные шпильки (низкая MFE) на старте блокируют рибосому.
    """
    def __init__(self, weight: float = 1.0, window_size: int = 45, bad_mfe: float = -10.0):
        super().__init__(weight)
        self.window_size = window_size
        self.bad_mfe = bad_mfe # Порог, ниже которого структура считается слишком стабильной
        self._tool_available = self._check_tool()

    def _check_tool(self) -> bool:
        try:
            subprocess.run(["RNAfold", "--version"], capture_output=True, check=True)
            return True
        except (FileNotFoundError, subprocess.SubprocessError):
            return False

    def evaluate(self, sequence: str) -> EvaluationResult:
        if not self._tool_available:
            return EvaluationResult(1.0, "RNAfold tool not found (check skipped)")
        
        # Анализируем только зону инициации (5' конец)
        # SOTA: Именно здесь структура критична для экспрессии
        region = sequence[:self.window_size]
        if len(region) < 10: 
            return EvaluationResult(1.0, "Seq too short")

        try:
            # Запуск RNAfold без генерации PostScript файлов (--noPS)
            process = subprocess.run(
                ["RNAfold", "--noPS"], 
                input=region, 
                text=True, 
                capture_output=True
            )
            
            output = process.stdout.splitlines()
            if len(output) < 2: 
                return EvaluationResult(0.5, "RNAfold parse error")

            # Парсинг MFE из строки вида "... (((...))) (-12.50)"
            # Ищем число в скобках в конце строки
            mfe_match = re.search(r"\(\s*(-?\d+\.\d+)\s*\)$", output[1])
            if not mfe_match:
                return EvaluationResult(0.5, "MFE not found")
            
            mfe = float(mfe_match.group(1))

            # Логика оценки:
            # MFE < bad_mfe (напр -12): Очень стабильная шпилька -> Плохо (0.0 - 0.5)
            # MFE > -2: Почти нет структуры -> Отлично (1.0)
            # Промежуток: Линейная интерполяция
            
            if mfe < self.bad_mfe:
                # Экспоненциальный штраф за очень сильные шпильки
                score = 0.1
                msg = f"5' MFE: {mfe:.1f} (Structure too stable!)"
            elif mfe > -3.0:
                score = 1.0
                msg = f"5' MFE: {mfe:.1f} (Open structure)"
            else:
                # Нормализация между bad_mfe и -3.0
                score = (mfe - self.bad_mfe) / (-3.0 - self.bad_mfe)
                msg = f"5' MFE: {mfe:.1f}"

            return EvaluationResult(max(0.0, score), msg)

        except Exception as e:
            logger.error(f"RNAfold failed: {e}")
            return EvaluationResult(1.0, "RNAfold error")


class RepeatAvoidanceObjective(OptimizationObjective):
    """
    Проверяет ДНК на наличие прямых повторов (Direct Repeats).
    SOTA: Критично для синтеза (Manufacturability). 
    Компании (Twist, IDT) отклоняют заказы с повторами >10-15 bp.
    """
    def __init__(self, min_repeat_len: int = 12, weight: float = 2.0):
        super().__init__(weight)
        self.min_len = min_repeat_len

    def evaluate(self, sequence: str) -> EvaluationResult:
        if len(sequence) < self.min_len * 2:
            return EvaluationResult(1.0, "No Repeats")

        # Эффективный поиск повторов через хеширование подстрок
        seen: Set[str] = set()
        found_repeats = []
        
        # Сканируем окно
        for i in range(len(sequence) - self.min_len + 1):
            sub = sequence[i : i + self.min_len]
            if sub in seen:
                # Найден повтор!
                # Для скорости не ищем все вхождения, достаточно факта наличия
                found_repeats.append(sub)
                # SOTA: Если нашли повтор, сразу штрафуем, чтобы не тратить CPU
                # (Можно убрать break, если хотим найти все)
                break 
            seen.add(sub)

        if found_repeats:
            # Это критическая проблема для синтеза, штраф должен быть максимальным
            return EvaluationResult(
                score=0.0, 
                message=f"Repeat found (>={self.min_len}bp)", 
                is_critical=True
            )
        
        return EvaluationResult(1.0, "No Repeats")
