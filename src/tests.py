import unittest
import subprocess
import logging
from codon_optimizer import *

# Настройка логирования
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class TestCodonOptimization(unittest.TestCase):
    def setUp(self):
        self.codon_usage = DEFAULT_CODON_USAGE
        self.aa_sequence = "MAGWSR"
        self.config = ConfigurationManager().load_config()
        # Проверяем, доступен ли RNAfold, чтобы избежать ошибок
        try:
            subprocess.run(["RNAfold", "--version"], capture_output=True, check=True)
            rna_folding_available = True
        except (subprocess.SubprocessError, FileNotFoundError):
            rna_folding_available = False
            logger.warning("RNAfold не найден. Тесты для RnaFoldingSpecification будут пропущены.")

        self.specifications = [
            CodonUsageSpecification(self.codon_usage),
            GcContentSpecification(self.config["target_gc_range"]),
            AvoidPatternSpecification(self.config["avoid_motifs"]),
            RbsSpecification(),
        ]
        # Добавляем RnaFoldingSpecification только если RNAfold доступен
        if rna_folding_available:
            self.specifications.append(RnaFoldingSpecification(**self.config["rna_folding"]))
        else:
            # Удаляем вес для RnaFolding5Prime, если спецификация недоступна
            self.config["spec_weights"] = {
                k: v for k, v in self.config["spec_weights"].items() if k != "RnaFolding5Prime"
            }

        self.weights = self.config["spec_weights"]
        self.output = OutputManager()

    def test_codon_usage_spec(self):
        spec = CodonUsageSpecification(self.codon_usage)
        dna = "ATGGCAGGTTGGAGCCGT"  # MAGWSR
        eval_result = spec.evaluate(dna)
        self.assertGreaterEqual(eval_result.score, 0.0)
        self.assertLessEqual(eval_result.score, 1.0)
        print(f"CodonUsage оценка: {eval_result.score:.4f}")

    def test_gc_content_spec(self):
        spec = GcContentSpecification([0.45, 0.55])
        dna = "ATGGCAGGTTGGAGCCGT"  # GC = 11/18 ~ 0.611
        eval_result = spec.evaluate(dna)
        self.assertLess(eval_result.score, 1.0)  # Ожидаем штраф
        print(f"GC Content оценка: {eval_result.score:.4f}")

    def test_optimizer(self):
        optimizer = CodonOptimizer(
            aa_sequence=self.aa_sequence,
            codon_usage=self.codon_usage,
            specifications=self.specifications,
            weights=self.weights,
            num_processes=2
        )
        dna, metrics = optimizer.optimize_ma(**self.config["ma_parameters"])
        self.assertEqual(len(dna), len(self.aa_sequence) * 3)
        self.assertGreater(metrics["fitness"], 0.0)
        print(f"Оптимизированная ДНК: {dna}")
        print(f"Приспособленность: {metrics['fitness']:.4f}, GC: {metrics['gc_content']:.3f}, CAI: {metrics['cai']:.4f}")
        self.output.save_metrics_plot(metrics, "Тестовые метрики", "test_metrics.png")

    def test_app_run_ecoli(self):
        app = CodonOptimizationApp(organism="Escherichia coli K-12", organism_id="83333")
        # Уменьшаем количество последовательностей для тестов
        app.test_sequences = {"ShortPeptide": "MAGWSR"}
        app.run()

    def test_app_run_bsubtilis(self):
        app = CodonOptimizationApp(organism="Bacillus subtilis", organism_id="1423")
        # Уменьшаем количество последовательностей для тестов
        app.test_sequences = {"ShortPeptide": "MAGWSR"}
        app.run()

# Запуск тестов
if __name__ == "__main__":
    # Для Jupyter Notebook и консоли используем TestLoader
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCodonOptimization)
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)