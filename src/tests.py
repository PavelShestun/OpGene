import unittest
from codon_optimizer import *

class TestCodonOptimization(unittest.TestCase):
    def setUp(self):
        self.codon_usage = DEFAULT_CODON_USAGE
        self.aa_sequence = "MAGWSR"
        self.config = ConfigurationManager().load_config()
        self.specifications = [
            CodonUsageSpecification(self.codon_usage),
            GcContentSpecification(self.config["target_gc_range"]),
            AvoidPatternSpecification(self.config["avoid_motifs"]),
            RbsSpecification(),
            RnaFoldingSpecification(**self.config["rna_folding"])
        ]
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

    def test_app_run(self):
        app = CodonOptimizationApp()
        app.run()

if __name__ == "__main__":
    unittest.main()
