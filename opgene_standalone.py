import json
import os
import random
import gzip
import shutil
import urllib.request
import matplotlib.pyplot as plt
import math
import time
import re
import subprocess
import statistics
import logging
from typing import List, Tuple, Dict, Union, Optional
from math import log, exp
from Bio import SeqIO, Entrez
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool, cpu_count
from concurrent.futures import ProcessPoolExecutor
from collections import OrderedDict
import functools
import psutil

# --- Настройка логирования ---
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# --- Константы ---
FALLBACK_ORGANISM = "Escherichia coli K-12"
FALLBACK_ORGANISM_ID = "83333"
FALLBACK_FASTA_FILENAME = "E_coli_K12_MG1655.fna"
CPS_TABLE_FILENAME = "cps_table.json"
OUTPUT_DIR = "optimization_results"
DATA_DIR = "data"
CONFIG_FILE = "config.json"

# Таблица использования кодонов по умолчанию (E. coli K-12 MG1655)
DEFAULT_CODON_USAGE = {
    'TTT': 0.026, 'TTC': 0.016, 'TTA': 0.013, 'TTG': 0.013, 'CTT': 0.011, 'CTC': 0.010,
    'CTA': 0.004, 'CTG': 0.050, 'ATT': 0.030, 'ATC': 0.024, 'ATA': 0.007, 'ATG': 0.027,
    'GTT': 0.018, 'GTC': 0.015, 'GTA': 0.011, 'GTG': 0.023, 'TCT': 0.015, 'TCC': 0.014,
    'TCA': 0.012, 'TCG': 0.013, 'CCT': 0.017, 'CCC': 0.013, 'CCA': 0.020, 'CCG': 0.023,
    'ACT': 0.019, 'ACC': 0.023, 'ACA': 0.014, 'ACG': 0.015, 'GCT': 0.016, 'GCC': 0.026,
    'GCA': 0.021, 'GCG': 0.031, 'TAT': 0.018, 'TAC': 0.012, 'CAT': 0.013, 'CAC': 0.010,
    'CAA': 0.014, 'CAG': 0.029, 'AAT': 0.016, 'AAC': 0.021, 'AAA': 0.032, 'AAG': 0.011,
    'GAT': 0.035, 'GAC': 0.019, 'GAA': 0.039, 'GAG': 0.019, 'TGT': 0.005, 'TGC': 0.006,
    'TGG': 0.013, 'CGT': 0.020, 'CGC': 0.021, 'CGA': 0.006, 'CGG': 0.010, 'AGT': 0.009,
    'AGC': 0.015, 'AGA': 0.004, 'AGG': 0.002, 'GGT': 0.018, 'GGC': 0.027, 'GGA': 0.011,
    'GGG': 0.012
}

# --- Менеджер конфигурации ---
class ConfigurationManager:
    def __init__(self, config_file: str = CONFIG_FILE):
        self.config_file = config_file
        self.default_config = {
            "spec_weights": {
                "CodonUsage": 1.2,
                "GcContent": 1.0,
                "AvoidPattern": 1.5,
                "RnaFolding5Prime": 0.8,
                "CodonPairBias": 1.0,
                "RbsSpecification": 0.5
            },
            "ma_parameters": {
                "population_size": 10,
                "num_generations": 5,
                "crossover_rate": 0.85,
                "mutation_rate": 0.06,
                "tournament_size": 4,
                "elitism_count": 2,
                "local_search_steps": 3
            },
            "rna_folding": {
                "window_size": 40,
                "bad_mfe": -12.0,
                "ideal_mfe": -4.0
            },
            "cpb_parameters": {
                "sigmoid_k": 3.0,
                "sigmoid_center": 0.0
            },
            "target_gc_range": [0.45, 0.55],
            "avoid_motifs": [
                "AATAAA", "GATC", "TATAAT", "GGGGGG", "CCCCCC", "AAAAAA", "TTTTTT"
            ],
            "default_codon_table": 11
        }

    def load_config(self) -> Dict:
        if not os.path.exists(self.config_file):
            logger.info(f"Файл конфигурации не найден по пути {self.config_file}. Используются значения по умолчанию.")
            return self.default_config
        try:
            with open(self.config_file, 'r') as f:
                config = json.load(f)
                def merge_dicts(default, custom):
                    for key, value in custom.items():
                        if isinstance(value, dict) and key in default:
                            # Переименование устаревших ключей в rna_folding
                            if key == "rna_folding":
                                if "bad_mfe_threshold" in value:
                                    value["bad_mfe"] = value.pop("bad_mfe_threshold")
                                    logger.warning("Ключ 'bad_mfe_threshold' устарел. Используйте 'bad_mfe'.")
                                if "ideal_mfe_threshold" in value:
                                    value["ideal_mfe"] = value.pop("ideal_mfe_threshold")
                                    logger.warning("Ключ 'ideal_mfe_threshold' устарел. Используйте 'ideal_mfe'.")
                            default[key] = merge_dicts(default[key], value)
                        else:
                            default[key] = value
                    return default
                merged_config = merge_dicts(self.default_config.copy(), config)
                if not isinstance(merged_config["target_gc_range"], list) or len(merged_config["target_gc_range"]) != 2:
                    logger.warning(f"Неверный формат target_gc_range: {merged_config['target_gc_range']}. Используется значение по умолчанию [0.45, 0.55].")
                    merged_config["target_gc_range"] = [0.45, 0.55]
                return merged_config
        except Exception as e:
            logger.warning(f"Не удалось загрузить файл конфигурации {self.config_file}: {e}. Используются значения по умолчанию.")
            return self.default_config

# --- Менеджер ресурсов ---
class ResourceManager:
    def __init__(self, organism: str = FALLBACK_ORGANISM, organism_id: str = FALLBACK_ORGANISM_ID,
                 data_dir: str = DATA_DIR, codon_table_id: int = 11):
        self.organism = organism
        self.organism_id = organism_id
        self.data_dir = data_dir
        self.codon_table_id = codon_table_id
        self.fasta_path = os.path.join(data_dir, f"{organism.replace(' ', '_')}.fna")
        self.cps_path = os.path.join(data_dir, CPS_TABLE_FILENAME)
        os.makedirs(data_dir, exist_ok=True)
        self.entrez_email = self._get_entrez_email()
        self._download_fasta_if_needed()

    def _get_entrez_email(self) -> str:
        email = os.getenv("ENTREZ_EMAIL", None)
        if not email:
            logger.warning("Email для Entrez не указан. Запросы к NCBI могут быть ограничены.")
            return "example@example.com"
        return email

    def _download_fasta_if_needed(self):
        if os.path.exists(self.fasta_path):
            logger.info(f"FASTA файл для {self.organism} уже существует по пути {self.fasta_path}.")
            return
        logger.info(f"Загрузка FASTA файла для {self.organism} (TaxID: {self.organism_id})...")
        try:
            Entrez.email = self.entrez_email
            search_handle = Entrez.esearch(db="nucleotide", term=f"txid{self.organism_id}[Organism] AND chromosome[Title]", retmax=1)
            search_results = Entrez.read(search_handle)
            search_handle.close()
            if not search_results["IdList"]:
                logger.warning(f"Не удалось найти геном для {self.organism} (TaxID: {self.organism_id}).")
                if self.organism_id == FALLBACK_ORGANISM_ID:
                    url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
                    self._download_fallback_fasta(url)
                return
            genome_id = search_results["IdList"][0]
            fetch_handle = Entrez.efetch(db="nucleotide", id=genome_id, rettype="fasta", retmode="text")
            raw_data = fetch_handle.read()
            fetch_handle.close()

            # Предобработка: удаляем комментарии перед первой последовательностью
            lines = raw_data.splitlines()
            cleaned_lines = []
            in_sequence = False
            for line in lines:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    in_sequence = True
                if in_sequence:
                    cleaned_lines.append(line)
                elif not (line.startswith("#") or line.startswith(";") or line.startswith("!")):
                    logger.warning(f"Неожиданная строка в начале FASTA файла: {line}")
            cleaned_data = "\n".join(cleaned_lines)

            with open(self.fasta_path, "w") as f:
                f.write(cleaned_data)
            logger.info(f"FASTA файл успешно загружен по пути {self.fasta_path}")
        except Exception as e:
            logger.warning(f"Ошибка загрузки FASTA файла через Entrez: {e}")
            if self.organism_id == FALLBACK_ORGANISM_ID:
                url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
                self._download_fallback_fasta(url)

    def _download_fallback_fasta(self, url: str):
        logger.info(f"Загрузка резервного FASTA файла для {self.organism}...")
        gz_path = self.fasta_path + ".gz"
        try:
            urllib.request.urlretrieve(url, gz_path)
            with gzip.open(gz_path, 'rb') as f_in, open(self.fasta_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(gz_path)
            logger.info(f"Резервный FASTA файл успешно загружен по пути {self.fasta_path}")
        except Exception as e:
            logger.warning(f"Не удалось загрузить резервный FASTA файл: {e}")

    def _process_sequence(self, seq: str) -> Tuple[Dict[str, int], int]:
        """Обработка одной последовательности для подсчёта кодонов (для параллельного выполнения)."""
        codon_counts = {}
        total_codons = 0
        if len(seq) % 3 != 0:
            return codon_counts, total_codons
        for i in range(0, len(seq) - 3, 3):
            codon = seq[i:i+3].upper()
            if 'N' not in codon:
                codon_counts[codon] = codon_counts.get(codon, 0) + 1
                total_codons += 1
        return codon_counts, total_codons

    def _merge_codon_counts(self, results: List[Tuple[Dict[str, int], int]]) -> Tuple[Dict[str, int], int]:
        """Объединение результатов подсчёта кодонов из параллельных процессов."""
        merged_codon_counts = {}
        total_codons = 0
        for codon_counts, count in results:
            for codon, freq in codon_counts.items():
                merged_codon_counts[codon] = merged_codon_counts.get(codon, 0) + freq
            total_codons += count
        return merged_codon_counts, total_codons

    def _load_codon_usage_from_file(self) -> Optional[Dict[str, float]]:
        """Попытка загрузить таблицу кодонов из локального файла."""
        local_file = os.path.join(self.data_dir, f"{self.organism.replace(' ', '_')}_codon_usage.txt")
        if not os.path.exists(local_file):
            return None
        try:
            codon_usage = {}
            with open(local_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split()
                    if len(parts) >= 2:
                        codon, freq = parts[0], float(parts[1])
                        codon_usage[codon.upper()] = freq
            total = sum(codon_usage.values())
            if total > 0:
                return {codon: freq / total for codon, freq in codon_usage.items()}
            return None
        except Exception as e:
            logger.warning(f"Ошибка загрузки локальной таблицы кодонов из {local_file}: {e}")
            return None

    def load_codon_usage(self, force_recalculate: bool = False) -> Dict[str, float]:
        codon_usage_path = os.path.join(self.data_dir, f"{self.organism.replace(' ', '_')}_codon_usage.json")
        if os.path.exists(codon_usage_path) and not force_recalculate:
            logger.info(f"Загрузка сохранённой таблицы кодонов для {self.organism} из {codon_usage_path}")
            with open(codon_usage_path, 'r') as f:
                return json.load(f)

        # Попытка загрузить из локального файла
        local_codon_usage = self._load_codon_usage_from_file()
        if local_codon_usage:
            logger.info(f"Используется локальная таблица кодонов для {self.organism}")
            with open(codon_usage_path, 'w') as f:
                json.dump(local_codon_usage, f, indent=2)
            return local_codon_usage

        logger.info(f"Расчёт таблицы использования кодонов для {self.organism} (TaxID: {self.organism_id})...")
        try:
            Entrez.email = self.entrez_email
            # Поиск генов для организма (уменьшили retmax для ускорения)
            search_handle = Entrez.esearch(db="nucleotide", term=f"txid{self.organism_id}[Organism] AND gene[Feature]", retmax=100)
            search_results = Entrez.read(search_handle)
            search_handle.close()
            if not search_results["IdList"]:
                logger.warning(f"Не удалось найти гены для {self.organism}. Используется таблица по умолчанию.")
                return DEFAULT_CODON_USAGE

            # Пакетная загрузка генов (по 50 ID за раз)
            gene_ids = search_results["IdList"]
            batch_size = 50
            sequences = []
            for i in range(0, len(gene_ids), batch_size):
                batch_ids = gene_ids[i:i + batch_size]
                try:
                    fetch_handle = Entrez.efetch(
                        db="nucleotide",
                        id=",".join(batch_ids),  # Пакетная загрузка
                        rettype="fasta_cds_na",
                        retmode="text"
                    )
                    for record in SeqIO.parse(fetch_handle, "fasta"):
                        seq = str(record.seq)
                        if len(seq) < 30 or len(seq) % 3 != 0:  # Пропускаем короткие или некорректные последовательности
                            continue
                        sequences.append(seq)
                    fetch_handle.close()
                except Exception as e:
                    logger.warning(f"Ошибка загрузки пакета генов {batch_ids}: {e}")

            if not sequences:
                logger.warning(f"Не удалось собрать данные о кодонах для {self.organism}. Используется таблица по умолчанию.")
                return DEFAULT_CODON_USAGE

            # Параллельная обработка последовательностей
            num_processes = min(cpu_count(), 4)  # Ограничиваем количество процессов
            with Pool(processes=num_processes) as pool:
                results = pool.map(self._process_sequence, sequences)

            # Объединение результатов
            codon_counts, total_codons = self._merge_codon_counts(results)

            if total_codons < 1000:  # Минимальный порог для надёжной статистики
                logger.warning(f"Собрано слишком мало кодонов ({total_codons}) для {self.organism}. Используется таблица по умолчанию.")
                return DEFAULT_CODON_USAGE

            # Нормализация частот
            codon_usage = {codon: count / total_codons for codon, count in codon_counts.items()}
            with open(codon_usage_path, 'w') as f:
                json.dump(codon_usage, f, indent=2)
            logger.info(f"Таблица использования кодонов сохранена по пути {codon_usage_path}")
            return codon_usage
        except Exception as e:
            logger.warning(f"Ошибка загрузки данных через Entrez: {e}. Используется таблица по умолчанию.")
            return DEFAULT_CODON_USAGE

    def calculate_codon_pair_scores(self) -> Optional[Dict[str, float]]:
        if not os.path.exists(self.fasta_path):
            logger.info("FASTA файл недоступен. Расчёт CPB пропущен.")
            return None
        try:
            cps_table = {}
            total_pairs = 0
            # Try fasta-blast first, fallback to fasta
            try:
                records = list(SeqIO.parse(self.fasta_path, "fasta-blast"))
            except:
                records = list(SeqIO.parse(self.fasta_path, "fasta"))

            for record in records:
                seq = str(record.seq)
                for i in range(0, len(seq) - 6, 3):
                    pair = seq[i:i+6].upper()
                    if 'N' not in pair and len(pair) == 6:
                        cps_table[pair] = cps_table.get(pair, 0) + 1
                        total_pairs += 1
            if total_pairs == 0:
                return None
            for pair in cps_table:
                cps_table[pair] = cps_table[pair] / total_pairs
            with open(self.cps_path, 'w') as f:
                json.dump(cps_table, f, indent=2)
            logger.info(f"Таблица CPB сохранена по пути {self.cps_path}")
            return cps_table
        except Exception as e:
            logger.warning(f"Ошибка расчёта CPB: {e}. CPB отключен.")
            return None

# --- Классы спецификаций ---
class Pattern:
    def __init__(self, motif: str):
        self.motif = motif.upper()
        self.revcomp = self._reverse_complement(motif)
        self._regex = re.compile(f"(?={self.motif})") if motif else None
        self._revcomp_regex = re.compile(f"(?={self.revcomp})") if self.revcomp != self.motif else None

    def _reverse_complement(self, seq: str) -> str:
        complement = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
        return seq.translate(complement)[::-1]

    def find_matches(self, sequence: str) -> List[Tuple[int, int]]:
        sequence = sequence.upper()
        matches = set()
        if self._regex:
            for m in self._regex.finditer(sequence):
                matches.add((m.start(), m.start() + len(self.motif)))
        if self._revcomp_regex:
            for m in self._revcomp_regex.finditer(sequence):
                matches.add((m.start(), m.start() + len(self.revcomp)))
        return sorted(list(matches))

class SpecEvaluation:
    def __init__(self, score: float, message: str = "", locations: Optional[List[Tuple[int, int]]] = None):
        self.score = max(0.0, min(1.0, float(score) if not math.isnan(score) else 0.0))
        self.message = str(message)
        self.locations = locations or []

    def as_dict(self) -> Dict:
        return {"score": self.score, "message": self.message, "locations": self.locations}

class Specification:
    def evaluate(self, sequence: str) -> SpecEvaluation:
        raise NotImplementedError

    def get_name(self) -> str:
        return self.__class__.__name__.replace("Specification", "")

class CodonUsageSpecification(Specification):
    def __init__(self, codon_usage: Dict[str, float], codon_table_id: int = 11, mode: str = "optimization"):
        self.codon_usage = codon_usage or DEFAULT_CODON_USAGE
        self.codon_table_id = codon_table_id
        self.mode = mode # "optimization" or "harmonization"
        self.codon_weights = self._compute_codon_weights()

    def _compute_codon_weights(self) -> Dict[str, float]:
        codon_weights = {}
        table = CodonTable.unambiguous_dna_by_id.get(self.codon_table_id, CodonTable.standard_dna_table)
        processed_aas = set()
        for codon, aa in table.forward_table.items():
            if aa in processed_aas:
                continue
            synonymous_codons = [c for c, a in table.forward_table.items() if a == aa]
            freq_list = [self.codon_usage.get(c, 1e-9) for c in synonymous_codons]
            max_freq = max(freq_list) if freq_list else 1e-9
            for c in synonymous_codons:
                codon_weights[c] = max(1e-9, self.codon_usage.get(c, 1e-9) / max_freq) if max_freq > 1e-10 else 1e-9
            processed_aas.add(aa)
        for codon in table.forward_table:
            if len([c for c, a in table.forward_table.items() if a == table.forward_table[codon]]) == 1:
                codon_weights[codon] = 1.0
        return codon_weights

    def evaluate(self, sequence: str) -> SpecEvaluation:
        if not sequence:
            return SpecEvaluation(0.0, f"{self.get_name()}: N/A (пустая последовательность)")
        if len(sequence) % 3 != 0:
            return SpecEvaluation(0.0, f"{self.get_name()}: Длина последовательности не кратна 3")

        if self.mode == "harmonization":
            score = self._calculate_harmonization(sequence)
            return SpecEvaluation(score, f"Harmonization: {score:.4f}")
        else:
            score = self._calculate_cai(sequence)
            return SpecEvaluation(score, f"CAI: {score:.4f}")

    def _calculate_harmonization(self, dna_seq: str) -> float:
        """Calculates how well the sequence matches the host codon usage distribution."""
        table = CodonTable.unambiguous_dna_by_id.get(self.codon_table_id, CodonTable.standard_dna_table)

        # Count codons in sequence
        counts = {}
        total_by_aa = {}
        for i in range(0, len(dna_seq), 3):
            codon = dna_seq[i:i+3].upper()
            if codon in table.stop_codons or codon not in self.codon_weights:
                continue
            aa = table.forward_table.get(codon)
            if not aa: continue
            counts[codon] = counts.get(codon, 0) + 1
            total_by_aa[aa] = total_by_aa.get(aa, 0) + 1

        if not total_by_aa:
            return 0.0

        # Calculate MSE between sequence distribution and host distribution
        errors = []
        for codon, count in counts.items():
            aa = table.forward_table[codon]
            seq_freq = count / total_by_aa[aa]

            # Host relative frequency
            synonymous = [c for c, a in table.forward_table.items() if a == aa]
            host_total = sum(self.codon_usage.get(c, 0) for c in synonymous)
            host_freq = self.codon_usage.get(codon, 0) / host_total if host_total > 0 else 0

            errors.append((seq_freq - host_freq) ** 2)

        if not errors:
            return 1.0

        avg_mse = sum(errors) / len(errors)
        return max(0.0, 1.0 - avg_mse * 5) # Scale to 0-1

    def _calculate_cai(self, dna_seq: str) -> float:
        table = CodonTable.unambiguous_dna_by_id.get(self.codon_table_id, CodonTable.standard_dna_table)
        log_weights_sum = 0.0
        valid_codon_count = 0
        for i in range(0, len(dna_seq), 3):
            codon = dna_seq[i:i+3].upper()
            if len(codon) != 3 or codon in table.stop_codons or codon not in self.codon_weights:
                continue
            weight = self.codon_weights.get(codon, 1e-9)
            if weight > 1e-10:
                log_weights_sum += log(weight)
                valid_codon_count += 1
        if valid_codon_count == 0:
            return 0.0
        return exp(log_weights_sum / valid_codon_count)

class GcContentSpecification(Specification):
    def __init__(self, target_gc: Union[float, Tuple[float, float]], window: Optional[float] = None):
        try:
            if isinstance(target_gc, (float, int)) and window is not None:
                target = max(0.0, min(1.0, float(target_gc)))
                self.min_gc = max(0.0, target - abs(float(window)))
                self.max_gc = min(1.0, target + abs(float(window)))
                self.target_repr = f"{target:.2f} (±{abs(float(window)):.2f})"
            elif isinstance(target_gc, (list, tuple)) and len(target_gc) == 2:
                self.min_gc = max(0.0, min(float(target_gc[0]), float(target_gc[1])))
                self.max_gc = min(1.0, max(float(target_gc[0]), float(target_gc[1])))
                self.target_repr = f"[{self.min_gc:.2f} - {self.max_gc:.2f}]"
            else:
                logger.warning("Неверный формат target_gc. Используется значение по умолчанию [0.45, 0.55].")
                self.min_gc, self.max_gc = 0.45, 0.55
                self.target_repr = "[0.45 - 0.55]"
        except (TypeError, ValueError) as e:
            logger.warning(f"Ошибка парсинга target_gc: {e}. Используется значение по умолчанию [0.45, 0.55].")
            self.min_gc, self.max_gc = 0.45, 0.55
            self.target_repr = "[0.45 - 0.55]"

    def evaluate(self, sequence: str) -> SpecEvaluation:
        if not sequence:
            return SpecEvaluation(0.0, f"GC: N/A (пустая последовательность), Цель: {self.target_repr}")
        gc_count = sum(sequence.upper().count(c) for c in 'GC')
        gc_content = gc_count / len(sequence) if len(sequence) > 0 else 0.0
        if self.min_gc <= gc_content <= self.max_gc:
            score = 1.0
            message = f"GC: {gc_content:.3f} (в пределах целевого диапазона {self.target_repr})"
        else:
            deviation = min(abs(gc_content - self.min_gc), abs(gc_content - self.max_gc))
            score = max(0.0, 1.0 - (deviation ** 2) / 0.1 ** 2)
            message = f"GC: {gc_content:.3f} (вне целевого диапазона {self.target_repr})"
        return SpecEvaluation(score, message)

class AvoidPatternSpecification(Specification):
    def __init__(self, patterns: Union[List[Pattern], List[str]]):
        self.patterns = [Pattern(motif) for motif in patterns] if patterns and isinstance(patterns[0], str) else patterns

    def evaluate(self, sequence: str) -> SpecEvaluation:
        if not sequence or not self.patterns:
            return SpecEvaluation(1.0, "Запрещенные мотивы не найдены")
        all_locations = []
        for pattern in self.patterns:
            all_locations.extend(pattern.find_matches(sequence))
        unique_locations = sorted(set(all_locations))
        score = max(0.0, 1.0 - 0.2 * len(unique_locations))
        msg = f"Найдено {len(unique_locations)} запрещенных мотивов" if unique_locations else "Запрещенные мотивы не найдены"
        return SpecEvaluation(score, msg, unique_locations)

class RnaFoldingSpecification(Specification):
    _mfe_cache = {}

    def __init__(self, window_size: int, ideal_mfe: float, bad_mfe: float):
        self.window_size = window_size
        self.ideal_mfe = ideal_mfe
        self.bad_mfe = min(ideal_mfe, bad_mfe)
        self.has_rnafold = self._check_rnafold()

    def _check_rnafold(self) -> bool:
        try:
            subprocess.run(["RNAfold", "--version"], capture_output=True, check=True)
            return True
        except (subprocess.SubprocessError, FileNotFoundError):
            logger.warning("RNAfold не найден. Спецификация RnaFolding отключена.")
            return False

    def evaluate(self, sequence: str) -> SpecEvaluation:
        if not self.has_rnafold:
            return SpecEvaluation(1.0, "RNAfold недоступен, спецификация отключена")
        if not sequence:
            return SpecEvaluation(0.0, "MFE: N/A (пустая последовательность)")
        if len(sequence) < self.window_size:
            return SpecEvaluation(0.0, f"MFE: Последовательность короче окна ({len(sequence)} < {self.window_size})")

        prefix = sequence[:self.window_size].upper()
        if prefix in self._mfe_cache:
            mfe = self._mfe_cache[prefix]
        else:
            try:
                process = subprocess.run(
                    ["RNAfold", "--noPS"], input=prefix, text=True, capture_output=True, check=True
                )
                output = process.stdout.strip().split('\n')
                mfe_line = output[1] if len(output) > 1 else ''
                mfe_match = re.search(r'\(\s*([-+]?\d*\.?\d+)\)', mfe_line)
                if mfe_match:
                    mfe = float(mfe_match.group(1))
                    self._mfe_cache[prefix] = mfe
                else:
                    raise ValueError(f"Не удалось распарсить MFE из вывода: {mfe_line}")
            except Exception as e:
                logger.warning(f"Ошибка выполнения RNAfold: {e}")
                return SpecEvaluation(0.0, f"MFE: Ошибка ({str(e)})")

        if mfe <= self.bad_mfe:
            score = 0.0
        elif mfe >= self.ideal_mfe:
            score = 1.0
        else:
            score = (mfe - self.bad_mfe) / (self.ideal_mfe - self.bad_mfe)
        return SpecEvaluation(score, f"MFE: {mfe:.2f} (цель: {self.ideal_mfe}, порог: {self.bad_mfe})")

class CodonPairBiasSpecification(Specification):
    def __init__(self, cps_table: Dict[str, float], sigmoid_k: float = 3.0, sigmoid_center: float = 0.0):
        self.cps_table = cps_table
        self.k = sigmoid_k
        self.center = sigmoid_center

    def evaluate(self, sequence: str) -> SpecEvaluation:
        if not sequence or len(sequence) < 6:
            return SpecEvaluation(0.0, "CPB: N/A (последовательность слишком короткая)")
        scores = []
        for i in range(0, len(sequence) - 6, 3):
            pair = sequence[i:i+6].upper()
            if len(pair) != 6:
                continue
            cps = self.cps_table.get(pair, 1e-9)
            score = 1.0 / (1.0 + exp(-self.k * (log(cps + 1e-9) - self.center)))
            scores.append(score)
        if not scores:
            return SpecEvaluation(0.0, "CPB: Нет валидных пар кодонов")
        avg_score = statistics.mean(scores)
        return SpecEvaluation(avg_score, f"CPB: {avg_score:.4f}")

class RbsSpecification(Specification):
    def __init__(self, sd_sequence: str = "AGGAGG", optimal_spacing: Tuple[int, int] = (5, 10)):
        self.sd_sequence = sd_sequence.upper()
        self.min_spacing, self.max_spacing = optimal_spacing
        self._regex = re.compile(f"(?={self.sd_sequence})")

    def get_name(self) -> str:
        return "RbsSpecification"

    def evaluate(self, sequence: str) -> SpecEvaluation:
        if len(sequence) < len(self.sd_sequence) + self.min_spacing + 3:
            return SpecEvaluation(0.0, "RBS: Последовательность слишком короткая для оценки RBS")
        prefix = sequence[:50]
        matches = [(m.start(), m.start() + len(self.sd_sequence)) for m in self._regex.finditer(prefix)]
        if not matches:
            return SpecEvaluation(0.0, "RBS: Последовательность Shine-Dalgarno не найдена")
        best_score = 0.0
        best_message = ""
        for start, end in matches:
            spacing = len(sequence) - end - 3
            if self.min_spacing <= spacing <= self.max_spacing:
                score = 1.0
                message = f"RBS: SD найдена на {start}-{end}, расстояние {spacing} (оптимально)"
            else:
                penalty = min(abs(spacing - self.min_spacing), abs(spacing - self.max_spacing)) / self.max_spacing
                score = max(0.0, 1.0 - penalty)
                message = f"RBS: SD найдена на {start}-{end}, расстояние {spacing} (не оптимально)"
            if score > best_score:
                best_score = score
                best_message = message
        return SpecEvaluation(best_score, best_message)

def _evaluate_individual(dna_sequence: str, aa_sequence: str, specifications: List[Specification], weights: Dict[str, float]) -> Tuple[str, float, Dict[str, float]]:
    """Standalone function for parallel evaluation of an individual."""
    if not dna_sequence or len(dna_sequence) != len(aa_sequence) * 3:
        return dna_sequence, -float('inf'), {}

    spec_scores = {}
    for spec in specifications:
        try:
            eval_result = spec.evaluate(dna_sequence)
            spec_scores[spec.get_name()] = eval_result.score
        except Exception as e:
            logger.warning(f"Ошибка оценки {spec.get_name()}: {e}")
            spec_scores[spec.get_name()] = 0.0

    total_score = sum(weights.get(name, 0.0) * score for name, score in spec_scores.items())
    return dna_sequence, total_score, spec_scores

# --- Оптимизатор кодонов ---
class CodonOptimizer:
    def __init__(self, aa_sequence: str, codon_usage: Dict[str, float], specifications: List[Specification],
                 weights: Dict[str, float], codon_table_id: int = 11, num_processes: Optional[int] = None):
        self.aa_sequence = aa_sequence.upper().replace('*', '')
        self.codon_usage = codon_usage
        self.specifications = specifications
        self.weights = weights
        self.codon_table_id = codon_table_id
        self.num_processes = num_processes or min(cpu_count(), 8)
        self.back_table = self._get_backtranslation_table()
        self.fitness_cache = OrderedDict()
        self.spec_cache = OrderedDict()
        self.cache_limit = 50000
        self._executor = None
        if not self.aa_sequence:
            raise ValueError("Пустая аминокислотная последовательность")
        self._validate_weights()
        for i, aa in enumerate(self.aa_sequence):
            if aa not in self.back_table:
                raise ValueError(f"Недопустимая аминокислота '{aa}' на позиции {i}")

    def __enter__(self):
        self._executor = ProcessPoolExecutor(max_workers=self.num_processes)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._executor:
            self._executor.shutdown()
            self._executor = None

    def _get_backtranslation_table(self) -> Dict[str, List[str]]:
        table = CodonTable.unambiguous_dna_by_id.get(self.codon_table_id, CodonTable.standard_dna_table)
        result = {}
        for codon, aa in table.forward_table.items():
            if aa != '*':
                result.setdefault(aa, []).append(codon)
        return result

    def _validate_weights(self):
        total_weight = 0.0
        logger.info("--- Веса оптимизатора ---")
        for spec in self.specifications:
            spec_name = spec.get_name()
            if spec_name not in self.weights:
                self.weights[spec_name] = 1.0
                logger.info(f"Вес по умолчанию 1.0 для {spec_name}")
            if self.weights[spec_name] < 0:
                logger.warning(f"Отрицательный вес для {spec_name}. Установлен в 0.")
                self.weights[spec_name] = 0.0
            status = "ОТКЛЮЧЕНО" if self.weights[spec_name] == 0 else f"{self.weights[spec_name]:.2f}"
            logger.info(f"  - {spec_name}: {status}")
            total_weight += self.weights[spec_name]
        if total_weight <= 0:
            logger.warning("Все веса нулевые или отрицательные. Оптимизация может быть неэффективной.")

    def _calculate_fitness(self, dna_sequence: str) -> float:
        """Sequential fitness calculation for a single sequence (with caching)."""
        if dna_sequence in self.fitness_cache:
            return self.fitness_cache[dna_sequence]

        _, total_score, spec_scores = _evaluate_individual(
            dna_sequence, self.aa_sequence, self.specifications, self.weights
        )

        self.fitness_cache[dna_sequence] = total_score
        self.spec_cache[dna_sequence] = spec_scores

        if len(self.fitness_cache) > self.cache_limit:
            self.fitness_cache.popitem(last=False)
            self.spec_cache.popitem(last=False)
        return total_score

    def _calculate_population_fitness(self, population: List[str]) -> List[float]:
        """Parallel fitness calculation for a whole population."""
        # Check cache first
        needed_indices = []
        fitness_scores = [0.0] * len(population)

        for i, ind in enumerate(population):
            if ind in self.fitness_cache:
                fitness_scores[i] = self.fitness_cache[ind]
            else:
                needed_indices.append(i)

        if not needed_indices:
            return fitness_scores

        needed_sequences = [population[i] for i in needed_indices]

        if self._executor:
            # Parallel evaluation
            eval_func = functools.partial(
                _evaluate_individual,
                aa_sequence=self.aa_sequence,
                specifications=self.specifications,
                weights=self.weights
            )
            results = list(self._executor.map(eval_func, needed_sequences))
        else:
            # Fallback to sequential
            results = [
                _evaluate_individual(seq, self.aa_sequence, self.specifications, self.weights)
                for seq in needed_sequences
            ]

        for idx, (seq, score, specs) in zip(needed_indices, results):
            fitness_scores[idx] = score
            self.fitness_cache[seq] = score
            self.spec_cache[seq] = specs

        return fitness_scores

    def _evaluate_final_sequence(self, dna_sequence: str) -> Dict:
        evaluations = {spec.get_name(): spec.evaluate(dna_sequence).as_dict() for spec in self.specifications}
        return {
            'aa_sequence': self.aa_sequence,
            'optimized_dna_sequence': dna_sequence,
            'fitness': self._calculate_fitness(dna_sequence),
            'gc_content': sum(dna_sequence.upper().count(c) for c in 'GC') / len(dna_sequence) if dna_sequence else 0.0,
            'cai': CodonUsageSpecification(self.codon_usage, self.codon_table_id).evaluate(dna_sequence).score,
            'spec_evaluations': evaluations
        }

    def optimize_greedy(self) -> Tuple[str, Dict]:
        logger.info("Запуск жадной оптимизации...")
        start_time = time.time()
        dna = ""
        for i, aa in enumerate(self.aa_sequence):
            codons = self.back_table.get(aa, [])
            if not codons:
                logger.warning(f"Нет кодонов для аминокислоты '{aa}' на позиции {i}")
                dna += "NNN"
                continue
            best_codon = max(
                codons,
                key=lambda c: sum(
                    self.weights.get(spec.get_name(), 0.0) * spec.evaluate(dna + c).score
                    for spec in self.specifications
                ),
                default=codons[0]
            )
            dna += best_codon
        metrics = self._evaluate_final_sequence(dna)
        logger.info(f"Жадная оптимизация завершена за {time.time() - start_time:.2f} секунд")
        return dna, metrics

    def optimize_ma(self, population_size: int, num_generations: int, crossover_rate: float,
                    mutation_rate: float, tournament_size: int, elitism_count: int,
                    local_search_steps: int) -> Tuple[str, Dict]:
        logger.info("Запуск меметического алгоритма...")
        logger.info(f"Параметры: Pop={population_size}, Gens={num_generations}, XRate={crossover_rate:.2f}, "
                    f"MutRate={mutation_rate:.3f}, TournSize={tournament_size}, Elite={elitism_count}, LS_Steps={local_search_steps}")
        start_time = time.time()

        population = [self._random_individual() for _ in range(population_size - 1)]
        greedy_dna, _ = self.optimize_greedy()
        population.append(greedy_dna)

        population_fitness = self._calculate_population_fitness(population)
        best_dna, best_fitness = max(zip(population, population_fitness), key=lambda x: x[1])
        logger.info(f"Начальная лучшая приспособленность: {best_fitness:.4f}")

        for gen in range(num_generations):
            next_population = []
            # Elite stay
            elite_indices = sorted(range(len(population_fitness)), key=lambda i: population_fitness[i], reverse=True)[:elitism_count]
            next_population.extend([population[i] for i in elite_indices])

            while len(next_population) < population_size:
                parent1, parent2 = self._select_parents(population, population_fitness, tournament_size)
                if random.random() < crossover_rate:
                    child1, child2 = self._crossover(parent1, parent2)
                else:
                    child1, child2 = parent1, parent2

                child1 = self._mutate(child1, mutation_rate)
                child2 = self._mutate(child2, mutation_rate)

                effective_ls_steps = 0 if len(self.aa_sequence) < 10 else local_search_steps
                if effective_ls_steps > 0:
                    child1 = self._local_search(child1, effective_ls_steps)
                    if len(next_population) < population_size - 1:
                        child2 = self._local_search(child2, effective_ls_steps)

                next_population.extend([child1, child2][:population_size - len(next_population)])

            population = next_population
            population_fitness = self._calculate_population_fitness(population)

            current_best_idx = max(range(len(population_fitness)), key=lambda i: population_fitness[i])
            if population_fitness[current_best_idx] > best_fitness:
                best_dna = population[current_best_idx]
                best_fitness = population_fitness[current_best_idx]

            if (gen + 1) % max(1, num_generations // 5) == 0:
                logger.info(f"Поколение: {gen+1}/{num_generations}, Лучшая: {best_fitness:.4f}, Средняя: {statistics.mean(population_fitness):.4f}")

        metrics = self._evaluate_final_sequence(best_dna)
        logger.info(f"Меметический алгоритм завершен за {time.time() - start_time:.2f} секунд")
        return best_dna, metrics

    def _random_individual(self) -> str:
        return "".join(random.choice(self.back_table.get(aa, ['NNN'])) for aa in self.aa_sequence)

    def _calculate_diversity(self, ind: str, population: List[str]) -> float:
        """Calculates average Hamming distance (at codon level) from individual to population."""
        if not population: return 0.0
        total_dist = 0
        num_codons = len(ind) // 3
        if num_codons == 0: return 0.0

        # Sample population if it's too large for performance
        sample_size = min(len(population), 20)
        sample = random.sample(population, sample_size)

        for other in sample:
            if len(ind) != len(other): continue
            dist = sum(1 for i in range(0, len(ind), 3) if ind[i:i+3] != other[i:i+3])
            total_dist += dist / num_codons
        return total_dist / len(sample)

    def _select_parents(self, population: List[str], fitness_scores: List[float], tournament_size: int) -> Tuple[str, str]:
        """Tournament selection with diversity awareness."""
        def select_one():
            competitors = random.sample(range(len(population)), min(tournament_size, len(population)))
            diversity_scores = []
            for i in competitors:
                diversity = self._calculate_diversity(population[i], population)
                # Boost fitness of diverse individuals
                adjusted_fitness = fitness_scores[i] * (1 + 0.2 * diversity)
                diversity_scores.append((i, adjusted_fitness))
            return population[max(diversity_scores, key=lambda x: x[1])[0]]

        parent1 = select_one()
        parent2 = select_one()
        # Try to find a different parent for crossover
        for _ in range(5):
            if parent2 != parent1: break
            parent2 = select_one()
        return parent1, parent2

    def _crossover(self, parent1: str, parent2: str) -> Tuple[str, str]:
        if len(self.aa_sequence) < 2:
            return parent1, parent2
        cut = random.randint(1, len(self.aa_sequence) - 1) * 3
        return parent1[:cut] + parent2[cut:], parent2[:cut] + parent1[cut:]

    def _mutate(self, dna: str, mutation_rate: float) -> str:
        codons = [dna[i:i+3] for i in range(0, len(dna), 3)]
        for i, aa in enumerate(self.aa_sequence):
            if random.random() < mutation_rate:
                alternatives = [c for c in self.back_table.get(aa, []) if c != codons[i]]
                if alternatives:
                    codons[i] = random.choice(alternatives)
        return "".join(codons)

    def _local_search(self, dna: str, max_steps: int) -> str:
        """Improved local search targeting specific problem areas."""
        current_dna = dna
        current_fitness = self._calculate_fitness(dna)

        # Get problem areas from specifications
        spec_evals = self.spec_cache.get(dna, {})
        problem_indices = set()

        # Add indices where patterns are found
        avoid_spec = next((s for s in self.specifications if isinstance(s, AvoidPatternSpecification)), None)
        if avoid_spec:
            eval_res = avoid_spec.evaluate(dna)
            for start, end in eval_res.locations:
                # Convert nucleotide indices to codon indices
                for i in range(start // 3, (end + 2) // 3):
                    if i < len(self.aa_sequence):
                        problem_indices.add(i)

        # Add 5' end if RNA folding is poor
        rna_spec = next((s for s in self.specifications if isinstance(s, RnaFoldingSpecification)), None)
        if rna_spec and spec_evals.get("RnaFolding", 1.0) < 0.8:
            for i in range(min(len(self.aa_sequence), rna_spec.window_size // 3 + 1)):
                problem_indices.add(i)

        # Fallback to random mutable indices if no specific problems
        mutable_indices = [i for i, aa in enumerate(self.aa_sequence) if len(self.back_table.get(aa, [])) > 1]
        if not problem_indices:
            problem_indices = set(random.sample(mutable_indices, min(len(mutable_indices), max_steps * 2)))

        # Limit problem indices to mutable ones
        problem_indices = [idx for idx in problem_indices if idx in mutable_indices]
        random.shuffle(problem_indices)

        for step, idx in enumerate(problem_indices):
            if step >= max_steps:
                break

            original_codon = current_dna[idx*3:idx*3+3]
            alternatives = [c for c in self.back_table.get(self.aa_sequence[idx], []) if c != original_codon]

            best_local_dna = current_dna
            best_local_fitness = current_fitness

            dna_list = list(current_dna)
            for new_codon in alternatives:
                dna_list[idx*3:idx*3+3] = list(new_codon)
                temp_dna = "".join(dna_list)
                fitness = self._calculate_fitness(temp_dna)
                if fitness > best_local_fitness:
                    best_local_fitness = fitness
                    best_local_dna = temp_dna

            if best_local_dna != current_dna:
                current_dna = best_local_dna
                current_fitness = best_local_fitness

        return current_dna

# --- Менеджер вывода ---
class OutputManager:
    def __init__(self, output_dir: str = OUTPUT_DIR):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

    def save_fasta(self, seq: str, description: str, filename: str):
        record = SeqRecord(Seq(seq), id=os.path.basename(filename).split('.')[0], description=description)
        filepath = os.path.join(self.output_dir, filename)
        try:
            SeqIO.write(record, filepath, "fasta")
            logger.info(f"Последовательность сохранена по пути {filepath}")
        except Exception as e:
            logger.error(f"Ошибка сохранения FASTA {filepath}: {e}")

    def save_metrics_plot(self, metrics: Dict, title: str, filename: str, show: bool = False):
        if "spec_evaluations" not in metrics:
            return
        labels, scores = zip(*sorted((k, v["score"]) for k, v in metrics["spec_evaluations"].items()))
        plt.figure(figsize=(7, 5))
        plt.bar(labels, scores, color="skyblue")
        plt.ylim(0, 1.1)
        plt.title(title)
        plt.ylabel("Оценка")
        plt.xticks(rotation=25, ha='right')
        plt.tight_layout()
        filepath = os.path.join(self.output_dir, filename)
        try:
            plt.savefig(filepath)
            logger.info(f"График метрик сохранен по пути {filepath}")
            if show:
                plt.show()
        except Exception as e:
            logger.error(f"Ошибка сохранения графика {filepath}: {e}")
        plt.close()

    def save_results_summary(self, results: Dict, filename: str = "all_optimization_results_summary.json"):
        filepath = os.path.join(self.output_dir, filename)
        try:
            with open(filepath, 'w') as f:
                json.dump(results, f, indent=2)
            logger.info(f"Сводка результатов сохранена по пути {filepath}")
        except Exception as e:
            logger.error(f"Ошибка сохранения сводки {filepath}: {e}")

# --- Основное приложение ---
class CodonOptimizationApp:
    def __init__(self, organism: str = FALLBACK_ORGANISM, organism_id: str = FALLBACK_ORGANISM_ID,
                 codon_table_id: Optional[int] = None):
        self.config = ConfigurationManager().load_config()
        self.codon_table_id = codon_table_id if codon_table_id is not None else self.config["default_codon_table"]
        self.resources = ResourceManager(organism=organism, organism_id=organism_id, codon_table_id=self.codon_table_id)
        self.output = OutputManager()
        self.test_sequences = {
            "ShortPeptide": "MAGWSR",
            "HemoglobinSubunit": "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFDSFGNLSSPTAILGNPM",
            "HypotheticalProtein1": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVFVPP",
            "HypotheticalProtein2": "MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDALPNISDAERIFAELLTGLAAAQPGFPLAQLKTFVDQEFAQIKHVLHGISLLGQCPDSINAALICRGEKMSIAIMAGVLEARGHNVTVIDPVEKLLAVGHYLLE",
            "GFP": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
        }

    def setup_specifications(self, codon_usage: Dict[str, float], cps_table: Optional[Dict[str, float]], mode: str = "optimization") -> List[Specification]:
        specs = [
            CodonUsageSpecification(codon_usage, self.codon_table_id, mode=mode),
            GcContentSpecification(self.config["target_gc_range"]),
            AvoidPatternSpecification(self.config["avoid_motifs"]),
            RbsSpecification(),
            RnaFoldingSpecification(**self.config["rna_folding"])
        ]
        if cps_table:
            specs.append(CodonPairBiasSpecification(cps_table, **self.config["cpb_parameters"]))
        else:
            logger.info("Отключение CodonPairBias (таблица CPS отсутствует)")
            self.config["spec_weights"] = {k: v for k, v in self.config["spec_weights"].items() if k != "CodonPairBias"}
        return specs

    def log_resources(self):
        cpu_percent = psutil.cpu_percent()
        memory = psutil.virtual_memory()
        logger.info(f"Использование ресурсов: CPU: {cpu_percent}%, Память: {memory.percent}%")

    def run(self, show_plots: bool = False, mode: str = "optimization"):
        logger.info(f"--- Скрипт оптимизации кодонов ({mode}) ---")
        logger.info(f"Оптимизация выполняется для организма: {self.resources.organism} (TaxID: {self.resources.organism_id})")
        logger.info(f"Используется таблица кодонов: {self.codon_table_id}")
        try:
            self.log_resources()
            codon_usage = self.resources.load_codon_usage()
            cps_table = self.resources.calculate_codon_pair_scores()
            specifications = self.setup_specifications(codon_usage, cps_table, mode=mode)
            results = {}
            for i, (name, aa_seq) in enumerate(self.test_sequences.items(), 1):
                logger.info(f"--- Обработка последовательности {i}/{len(self.test_sequences)}: '{name}' ---")
                logger.info(f"Аминокислотная последовательность ({len(aa_seq)} aa): {aa_seq[:60]}{'...' if len(aa_seq) > 60 else ''}")
                try:
                    with CodonOptimizer(
                        aa_sequence=aa_seq,
                        codon_usage=codon_usage,
                        specifications=specifications,
                        weights=self.config["spec_weights"],
                        codon_table_id=self.codon_table_id
                    ) as optimizer:
                        dna, metrics = optimizer.optimize_ma(**self.config["ma_parameters"])
                        logger.info(f"Метрики: Приспособленность={metrics['fitness']:.4f}, GC={metrics['gc_content']:.3f}, CAI={metrics['cai']:.4f}")
                        self.output.save_fasta(dna, f"Оптимизированная CDS для {name}", f"optimized_{name}.fasta")
                        self.output.save_metrics_plot(metrics, f"Метрики ({name})", f"metrics_{name}.png", show=show_plots)
                        results[name] = metrics
                except Exception as e:
                    logger.error(f"Ошибка оптимизации '{name}': {e}")
            self.output.save_results_summary(results)
            self.log_resources()
            logger.info("--- Оптимизация завершена ---")
        except Exception as e:
            logger.error(f"Критическая ошибка: {e}")
            raise

if __name__ == "__main__":
    # Установите ваш email для Entrez
    if "ENTREZ_EMAIL" not in os.environ:
        os.environ["ENTREZ_EMAIL"] = "your_email@example.com"

    # Можно выбрать режим: "optimization" (максимизация CAI) или "harmonization" (соответствие распределению)
    mode = "optimization"

    app = CodonOptimizationApp(organism="Bacillus subtilis", organism_id="1423", codon_table_id=11)
    app.run(show_plots=False, mode=mode)
