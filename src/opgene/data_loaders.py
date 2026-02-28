# src/opgene/data_loaders.py
import os
import json
import logging
import math
from Bio import Entrez, SeqIO

logger = logging.getLogger(__name__)

class DataLoader:
    def __init__(self, email: str, data_dir="data"):
        self.email = email
        self.data_dir = data_dir
        Entrez.email = email
        os.makedirs(data_dir, exist_ok=True)

    def load_organism_data(self, organism_name: str, tax_id: str) -> dict:
        """Загружает частоты кодонов и CPS таблицу"""
        safe_name = organism_name.replace(' ', '_')
        filename = os.path.join(self.data_dir, f"{safe_name}_v2.json")
        
        if os.path.exists(filename):
            with open(filename, 'r') as f:
                return json.load(f)

        logger.info(f"Downloading genomic data for {organism_name}...")
        try:
            handle = Entrez.esearch(db="nucleotide", term=f"txid{tax_id}[Organism] AND RefSeq[filter]", retmax=30)
            ids = Entrez.read(handle)["IdList"]
            
            codon_counts = {}
            pair_counts = {}
            total_codons = 0
            total_pairs = 0
            
            handle = Entrez.efetch(db="nucleotide", id=",".join(ids), rettype="fasta_cds_na", retmode="text")
            for rec in SeqIO.parse(handle, "fasta"):
                seq = str(rec.seq).upper()
                if len(seq) < 9 or len(seq) % 3 != 0: continue
                
                codons = [seq[i:i+3] for i in range(0, len(seq), 3) if "N" not in seq[i:i+3]]
                for i in range(len(codons)):
                    c1 = codons[i]
                    codon_counts[c1] = codon_counts.get(c1, 0) + 1
                    total_codons += 1
                    
                    if i < len(codons) - 1:
                        pair = c1 + codons[i+1]
                        pair_counts[pair] = pair_counts.get(pair, 0) + 1
                        total_pairs += 1

            # Расчет частот и CPS
            usage = {k: v/total_codons for k, v in codon_counts.items()}
            cps = {}
            for pair, obs_count in pair_counts.items():
                c1, c2 = pair[:3], pair[3:]
                f_obs = obs_count / total_pairs
                f_exp = usage[c1] * usage[c2]
                cps[pair] = math.log(f_obs / f_exp) if f_exp > 0 else 0

            data = {"usage": usage, "cps": cps}
            with open(filename, 'w') as f:
                json.dump(data, f)
            return data
        except Exception as e:
            logger.error(f"Failed to fetch data: {e}")
            return {"usage": {}, "cps": {}}

    def _get_fallback_table(self):
        # Полная таблица E. coli K-12 MG1655
        # Это предотвратит CAI = 1.0 и ошибку Ramp too fast
        return {
            'TTT': 0.58, 'TTC': 0.42, 'TTA': 0.14, 'TTG': 0.13, 
            'CTT': 0.12, 'CTC': 0.10, 'CTA': 0.04, 'CTG': 0.47, 
            'ATT': 0.49, 'ATC': 0.39, 'ATA': 0.11, 'ATG': 1.00,
            'GTT': 0.28, 'GTC': 0.20, 'GTA': 0.17, 'GTG': 0.35, 
            'TCT': 0.17, 'TCC': 0.15, 'TCA': 0.14, 'TCG': 0.14, 
            'CCT': 0.18, 'CCC': 0.13, 'CCA': 0.20, 'CCG': 0.49,
            'ACT': 0.19, 'ACC': 0.40, 'ACA': 0.17, 'ACG': 0.25, 
            'GCT': 0.18, 'GCC': 0.26, 'GCA': 0.23, 'GCG': 0.33, 
            'TAT': 0.59, 'TAC': 0.41, 'CAT': 0.57, 'CAC': 0.43,
            'CAA': 0.34, 'CAG': 0.66, 'AAT': 0.49, 'AAC': 0.51, 
            'AAA': 0.74, 'AAG': 0.26, 'GAT': 0.63, 'GAC': 0.37, 
            'GAA': 0.68, 'GAG': 0.32, 'TGT': 0.46, 'TGC': 0.54,
            'TGG': 1.00, 'CGT': 0.36, 'CGC': 0.36, 'CGA': 0.07, 
            'CGG': 0.11, 'AGT': 0.16, 'AGC': 0.25, 'AGA': 0.07, 
            'AGG': 0.04, 'GGT': 0.35, 'GGC': 0.37, 'GGA': 0.13,
            'GGG': 0.15
        }
