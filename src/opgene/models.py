from dataclasses import dataclass, field
from typing import Dict, List, Optional
from enum import Enum

class OrganismType(Enum):
    PROKARYOTE = "prokaryote"
    EUKARYOTE = "eukaryote"
    MAMMALIAN = "mammalian"

@dataclass
class OrganismProfile:
    name: str
    tax_id: str
    org_type: OrganismType
    codon_usage: Dict[str, float]
    cps_table: Dict[str, float] = field(default_factory=dict) # ДОБАВЛЕНО
    codon_table_id: int = 11
    use_ramp_hypothesis: bool = True
    suppress_cpg: bool = False
    ideal_gc_range: tuple = (0.4, 0.6)
    forbidden_motifs: List[str] = field(default_factory=list)

    def __post_init__(self):
        if self.org_type == OrganismType.MAMMALIAN:
            self.suppress_cpg = True
            self.ideal_gc_range = (0.45, 0.65)
        elif self.org_type == OrganismType.PROKARYOTE:
            self.use_ramp_hypothesis = True
            self.ideal_gc_range = (0.35, 0.60)
