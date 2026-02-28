# src/opgene/factory.py
from .models import OrganismProfile, OrganismType
from .data_loaders import DataLoader

class OrganismFactory:
    def __init__(self, email: str):
        self.loader = DataLoader(email)

    @staticmethod
    def create(name: str, tax_id: str, email: str = "test@ex.com") -> OrganismProfile:
        name_lower = name.lower()
        if any(x in name_lower for x in ["coli", "subtilis", "bacterium"]):
            otype = OrganismType.PROKARYOTE
        elif any(x in name_lower for x in ["sapiens", "human", "mus", "cho", "hamster"]):
            otype = OrganismType.MAMMALIAN
        else:
            otype = OrganismType.EUKARYOTE
            
        loader = DataLoader(email)
        data = loader.load_organism_data(name, tax_id)
        
        return OrganismProfile(
            name=name, 
            tax_id=tax_id, 
            org_type=otype, 
            codon_usage=data["usage"],
            cps_table=data["cps"] # ТЕПЕРЬ ПЕРЕДАЕТСЯ
        )
