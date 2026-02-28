from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from io import StringIO

class GenBankExporter:
    @staticmethod
    def create_record(dna_seq, org_name, aa_seq, profile):
        seq_obj = Seq(dna_seq)
        record = SeqRecord(
            seq_obj,
            id="OpGene_Opt",
            name="Optimized",
            description=f"Codon optimized for {org_name}",
            annotations={"molecule_type": "DNA"}
        )

        # 1. CDS Feature
        cds_feature = SeqFeature(
            FeatureLocation(0, len(dna_seq)),
            type="CDS",
            qualifiers={
                "translation": aa_seq,
                "product": "Optimized protein",
                "note": f"Host: {org_name}"
            }
        )
        record.features.append(cds_feature)

        # 2. Ramp Zone Feature
        if profile.use_ramp_hypothesis:
            ramp_feature = SeqFeature(
                FeatureLocation(0, 45),
                type="misc_feature",
                qualifiers={"note": "Translation Ramp (Controlled speed zone)"}
            )
            record.features.append(ramp_feature)

        output = StringIO()
        SeqIO.write(record, output, "genbank")
        return output.getvalue()
