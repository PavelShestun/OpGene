from fpdf import FPDF
import datetime

class ReportGenerator:
    @staticmethod
    def create_pdf(org_name, aa_seq, dna_seq, metrics, fitness):
        pdf = FPDF()
        pdf.add_page()
        
        # Header
        pdf.set_font("Arial", "B", 16)
        pdf.cell(190, 10, "OpGene Optimization Report", ln=True, align="C")
        pdf.set_font("Arial", "", 10)
        pdf.cell(190, 10, f"Generated on: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}", ln=True, align="C")
        pdf.ln(10)

        # Project Info
        pdf.set_fill_color(240, 240, 240)
        pdf.set_font("Arial", "B", 12)
        pdf.cell(190, 8, "Project Summary", ln=True, fill=True)
        pdf.set_font("Arial", "", 10)
        pdf.cell(95, 8, f"Target Organism: {org_name}")
        pdf.cell(95, 8, f"Final Fitness Score: {fitness:.2f}", ln=True)
        pdf.ln(5)

        # Metrics Table
        pdf.set_font("Arial", "B", 11)
        pdf.cell(100, 8, "Objective", border=1)
        pdf.cell(90, 8, "Status/Value", border=1, ln=True)
        
        pdf.set_font("Arial", "", 10)
        for obj, val in metrics.items():
            pdf.cell(100, 8, str(obj), border=1)
            pdf.cell(90, 8, str(val), border=1, ln=True)
        pdf.ln(10)

        # Sequences
        pdf.set_font("Arial", "B", 12)
        pdf.cell(190, 8, "Amino Acid Sequence", ln=True, fill=True)
        pdf.set_font("Courier", "", 8)
        pdf.multi_cell(190, 5, aa_seq)
        pdf.ln(5)

        pdf.set_font("Arial", "B", 12)
        pdf.cell(190, 8, "Optimized DNA Sequence", ln=True, fill=True)
        pdf.set_font("Courier", "", 8)
        pdf.multi_cell(190, 5, dna_seq)
        
        return bytes(pdf.output())
