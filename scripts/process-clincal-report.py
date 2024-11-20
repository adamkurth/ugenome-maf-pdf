import pandas as pd
import json
from pathlib import Path
import logging
from datetime import datetime
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, Image
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.barcharts import VerticalBarChart
from reportlab.graphics.charts.piecharts import Pie
import matplotlib.pyplot as plt
import seaborn as sns
import io

class ClinicalReportGenerator:
    def __init__(self):
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)
        self.styles = getSampleStyleSheet()
        
        # Custom styles
        self.styles.add(ParagraphStyle(
            name='CustomTitle',
            parent=self.styles['Heading1'],
            fontSize=24,
            spaceAfter=30,
            textColor=colors.HexColor('#2C3E50')
        ))
        
        self.styles.add(ParagraphStyle(
            name='SectionHeader',
            parent=self.styles['Heading2'],
            fontSize=16,
            spaceBefore=20,
            spaceAfter=12,
            textColor=colors.HexColor('#34495E')
        ))
        
        self.styles.add(ParagraphStyle(
            name='SubSection',
            parent=self.styles['Heading3'],
            fontSize=14,
            spaceBefore=15,
            spaceAfter=10,
            textColor=colors.HexColor('#7F8C8D')
        ))

    def load_data(self, results_dir: Path) -> dict:
        """Load analysis results and processed data"""
        try:
            # Load complete analysis
            with open(results_dir / 'complete_analysis.json', 'r') as f:
                analysis_data = json.load(f)
            
            # Load processed data files
            processed_data = {}
            processed_dir = results_dir / 'processed_data'
            for file_path in processed_dir.glob('*.csv'):
                sample_name = file_path.stem.split('_')[0]
                if sample_name not in processed_data:
                    processed_data[sample_name] = {}
                variant_type = file_path.stem.split('_')[1]
                processed_data[sample_name][variant_type] = pd.read_csv(file_path)
            
            return {'analysis': analysis_data, 'processed': processed_data}
            
        except Exception as e:
            self.logger.error(f"Error loading data: {str(e)}")
            raise

    def create_variant_summary_table(self, sample_data: dict) -> Table:
        """Create a summary table of variant counts"""
        data = [['Variant Type', 'Count', 'Percentage']]
        
        variant_counts = sample_data['variant_metrics']['variant_counts']['by_type']
        total = sum(variant_counts.values())
        
        for variant_type, count in variant_counts.items():
            percentage = (count / total * 100) if total > 0 else 0
            data.append([
                variant_type,
                str(count),
                f"{percentage:.1f}%"
            ])
            
        table = Table(data, colWidths=[2*inch, 1.5*inch, 1.5*inch])
        table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#34495E')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 12),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.whitesmoke),
            ('TEXTCOLOR', (0, 1), (-1, -1), colors.black),
            ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
            ('FONTSIZE', (0, 1), (-1, -1), 10),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('BOX', (0, 0), (-1, -1), 2, colors.black),
            ('LINEBELOW', (0, 0), (-1, 0), 2, colors.black),
        ]))
        
        return table

    def create_impact_chart(self, sample_data: dict) -> Drawing:
        """Create a pie chart of variant impact distribution"""
        impact_data = sample_data['variant_metrics'].get('impact_distribution', {})
        
        drawing = Drawing(400, 200)
        pie = Pie()
        pie.x = 100
        pie.y = 25
        pie.width = 200
        pie.height = 150
        
        # Prepare data
        pie.data = list(impact_data.values())
        pie.labels = list(impact_data.keys())
        
        # Style
        pie.slices.strokeWidth = 0.5
        colors_list = [colors.HexColor('#3498DB'), colors.HexColor('#E74C3C'),
                      colors.HexColor('#2ECC71'), colors.HexColor('#F1C40F')]
        for i, slice in enumerate(pie.slices):
            slice.fillColor = colors_list[i % len(colors_list)]
        
        drawing.add(pie)
        return drawing

    def generate_gene_summary(self, sample_data: dict) -> list:
        """Generate gene summary content"""
        elements = []
        
        # Most mutated genes
        elements.append(Paragraph('Most Frequently Mutated Genes', self.styles['SubSection']))
        
        data = [['Gene', 'Mutation Count']]
        most_mutated = list(sample_data['gene_summary']['most_mutated_genes'].items())[:10]
        for gene, count in most_mutated:
            data.append([gene, str(count)])
        
        table = Table(data, colWidths=[2*inch, 1.5*inch])
        table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#34495E')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 12),
            ('GRID', (0, 0), (-1, -1), 1, colors.black)
        ]))
        
        elements.append(table)
        elements.append(Spacer(1, 20))
        
        return elements

    def create_clinical_report(self, sample_name: str, data: dict, output_dir: Path):
        """Generate complete clinical report"""
        doc = SimpleDocTemplate(
            output_dir / f"{sample_name}_clinical_report.pdf",
            pagesize=letter,
            rightMargin=72,
            leftMargin=72,
            topMargin=72,
            bottomMargin=72
        )
        
        elements = []
        
        # Title
        elements.append(Paragraph(f"Clinical Genomic Analysis Report", self.styles['CustomTitle']))
        elements.append(Paragraph(f"Sample ID: {sample_name}", self.styles['SectionHeader']))
        elements.append(Paragraph(f"Report Date: {datetime.now().strftime('%Y-%m-%d')}", self.styles['Normal']))
        elements.append(Spacer(1, 30))
        
        # Sample Information
        elements.append(Paragraph("Sample Overview", self.styles['SectionHeader']))
        sample_data = data['analysis']['sample_analyses'][sample_name]
        elements.append(Paragraph(
            f"Total Variants: {sample_data['file_stats']['total_variants']}", 
            self.styles['Normal']
        ))
        elements.append(Spacer(1, 20))
        
        # Variant Summary
        elements.append(Paragraph("Variant Summary", self.styles['SectionHeader']))
        elements.append(self.create_variant_summary_table(sample_data))
        elements.append(Spacer(1, 20))
        
        # Impact Distribution
        elements.append(Paragraph("Variant Impact Distribution", self.styles['SectionHeader']))
        elements.append(self.create_impact_chart(sample_data))
        elements.append(Spacer(1, 20))
        
        # Gene Analysis
        elements.append(Paragraph("Gene Analysis", self.styles['SectionHeader']))
        elements.extend(self.generate_gene_summary(sample_data))
        
        # Clinical Significance
        if 'clinical_significance' in sample_data['variant_metrics']:
            elements.append(PageBreak())
            elements.append(Paragraph("Clinical Significance", self.styles['SectionHeader']))
            clin_sig_data = sample_data['variant_metrics']['clinical_significance']
            for significance, count in clin_sig_data.items():
                elements.append(Paragraph(
                    f"{significance}: {count} variants",
                    self.styles['Normal']
                ))
        
        # Build the PDF
        doc.build(elements)

def main():
    # Setup
    results_dir = Path('data/results')
    output_dir = results_dir / 'reports'
    output_dir.mkdir(exist_ok=True)
    
    # Initialize report generator
    generator = ClinicalReportGenerator()
    
    try:
        # Load data
        data = generator.load_data(results_dir)
        
        # Generate report for each sample
        for sample_name in data['analysis']['sample_analyses'].keys():
            generator.create_clinical_report(sample_name, data, output_dir)
            logging.info(f"Generated clinical report for {sample_name}")
            
    except Exception as e:
        logging.error(f"Error generating reports: {str(e)}")
        raise

if __name__ == "__main__":
    main()