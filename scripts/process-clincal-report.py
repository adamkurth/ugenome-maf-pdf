#!/usr/bin/env python3
import pandas as pd
import json
from pathlib import Path
import logging
from datetime import datetime
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table, 
                              TableStyle, PageBreak)

class ClinicalReportGenerator:
    def __init__(self):
        logging.basicConfig(level=logging.INFO,
                          format='%(asctime)s - %(levelname)s - %(message)s')
        self.logger = logging.getLogger(__name__)
        self.styles = getSampleStyleSheet()
        
        # Custom styles setup
        self.setup_styles()

    def setup_styles(self):
        """Set up all custom styles for the report"""
        self.styles.add(ParagraphStyle(
            name='ReportTitle',
            parent=self.styles['Heading1'],
            fontSize=14,
            spaceAfter=20,
            textColor=colors.black
        ))
        
        self.styles.add(ParagraphStyle(
            name='SectionHeader',
            parent=self.styles['Heading2'],
            fontSize=12,
            spaceBefore=15,
            spaceAfter=10,
            textColor=colors.black
        ))
        
        self.styles.add(ParagraphStyle(
            name='TableHeader',
            parent=self.styles['Normal'],
            fontSize=10,
            textColor=colors.black,
            alignment=1
        ))

    def load_processed_data(self, results_dir: Path, sample_name: str) -> dict:
        """Load processed data for the given sample"""
        processed_dir = results_dir / 'processed_data'
        data = {}
        
        try:
            # Load complete data
            complete_file = processed_dir / f'{sample_name}_complete.csv'
            if complete_file.exists():
                self.logger.info(f"Loading complete data from {complete_file}")
                data['complete'] = pd.read_csv(complete_file)
            else:
                self.logger.warning(f"Complete data file not found: {complete_file}")
            
            # Load variant-specific data
            for variant_type in ['snp', 'del', 'ins']:
                variant_file = processed_dir / f'{sample_name}_{variant_type}.csv'
                if variant_file.exists():
                    self.logger.info(f"Loading {variant_type} data from {variant_file}")
                    data[variant_type] = pd.read_csv(variant_file)
                else:
                    self.logger.info(f"No {variant_type} data file found: {variant_file}")
            
            if not data:
                raise FileNotFoundError(f"No data files found for sample {sample_name}")
                
            return data
            
        except Exception as e:
            self.logger.error(f"Error loading data for {sample_name}: {str(e)}")
            raise

    def create_header_section(self) -> list:
        """Create the report header section"""
        elements = []
        
        report_id = f"MAF-{datetime.now().strftime('%Y%m%d')}"
        
        # Header table data
        header_data = [
            ['Report N°', report_id],
            ['Date', datetime.now().strftime('%Y-%m-%d')]
        ]
        
        header_table = Table(header_data, colWidths=[2*inch, 4*inch])
        header_table.setStyle(TableStyle([
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, -1), 'Helvetica'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 12),
        ]))
        
        elements.append(header_table)
        elements.append(Spacer(1, 20))
        
        return elements

    def create_patient_section(self, sample_name: str) -> list:
        """Create the patient information section"""
        elements = []
        
        patient_data = [
            ['PATIENT', ''],
            ['ID N°', 'XXXXX'],
            ['Clinical Diagnosis', 'mCRC'],
            ['Clinical Trial', '_____________']
        ]
        
        patient_table = Table(patient_data, colWidths=[2*inch, 4*inch])
        patient_table.setStyle(TableStyle([
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, -1), 'Helvetica'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 12),
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#E8E8E8'))
        ]))
        
        elements.append(patient_table)
        elements.append(Spacer(1, 20))
        
        return elements

    def create_variant_table(self, df: pd.DataFrame, variant_type: str) -> Table:
        """Create a table for variants"""
        # Select relevant columns
        columns = ['Hugo_Symbol', 'Chromosome', 'Start_Position', 
                  'Reference_Allele', 'Tumor_Seq_Allele2', 'HGVSp_Short']
        
        # Filter and prepare data
        table_data = [columns]  # Header row
        for _, row in df.iterrows():
            table_data.append([str(row[col]) for col in columns])
            
        # Create table
        table = Table(table_data, repeatRows=1)
        table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#F5F5F5')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.black),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 8),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('TEXTCOLOR', (0, 1), (-1, -1), colors.black),
            ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
            ('FONTSIZE', (0, 1), (-1, -1), 8),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('BOX', (0,0), (-1,-1), 2, colors.black),
            ('LINEBELOW', (0,0), (-1,0), 2, colors.black),
        ]))
        
        return table

    def generate_report(self, sample_name: str, data: dict, output_dir: Path):
        """Generate the complete clinical report"""
        doc = SimpleDocTemplate(
            str(output_dir / f"{sample_name}_clinical_report.pdf"),
            pagesize=letter,
            rightMargin=72,
            leftMargin=72,
            topMargin=72,
            bottomMargin=72
        )
        
        elements = []
        
        # Add header section
        elements.extend(self.create_header_section())
        
        # Add patient section
        elements.extend(self.create_patient_section(sample_name))
        
        # Add assay section
        elements.append(Paragraph("ASSAY", self.styles['SectionHeader']))
        
        complete_df = data['complete']
        total_variants = len(complete_df)
        
        assay_data = [
            ['Genomic Target', 'WES', 'Sequencer', 'HiSeq2000'],
            ['Target size', '33,000,000 bp', 'Run QC outcome:', 'Passed'],
            ['Total Variants', str(total_variants), '', ''],
            ['SNVs', str(len(data.get('snp', pd.DataFrame()))), 
             'INDELs', str(len(data.get('del', pd.DataFrame())) + len(data.get('ins', pd.DataFrame())))]
        ]
        
        assay_table = Table(assay_data, colWidths=[1.5*inch, 1.5*inch, 1.5*inch, 1.5*inch])
        assay_table.setStyle(TableStyle([
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, -1), 'Helvetica'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('GRID', (0, 0), (-1, -1), 1, colors.grey),
        ]))
        
        elements.append(assay_table)
        elements.append(Spacer(1, 20))
        
        # Add variant tables
        elements.append(Paragraph("SNVs identified", self.styles['SectionHeader']))
        if 'snp' in data:
            elements.append(self.create_variant_table(data['snp'], 'SNVs'))
        elements.append(Spacer(1, 20))
        
        elements.append(Paragraph("INDELs identified", self.styles['SectionHeader']))
        indels_df = pd.concat([data.get('del', pd.DataFrame()), 
                             data.get('ins', pd.DataFrame())], 
                            ignore_index=True)
        if not indels_df.empty:
            elements.append(self.create_variant_table(indels_df, 'INDELs'))
        
        # Build the PDF
        doc.build(elements)
        self.logger.info(f"Generated clinical report for {sample_name}")

def main():
    try:
        # Setup paths
        root = Path(__file__).resolve().parents[1]
        results_dir = root / "data" / "results"
        
        # Ensure results directory exists
        if not results_dir.exists():
            raise FileNotFoundError(f"Results directory not found: {results_dir}")
            
        output_dir = results_dir / "reports"
        output_dir.mkdir(exist_ok=True, parents=True)
        
        # Verify processed data directory exists
        processed_dir = results_dir / "processed_data"
        if not processed_dir.exists():
            raise FileNotFoundError(f"Processed data directory not found: {processed_dir}")
            
        logging.info(f"Output directory created at: {output_dir}")
    
    except Exception as e:
        logging.error(f"Error setting up directories: {str(e)}")
        raise
    
    # Initialize report generator
    generator = ClinicalReportGenerator()
    
    try:
        # Process each sample
        samples = ["MMCID-26B", "MMCID-30B"]  # Add other samples as needed
        
        for sample_name in samples:
            # Load processed data
            data = generator.load_processed_data(results_dir, sample_name)
            if not data:
                logging.warning(f"No processed data found for {sample_name}")
                continue
                
            # Generate report
            generator.generate_report(
                sample_name=sample_name,
                data=data,
                output_dir=output_dir
            )
            
    except Exception as e:
        logging.error(f"Error generating reports: {str(e)}")
        raise

if __name__ == "__main__":
    main()