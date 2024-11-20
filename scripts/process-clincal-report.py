#!/usr/bin/env python3
import pandas as pd
import json
import shutil
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

    def analyze_processed_data(self, processed_dir: Path) -> dict:
        """Analyze all processed data files"""
        analysis = {
            'total_variants': 0,
            'variant_types': {},
            'gene_frequency': {},
            'chromosome_distribution': {}
        }
        
        try:
            # Read all complete.csv files
            for complete_file in processed_dir.glob('*_complete.csv'):
                df = pd.read_csv(complete_file)
                
                # Update total variants
                analysis['total_variants'] += len(df)
                
                # Update variant types
                variant_counts = df['Variant_Type'].value_counts()
                for variant_type, count in variant_counts.items():
                    analysis['variant_types'][variant_type] = analysis['variant_types'].get(variant_type, 0) + count
                
                # Update gene frequency
                gene_counts = df['Hugo_Symbol'].value_counts()
                for gene, count in gene_counts.items():
                    analysis['gene_frequency'][gene] = analysis['gene_frequency'].get(gene, 0) + count
                
                # Update chromosome distribution
                chrom_counts = df['Chromosome'].value_counts()
                for chrom, count in chrom_counts.items():
                    analysis['chromosome_distribution'][chrom] = analysis['chromosome_distribution'].get(chrom, 0) + count
                    
            # Sort dictionaries by value
            analysis['gene_frequency'] = dict(sorted(analysis['gene_frequency'].items(), 
                                                   key=lambda x: x[1], reverse=True))
            analysis['chromosome_distribution'] = dict(sorted(analysis['chromosome_distribution'].items(), 
                                                           key=lambda x: x[1], reverse=True))
            
            return analysis
            
        except Exception as e:
            self.logger.error(f"Error analyzing processed data: {str(e)}")
            raise

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
                
                # Create a copy of the complete file in the reports directory
                report_data_dir = results_dir / 'reports' / 'data'
                report_data_dir.mkdir(exist_ok=True)
                shutil.copy2(complete_file, report_data_dir / f'{sample_name}_complete.csv')
            
            # Load variant-specific data and get top 5 for each type
            for variant_type in ['snp', 'del', 'ins']:
                variant_file = processed_dir / f'{sample_name}_{variant_type}.csv'
                if variant_file.exists():
                    df = pd.read_csv(variant_file)
                    # Keep only top 5 variants by impact or another relevant metric
                    if 'IMPACT' in df.columns:
                        df = df.sort_values('IMPACT', ascending=False).head(5)
                    else:
                        df = df.head(5)
                    data[variant_type] = df
                    
                    # Copy variant files to reports directory
                    shutil.copy2(variant_file, report_data_dir / f'{sample_name}_{variant_type}.csv')
                    
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
        """Create a table for top 5 variants"""
        columns = ['Hugo_Symbol', 'Chromosome', 'Start_Position', 
                  'Reference_Allele', 'Tumor_Seq_Allele2', 'HGVSp_Short']
        
        # Create table data with only first 5 rows
        table_data = [columns]
        for _, row in df.head(5).iterrows():
            table_data.append([str(row[col]) for col in columns])
            
        table = Table(table_data, repeatRows=1)
        table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#F5F5F5')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.black),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 8),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('GRID', (0, 0), (-1, -1), 1, colors.black)
        ]))
        
        return table

    def create_summary_section(self, analysis_data: dict) -> list:
        """Create summary section with overall analysis"""
        elements = []
        
        # Add summary header
        elements.append(Paragraph("Analysis Summary", self.styles['SectionHeader']))
        
        # Create summary table
        summary_data = [
            ['Total Variants', str(analysis_data['total_variants'])],
            ['Most Common Variant Type', max(analysis_data['variant_types'].items(), key=lambda x: x[1])[0]],
            ['Most Mutated Gene', next(iter(analysis_data['gene_frequency']))],
            ['Most Affected Chromosome', next(iter(analysis_data['chromosome_distribution']))]
        ]
        
        summary_table = Table(summary_data, colWidths=[2*inch, 3*inch])
        summary_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#F5F5F5')),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('FONTNAME', (0, 0), (-1, -1), 'Helvetica'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT')
        ]))
        
        elements.append(summary_table)
        elements.append(Spacer(1, 20))
        
        return elements

    def generate_report(self, sample_name: str, data: dict, analysis_data: dict, output_dir: Path):
        """Generate complete clinical report"""
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
        
        # Add overall analysis summary
        elements.extend(self.create_summary_section(analysis_data))
        
        # Add variant tables (top 5 only)
        elements.append(Paragraph("Top 5 SNVs", self.styles['SectionHeader']))
        if 'snp' in data:
            elements.append(self.create_variant_table(data['snp'], 'SNVs'))
        elements.append(Spacer(1, 20))
        
        elements.append(Paragraph("Top 5 INDELs", self.styles['SectionHeader']))
        indels_df = pd.concat([data.get('del', pd.DataFrame()), 
                             data.get('ins', pd.DataFrame())], 
                            ignore_index=True).head(5)
        if not indels_df.empty:
            elements.append(self.create_variant_table(indels_df, 'INDELs'))
        
        # Add data location note
        elements.append(Spacer(1, 30))
        elements.append(Paragraph(
            f"Complete data files available in: {output_dir}/data/",
            self.styles['Normal']
        ))
        
        # Build the PDF
        doc.build(elements)
        self.logger.info(f"Generated clinical report for {sample_name}")

def main():
    try:
        # Setup paths
        root = Path(__file__).resolve().parents[1]
        results_dir = root / "data" / "results"
        
        if not results_dir.exists():
            raise FileNotFoundError(f"Results directory not found: {results_dir}")
            
        output_dir = results_dir / "reports"
        output_dir.mkdir(exist_ok=True, parents=True)
        
        processed_dir = results_dir / "processed_data"
        if not processed_dir.exists():
            raise FileNotFoundError(f"Processed data directory not found: {processed_dir}")
            
        # Initialize report generator
        generator = ClinicalReportGenerator()
        
        # Analyze all processed data
        analysis_data = generator.analyze_processed_data(processed_dir)
        
        # Process each sample
        samples = ["MMCID-26B", "MMCID-30B"]
        
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
                analysis_data=analysis_data,
                output_dir=output_dir
            )
            
        logging.info("Report generation complete. Data files copied to reports/data/")
            
    except Exception as e:
        logging.error(f"Error in report generation: {str(e)}")
        raise

if __name__ == "__main__":
    main()