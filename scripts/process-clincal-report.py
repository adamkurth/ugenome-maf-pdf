#!/usr/bin/env python3
import pandas as pd
import json
import shutil
from pathlib import Path
import logging
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib.utils import ImageReader
from reportlab.platypus import Frame, Image, SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
from PyPDF2 import PdfWriter, PdfReader
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
import io

class ClinicalReportGenerator:
    def __init__(self):
        logging.basicConfig(level=logging.INFO,
                          format='%(asctime)s - %(levelname)s - %(message)s')
        self.logger = logging.getLogger(__name__)
        self.styles = getSampleStyleSheet()
        self.gene_visualizer = GeneVisualizer()
        self.setup_styles()
        self.setup_table_styles()

    def setup_styles(self):
        """Set up custom styles for the report"""
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
        
        self.styles.add(ParagraphStyle(
            name='VisualizationCaption',
            parent=self.styles['Normal'],
            fontSize=10,
            textColor=colors.black,
            alignment=1,  # Center alignment
            spaceAfter=20,
            spaceBefore=10
        ))

        self.styles.add(ParagraphStyle(
            name='VisualizationHeader',
            parent=self.styles['Heading3'],
            fontSize=11,
            spaceBefore=20,
            spaceAfter=10,
            textColor=colors.black
        ))
   
    def setup_table_styles(self):
        """Setup modern table styles to be used across all sections"""
        self.modern_table_style = TableStyle([
            # Header row styling
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#2C3E50')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('ALIGN', (0, 0), (-1, 0), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 10),
            
            # Data rows styling
            ('BACKGROUND', (0, 1), (-1, -1), colors.white),
            ('TEXTCOLOR', (0, 1), (0, -1), colors.HexColor('#2C3E50')),
            ('TEXTCOLOR', (1, 1), (-1, -1), colors.HexColor('#34495E')),
            ('ALIGN', (0, 1), (0, -1), 'LEFT'),
            ('ALIGN', (1, 1), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
            ('FONTSIZE', (0, 1), (-1, -1), 9),
            
            # Grid styling
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#BDC3C7')),
            ('LINEBELOW', (0, 0), (-1, 0), 1, colors.HexColor('#2C3E50')),
            
            # Cell padding
            ('TOPPADDING', (0, 0), (-1, -1), 6),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
            ('LEFTPADDING', (0, 0), (-1, -1), 8),
            ('RIGHTPADDING', (0, 0), (-1, -1), 8),
        ])

        # Function to add alternating row colors
        self.add_row_colors = lambda table_style, num_rows: [
            table_style.add(('BACKGROUND', (0, i), (-1, i), colors.HexColor('#F8F9F9')))
            for i in range(1, num_rows, 2)
        ]

    def analyze_processed_data(self, processed_dir: Path) -> dict:
        """Analyze processed data with enhanced gene analysis"""
        analysis = {
            'total_variants': 0,
            'variant_types': {},
            'gene_frequency': {},
            'chromosome_distribution': {},
            'high_impact_genes': set()
        }
        
        try:
            for complete_file in processed_dir.glob('*_complete.csv'):
                df = pd.read_csv(complete_file)
                
                analysis['total_variants'] += len(df)
                
                variant_counts = df['Variant_Type'].value_counts()
                for variant_type, count in variant_counts.items():
                    analysis['variant_types'][variant_type] = analysis['variant_types'].get(variant_type, 0) + count
                
                gene_counts = df['Hugo_Symbol'].value_counts()
                for gene, count in gene_counts.items():
                    analysis['gene_frequency'][gene] = analysis['gene_frequency'].get(gene, 0) + count
                    
                    if 'IMPACT' in df.columns:
                        high_impact_mutations = df[
                            (df['Hugo_Symbol'] == gene) & 
                            (df['IMPACT'].isin(['HIGH', 'MODERATE']))
                        ]
                        if not high_impact_mutations.empty:
                            analysis['high_impact_genes'].add(gene)
                
                chrom_counts = df['Chromosome'].value_counts()
                for chrom, count in chrom_counts.items():
                    analysis['chromosome_distribution'][chrom] = analysis['chromosome_distribution'].get(chrom, 0) + count
                    
            analysis['high_impact_genes'] = list(analysis['high_impact_genes'])
            
            return analysis
            
        except Exception as e:
            self.logger.error(f"Error analyzing processed data: {str(e)}")
            raise

    def load_processed_data(self, results_dir: Path, sample_name: str) -> dict:
        """Load processed data and prepare it for the report"""
        processed_dir = results_dir / 'processed_data'
        data = {}
        
        try:
            complete_file = processed_dir / f'{sample_name}_complete.csv'
            if complete_file.exists():
                self.logger.info(f"Loading complete data from {complete_file}")
                data['complete'] = pd.read_csv(complete_file)
                
                report_data_dir = results_dir / 'reports' / 'data'
                report_data_dir.mkdir(exist_ok=True)
                shutil.copy2(complete_file, report_data_dir / f'{sample_name}_complete.csv')
            
            for variant_type in ['snp', 'del', 'ins']:
                variant_file = processed_dir / f'{sample_name}_{variant_type}.csv'
                if variant_file.exists():
                    df = pd.read_csv(variant_file)
                    if 'IMPACT' in df.columns:
                        df = df.sort_values('IMPACT', ascending=False).head(5)
                    else:
                        df = df.head(5)
                    data[variant_type] = df
                    
                    shutil.copy2(variant_file, report_data_dir / f'{sample_name}_{variant_type}.csv')
                    
            return data
            
        except Exception as e:
            self.logger.error(f"Error loading data for {sample_name}: {str(e)}")
            raise

    def create_modern_table(self, data: list, col_widths: list = None) -> Table:
        """Create a table with modern styling, avoiding negative indexing."""
        if not data or not data[0]:
            raise ValueError("No valid data provided for the table.")

        row_count = len(data)
        col_count = len(data[0])

        if col_widths is None:
            col_widths = [2 * inch] * col_count

        table = Table(data, colWidths=col_widths)

        # Base style without negative indexing:
        style_commands = [
            # Header row styling
            ('BACKGROUND', (0, 0), (col_count - 1, 0), colors.HexColor('#2C3E50')),
            ('TEXTCOLOR', (0, 0), (col_count - 1, 0), colors.white),
            ('ALIGN', (0, 0), (col_count - 1, 0), 'CENTER'),
            ('FONTNAME', (0, 0), (col_count - 1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (col_count - 1, 0), 10),

            # If there are data rows beyond the header:
            # Apply these only if row_count > 1
            # (If there's just one row, it's only the header.)
        ]

        if row_count > 1:
            style_commands.extend([
                ('BACKGROUND', (0, 1), (col_count - 1, row_count - 1), colors.white),
                ('TEXTCOLOR', (0, 1), (0, row_count - 1), colors.HexColor('#2C3E50')),
                ('TEXTCOLOR', (1, 1), (col_count - 1, row_count - 1), colors.HexColor('#34495E')),
                ('ALIGN', (0, 1), (0, row_count - 1), 'LEFT'),
                ('ALIGN', (1, 1), (col_count - 1, row_count - 1), 'CENTER'),
                ('FONTNAME', (0, 1), (col_count - 1, row_count - 1), 'Helvetica'),
                ('FONTSIZE', (0, 1), (col_count - 1, row_count - 1), 9),
            ])

        # Grid and lines apply to entire table area
        style_commands.extend([
            ('GRID', (0, 0), (col_count - 1, row_count - 1), 0.5, colors.HexColor('#BDC3C7')),
            ('LINEBELOW', (0, 0), (col_count - 1, 0), 1, colors.HexColor('#2C3E50')),

            # Cell padding for all cells
            ('TOPPADDING', (0, 0), (col_count - 1, row_count - 1), 6),
            ('BOTTOMPADDING', (0, 0), (col_count - 1, row_count - 1), 6),
            ('LEFTPADDING', (0, 0), (col_count - 1, row_count - 1), 8),
            ('RIGHTPADDING', (0, 0), (col_count - 1, row_count - 1), 8),
        ])

        # Add alternating row background color starting from the second row (row index 1)
        # Only if we have multiple rows.
        if row_count > 2:
            for i in range(1, row_count, 2):
                style_commands.append(('BACKGROUND', (0, i), (col_count - 1, i), colors.HexColor('#F8F9F9')))

        style = TableStyle(style_commands)
        table.setStyle(style)
        return table

    def create_header_section(self) -> list:
        """Create the report header section with modern styling as text."""
        elements = []
        
        # Create a styled section heading
        elements.append(Paragraph("Report Header", self.styles['SectionHeader']))
        
        report_id = f"MAF-{datetime.now().strftime('%Y%m%d')}"
        report_date = datetime.now().strftime('%Y-%m-%d')
        
        # Create paragraphs for each header field:
        elements.append(Paragraph(f"<b>Report N°:</b> {report_id}", self.styles['Normal']))
        elements.append(Spacer(1, 4))
        elements.append(Paragraph(f"<b>Date:</b> {report_date}", self.styles['Normal']))
        elements.append(Spacer(1, 20))
        
        return elements

    def create_patient_section(self, sample_name: str) -> list:
        """Create the patient information section as styled text, maintaining a modern look."""
        elements = []
        
        # Create a styled section heading
        elements.append(Paragraph("Patient Information", self.styles['SectionHeader']))
        
        # Patient fields as styled paragraphs
        elements.append(Paragraph(f"<b>PATIENT:</b> {sample_name}", self.styles['Normal']))
        elements.append(Spacer(1, 4))
        elements.append(Paragraph(f"<b>ID N°:</b> XXXXX", self.styles['Normal']))
        elements.append(Spacer(1, 4))
        elements.append(Paragraph(f"<b>Clinical Diagnosis:</b> mCRC", self.styles['Normal']))
        elements.append(Spacer(1, 4))
        elements.append(Paragraph(f"<b>Clinical Trial:</b> _____________", self.styles['Normal']))
        elements.append(Spacer(1, 20))
        
        return elements

    def create_tmb_section(self, data: dict) -> list:
        """Create TMB summary section with modern styling"""
        elements = []
        
        if 'complete' in data:
            tmb = self.calculate_tmb(data['complete'])
            
            tmb_data = [
                ['Metric', 'Value', 'Details'],
                ['Tumor Mutation Burden (TMB)', f"{tmb:.2f}", 'mutations/Mb'],
                ['Total Variants', f"{len(data['complete']):,}", 'All variant types'],
                ['SNVs', f"{len(data['complete'][data['complete']['Variant_Type'] == 'SNP']):,}", 'Single nucleotide variants'],
                ['Insertions', f"{len(data['complete'][data['complete']['Variant_Type'] == 'INS']):,}", 'Insertion events'],
                ['Deletions', f"{len(data['complete'][data['complete']['Variant_Type'] == 'DEL']):,}", 'Deletion events']
            ]
            
            elements.append(Paragraph("Mutation Burden Analysis", self.styles['SectionHeader']))
            elements.append(self.create_modern_table(tmb_data, [2.5*inch, 1.5*inch, 2*inch]))
            elements.append(Spacer(1, 15))
        
        return elements

    def create_variant_table(self, df: pd.DataFrame, variant_type: str) -> Table:
        """Create an enhanced table for variants with modern styling"""
        header_row = ['Gene', 'Chr', 'Position', 'Ref', 'Alt', 'Protein Change']
        
        table_data = [header_row]
        for _, row in df.head(5).iterrows():
            formatted_row = [
                str(row['Hugo_Symbol']),
                str(row['Chromosome']),
                f"{row['Start_Position']:,}",
                str(row['Reference_Allele']),
                str(row['Tumor_Seq_Allele2']),
                str(row['HGVSp_Short']) if pd.notna(row['HGVSp_Short']) else '-'
            ]
            table_data.append(formatted_row)
        
        col_widths = [1.5*inch, 0.6*inch, 1.2*inch, 0.8*inch, 0.8*inch, 1.5*inch]
        return self.create_modern_table(table_data, col_widths)

    def calculate_tmb(self, df: pd.DataFrame) -> float:
        """
        Calculate Tumor Mutation Burden (mutations per megabase)
        """
        try:
            # Filter for relevant variants (typically SNVs and INDELs)
            relevant_variants = df[
                df['Variant_Type'].isin(['SNP', 'INS', 'DEL'])
            ]
            
            # Calculate TMB (assuming ~30 megabases for typical exome)
            exome_size_mb = 30  # Size of typical exome in megabases
            tmb = len(relevant_variants) / exome_size_mb
            
            return round(tmb, 2)
        
        except Exception as e:
            self.logger.error(f"Error calculating TMB: {str(e)}")
            return 0.0

    def create_summary_section(self, analysis_data: dict) -> list:
        """Create summary section with modern styling"""
        elements = []
        
        # Add summary header with enhanced styling
        elements.append(Paragraph("Analysis Summary", self.styles['SectionHeader']))
        
        # Calculate percentage distributions for variant types
        total_variants = sum(analysis_data['variant_types'].values())
        variant_percentages = {k: (v/total_variants)*100 
                            for k, v in analysis_data['variant_types'].items()}
        
        # Create enhanced summary table with more metrics and better formatting
        summary_data = [
            ['Metric', 'Value', 'Additional Info'],
            ['Total Variants', f"{analysis_data['total_variants']:,}", ''],
            ['Primary Variant Type', 
            max(analysis_data['variant_types'].items(), key=lambda x: x[1])[0],
            f"{max(variant_percentages.values()):.1f}% of total"],
            ['Most Mutated Gene',
            next(iter(analysis_data['gene_frequency'])),
            f"{max(analysis_data['gene_frequency'].values())} mutations"],
            ['Most Affected Chr',
            next(iter(analysis_data['chromosome_distribution'])),
            f"{max(analysis_data['chromosome_distribution'].values())} variants"],
            ['High Impact Genes',
            f"{len(analysis_data.get('high_impact_genes', []))}", 
            'Genes with significant mutations']
        ]
        
        # Create table with modern styling
        table = Table(summary_data, colWidths=[2*inch, 2*inch, 2.5*inch])
        table.setStyle(TableStyle([
            # Header row styling
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#2C3E50')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('ALIGN', (0, 0), (-1, 0), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 10),
            
            # Data rows styling
            ('BACKGROUND', (0, 1), (-1, -1), colors.white),
            ('TEXTCOLOR', (0, 1), (0, -1), colors.HexColor('#2C3E50')),  # Metric names
            ('TEXTCOLOR', (1, 1), (2, -1), colors.HexColor('#34495E')),  # Values
            ('ALIGN', (0, 0), (0, -1), 'LEFT'),   # Left align first column
            ('ALIGN', (1, 0), (-1, -1), 'CENTER'), # Center align other columns
            ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
            ('FONTSIZE', (0, 1), (-1, -1), 9),
            
            # Alternating row colors
            ('BACKGROUND', (0, 1), (-1, 1), colors.HexColor('#F8F9F9')),
            ('BACKGROUND', (0, 3), (-1, 3), colors.HexColor('#F8F9F9')),
            ('BACKGROUND', (0, 5), (-1, 5), colors.HexColor('#F8F9F9')),
            
            # Grid styling
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#BDC3C7')),
            ('LINEBELOW', (0, 0), (-1, 0), 1, colors.HexColor('#2C3E50')),
            
            # Cell padding
            ('TOPPADDING', (0, 0), (-1, -1), 6),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
            ('LEFTPADDING', (0, 0), (-1, -1), 8),
            ('RIGHTPADDING', (0, 0), (-1, -1), 8),
        ]))
        
        elements.append(table)
        elements.append(Spacer(1, 20))
        
        return elements 
    
    def create_gene_section(self, data: dict) -> list:
        """Create section for gene-specific analysis"""
        elements = []
        
        if 'gene_visualizations' in data:
            elements.append(Paragraph("Gene Analysis", self.styles['SectionHeader']))
            
            # Add gene-specific visualizations
            for gene, viz_path in data['gene_visualizations'].items():
                elements.append(Paragraph(f"Mutation Profile - {gene}", self.styles['TableHeader']))
                elements.append(Image(viz_path, width=6*inch, height=3*inch))
                elements.append(Spacer(1, 10))
            
            # Add impact distribution visualization
            if 'impact_visualization' in data:
                elements.append(Paragraph("Impact Distribution", self.styles['TableHeader']))
                elements.append(Image(data['impact_visualization'], width=6*inch, height=3*inch))
                
        return elements

    def create_visualization_section(self, sample_name: str, data: dict) -> list:
        """Create visualization section with embedded visualizations"""
        elements = []
        
        elements.append(PageBreak())
        elements.append(Paragraph("Visualizations", self.styles['SectionHeader']))
        
        viz_dir = Path("data") / "results" / "visualizations"
        
        if not viz_dir.exists():
            self.logger.warning(f"Visualization directory not found at: {viz_dir}")
            return elements
        
        # Add mutation signature plot if available
        mutation_signature = viz_dir / f"{sample_name}_mutation_signature.png"
        if mutation_signature.exists():
            elements.append(Paragraph("Mutation Signature", self.styles['TableHeader']))
            elements.append(Image(str(mutation_signature), width=6*inch, height=3*inch))
            elements.append(Spacer(1, 10))
            
        # Add variant distribution plot
        variant_dist = viz_dir / f"{sample_name}_variant_distribution.png"
        if variant_dist.exists():
            elements.append(Paragraph("Variant Distribution", self.styles['TableHeader']))
            elements.append(Image(str(variant_dist), width=6*inch, height=3*inch))
            elements.append(Spacer(1, 10))
            
        # Add summary dashboard
        dashboard = viz_dir / f"{sample_name}_summary_dashboard.png"
        if dashboard.exists():
            elements.append(Paragraph("Summary Dashboard", self.styles['TableHeader']))
            elements.append(Image(str(dashboard), width=7*inch, height=5*inch))
            elements.append(Spacer(1, 10))

        # Add comparison plots if they exist
        for viz_path in viz_dir.glob(f"comparison*{sample_name}*.png"):
            comparison_name = viz_path.stem.replace('comparison_', '').replace('_', ' vs ')
            elements.append(Paragraph(f"Sample Comparison: {comparison_name}", 
                                    self.styles['TableHeader']))
            elements.append(Image(str(viz_path), width=7*inch, height=4*inch))
            elements.append(Spacer(1, 10))

        return elements
                        
    def create_visualization_section(self, sample_name: str, data: dict) -> list:
        """Create visualization section with improved layout and compact visuals"""
        elements = []
        
        elements.append(PageBreak())
        elements.append(Paragraph("Visualizations", self.styles['SectionHeader']))
        
        viz_dir = Path("data") / "results" / "visualizations"
        
        if not viz_dir.exists():
            self.logger.warning(f"Visualization directory not found at: {viz_dir}")
            return elements

            
        def add_compact_visualization(image_path, title):
            """Helper function to add visualizations with compact scaling"""
            if image_path.exists():
                elements.append(Paragraph(title, self.styles['VisualizationHeader']))
                
                img = Image(str(image_path))
                aspect_ratio = img.imageHeight / img.imageWidth
                
                # Set standard width slightly smaller for compactness
                standard_width = 5 * inch  # Reduced from 6 inches
                target_height = standard_width * aspect_ratio
                
                # Scale height if too large
                max_height = 3.5 * inch  # Reduced max height
                if target_height > max_height:
                    target_height = max_height
                    standard_width = target_height / aspect_ratio
                
                elements.append(Image(str(image_path), width=standard_width, height=target_height))
                elements.append(Spacer(1, 15))  # Reduced spacing between visualizations
        
        # 1. Summary Dashboard
        dashboard = viz_dir / f"{sample_name}_summary_dashboard.png"
        if dashboard.exists():
            add_compact_visualization(dashboard, "Summary Analysis")
        
        # 2. Variant Distribution
        variant_dist = viz_dir / f"{sample_name}_variant_distribution.png"
        if variant_dist.exists():
            add_compact_visualization(variant_dist, "Variant Type Distribution")
        
        # 3. Mutation Signature
        mutation_signature = viz_dir / f"{sample_name}_mutation_signature.png"
        if mutation_signature.exists():
            add_compact_visualization(mutation_signature, "Mutation Signature Profile")
        
        # 4. Sample Comparisons
        for viz_path in viz_dir.glob(f"comparison*{sample_name}*.png"):
            comparison_name = viz_path.stem.replace('comparison_', '').replace('_', ' vs ')
            add_compact_visualization(viz_path, f"Sample Comparison: {comparison_name}")

        return elements


    def add_visualization_with_caption(self, image_path: Path, caption: str) -> list:
        """Add visualization with proper caption and layout"""
        elements = []
        
        # Add caption
        elements.append(Paragraph(caption, self.styles['VisualizationCaption']))
        
        # Add image with proper scaling
        img = Image(str(image_path))
        aspect_ratio = img.imageHeight / img.imageWidth
        target_width = 6.5 * inch  # Max width for letter size with margins
        target_height = target_width * aspect_ratio
        
        # Scale down if height is too large
        if target_height > 8 * inch:
            target_height = 8 * inch
            target_width = target_height / aspect_ratio
            
        elements.append(Image(str(image_path), width=target_width, height=target_height))
        elements.append(Spacer(1, 10))
        
        return elements

    def generate_report(self, sample_name: str, data: dict, analysis_data: dict, output_dir: Path):
        """Generate complete clinical report with correct section ordering"""
        doc = SimpleDocTemplate(
            str(output_dir / f"{sample_name}_clinical_report.pdf"),
            pagesize=letter,
            rightMargin=72,
            leftMargin=72,
            topMargin=72,
            bottomMargin=72
        )
        
        elements = []
        
        # 1. Basic header and patient info
        elements.extend(self.create_header_section())
        elements.extend(self.create_patient_section(sample_name))
        
        # 2. Analysis Summary
        elements.extend(self.create_summary_section(analysis_data))
        
        # Add TMB section before variant tables
        elements.extend(self.create_tmb_section(data))

        # Variant Analysis section with enhanced tables
        elements.append(Paragraph("Variant Analysis", self.styles['SectionHeader']))
        
        if 'snp' in data:
            elements.append(Paragraph("Single Nucleotide Variants (SNV) (Top 5)", self.styles['TableHeader']))
            elements.append(self.create_variant_table(data['snp'], 'SNVs'))
            elements.append(Spacer(1, 15))
        
        indels_df = pd.concat([data.get('del', pd.DataFrame()), 
                            data.get('ins', pd.DataFrame())], 
                            ignore_index=True).head(5)
        if not indels_df.empty:
            elements.append(Paragraph("Insertions and Deletions (INDELS) (Top 5)", self.styles['TableHeader']))
            elements.append(self.create_variant_table(indels_df, 'INDELs'))
            elements.append(Spacer(1, 15))

        # 4. Visualizations Section
        elements.append(PageBreak())
        elements.append(Paragraph("Visualizations", self.styles['SectionHeader']))
        
        viz_dir = Path("data") / "results" / "visualizations"
        if not viz_dir.exists():
            self.logger.warning(f"Visualization directory not found at: {viz_dir}")
        else:
            # Summary Dashboard first
            dashboard = viz_dir / f"{sample_name}_summary_dashboard.png"
            if dashboard.exists():
                elements.append(Paragraph("Summary Dashboard", self.styles['VisualizationHeader']))
                elements.append(Image(str(dashboard), width=6*inch, height=4*inch))
                elements.append(Spacer(1, 15))
            
            # Variant Distribution
            variant_dist = viz_dir / f"{sample_name}_variant_distribution.png"
            if variant_dist.exists():
                elements.append(Paragraph("Distribution of Variant Types", self.styles['VisualizationHeader']))
                elements.append(Image(str(variant_dist), width=6*inch, height=3*inch))
                elements.append(Spacer(1, 15))
            
            # Mutation Signature
            mutation_signature = viz_dir / f"{sample_name}_mutation_signature.png"
            if mutation_signature.exists():
                elements.append(Paragraph("Mutation Signature Analysis", self.styles['VisualizationHeader']))
                elements.append(Image(str(mutation_signature), width=6*inch, height=3*inch))
                elements.append(Spacer(1, 15))
            
            # Sample Comparisons
            for viz_path in viz_dir.glob(f"comparison*{sample_name}*.png"):
                comparison_name = viz_path.stem.replace('comparison_', '').replace('_', ' vs ')
                elements.append(Paragraph(f"Sample Comparison: {comparison_name}", self.styles['VisualizationHeader']))
                elements.append(Image(str(viz_path), width=6*inch, height=3*inch))
                elements.append(Spacer(1, 15))
        
        # 5. Data location note
        elements.append(Spacer(1, 30))
        elements.append(Paragraph(
            f"Complete data files available in: {output_dir}/data/",
            self.styles['Normal']
        ))
        
        doc.build(elements)
        self.logger.info(f"Generated clinical report for {sample_name}")
        
class GeneVisualizer:
    """Class for generating gene-specific visualizations"""
    def __init__(self):
        self.logger = logging.getLogger(__name__)

    def plot_gene_mutation_profile(self, df: pd.DataFrame, gene: str, output_path: str):
        """Create detailed mutation profile for a specific gene"""
        try:
            # Convert output path to PNG
            output_path = str(output_path).replace('.pdf', '.png')
            
            gene_data = df[df['Hugo_Symbol'] == gene].copy()
            if gene_data.empty:
                return None

            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

            # Plot 1: Mutation types for this gene
            mutation_counts = gene_data['Variant_Classification'].value_counts()
            sns.barplot(x=mutation_counts.values, y=mutation_counts.index, ax=ax1)
            ax1.set_title(f'Mutation Types in {gene}')
            ax1.set_xlabel('Count')

            # Plot 2: Mutation positions along the gene
            if 'Start_Position' in gene_data.columns:
                sns.histplot(data=gene_data, x='Start_Position', ax=ax2)
                ax2.set_title(f'Mutation Positions in {gene}')
                ax2.set_xlabel('Genomic Position')

            plt.tight_layout()
            # Save as PNG with high DPI for better quality
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            return output_path

        except Exception as e:
            self.logger.error(f"Error creating gene profile for {gene}: {str(e)}")
            return None

    def plot_gene_impact_distribution(self, df: pd.DataFrame, genes: list, output_path: str):
        """Plot impact distribution for specified genes"""
        try:
            # Convert output path to PNG
            output_path = str(output_path).replace('.pdf', '.png')
            
            gene_data = df[df['Hugo_Symbol'].isin(genes)].copy()
            if gene_data.empty:
                return None

            plt.figure(figsize=(10, 6))
            impact_data = pd.crosstab(gene_data['Hugo_Symbol'], gene_data['IMPACT'])
            
            # Create a stacked bar plot
            impact_data.plot(kind='bar', stacked=True)
            plt.title('Mutation Impact Distribution by Gene')
            plt.xlabel('Gene')
            plt.ylabel('Count')
            plt.xticks(rotation=45)
            plt.legend(title='Impact', bbox_to_anchor=(1.05, 1))
            plt.tight_layout()
            
            # Save as PNG with high DPI for better quality
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            return output_path

        except Exception as e:
            self.logger.error(f"Error creating impact distribution plot: {str(e)}")
            return None

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
        # samples = ["TCMK1-14B", "TMCK1-23B"]
        
        for sample_name in samples:
            # Load processed data
            data = generator.load_processed_data(results_dir, sample_name)
            analysis_data = generator.analyze_processed_data(results_dir / 'processed_data')
            
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