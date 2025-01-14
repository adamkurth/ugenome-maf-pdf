#!/usr/bin/env python3
import pandas as pd
import numpy as np
from pathlib import Path
import logging
import json
from typing import Dict, List, Tuple, Optional
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# from maf_processor import MAFProcessor
# from mutation_signature import MutationSignatureAnalyzer
# from maf_visualizer import MAFVisualizer

class MAFProcessor:
    """Base class for MAF file processing and validation"""
    def __init__(self):
        self._setup_logging()
        self._setup_column_definitions()
        self._setup_valid_classifications()

    def _setup_logging(self):
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)


    def _setup_column_definitions(self):
        """Define required MAF columns"""
        self.required_columns = [
            'Hugo_Symbol', 'Chromosome', 'Start_Position', 
            'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele2',
            'Tumor_Sample_Barcode'
        ]
        
        self.variant_info_columns = [
            'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'End_Position', 'Strand',
            'Variant_Classification', 'dbSNP_RS', 'dbSNP_Val_Status'
        ]
        
        self.consequence_columns = [
            'HGVSc', 'HGVSp', 'HGVSp_Short', 'Transcript_ID',
            'Exon_Number', 'all_effects', 'Consequence', 'IMPACT',
            'SIFT', 'PolyPhen', 'CLIN_SIG'
        ]

        self.coverage_columns = [
            't_depth', 't_ref_count', 't_alt_count',
            'n_depth', 'n_ref_count', 'n_alt_count'
        ]
    
    def _setup_valid_classifications(self):
        """Define valid mutation classifications"""
        self.valid_classifications = {
            'Missense_Mutation', 'Nonsense_Mutation', 'Silent', 'Splice_Site',
            'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins',
            'Intron', '3\'UTR', '5\'UTR', 'IGR', '3\'Flank', '5\'Flank', 'RNA',
            'Splice_Region', 'Targeted_Region'
        }

    def read_maf(self, filepath: str) -> Optional[pd.DataFrame]:
        """Read and validate MAF file"""
        try:
            df = pd.read_csv(filepath, sep='\t', comment='#', low_memory=False)
            
            # Basic validation
            required_cols = ['Hugo_Symbol', 'Chromosome', 'Start_Position', 'Variant_Type']
            missing = [col for col in required_cols if col not in df.columns]
            if missing:
                self.logger.warning(f"Missing required columns in {filepath}: {missing}")
                return None

            # Data cleaning
            df['Chromosome'] = df['Chromosome'].astype(str)
            df = df.fillna({
                'Variant_Classification': 'Unknown',
                'Reference_Allele': '-',
                'Tumor_Seq_Allele1': '-',
                'Tumor_Seq_Allele2': '-'
            })
            
            return df
            
        except Exception as e:
            self.logger.error(f"Error reading MAF file {filepath}: {str(e)}")
            return None
    
    def validate_maf_data(self, df: pd.DataFrame) -> Tuple[bool, List[str]]:
        """Validate MAF data quality"""
        issues = []
        
        try:
            # Validate chromosome format
            invalid_chroms = df[~df['Chromosome'].str.match(r'^(chr)?\d+$|^(chr)?[XYM]$|^MT$')]
            if not invalid_chroms.empty:
                issues.append(f"Invalid chromosome formats: {invalid_chroms['Chromosome'].unique()}")
            
            # Validate coordinates
            if (df['Start_Position'] < 0).any():
                issues.append("Found negative genomic coordinates")
            
            if 'End_Position' in df.columns and (df['End_Position'] < df['Start_Position']).any():
                issues.append("Found End_Position < Start_Position")
            
            # Validate variant classifications
            if 'Variant_Classification' in df.columns:
                invalid_class = df[~df['Variant_Classification'].isin(self.valid_classifications)]
                if not invalid_class.empty:
                    issues.append(f"Invalid classifications: {invalid_class['Variant_Classification'].unique()}")
            
            # Check for duplicates
            duplicates = df.groupby(['Chromosome', 'Start_Position', 'Reference_Allele', 
                                   'Tumor_Seq_Allele2']).size().reset_index(name='count')
            duplicates = duplicates[duplicates['count'] > 1]
            if not duplicates.empty:
                issues.append(f"Found {len(duplicates)} duplicate mutations")
            
            return len(issues) == 0, issues
            
        except Exception as e:
            issues.append(f"Validation error: {str(e)}")
            return False, issues

    def calculate_basic_metrics(self, df: pd.DataFrame) -> Dict:
        """Calculate basic metrics from MAF data"""
        try:
            metrics = {
                'total_variants': len(df),
                'variant_types': df['Variant_Type'].value_counts().to_dict(),
                'by_chromosome': df['Chromosome'].value_counts().to_dict(),
                'variant_classification': (df['Variant_Classification'].value_counts().to_dict() 
                                        if 'Variant_Classification' in df.columns else {})
            }
            
            # Add variant allele frequencies if coverage data available
            if all(col in df.columns for col in ['t_depth', 't_alt_count']):
                df['VAF'] = df['t_alt_count'] / df['t_depth']
                metrics['vaf_stats'] = {
                    'mean': df['VAF'].mean(),
                    'median': df['VAF'].median(),
                    'std': df['VAF'].std()
                }
                
            return metrics
            
        except Exception as e:
            self.logger.error(f"Error calculating basic metrics: {str(e)}")
            return {}
            
    def calculate_gene_metrics(self, df: pd.DataFrame) -> Dict:
        """Calculate gene-level metrics"""
        try:
            metrics = {
                'genes_mutated': df['Hugo_Symbol'].nunique(),
                'top_mutated_genes': df['Hugo_Symbol'].value_counts().head(20).to_dict()
            }
            
            if 'Variant_Classification' in df.columns:
                gene_variant_counts = df.groupby('Hugo_Symbol')['Variant_Classification'].value_counts()
                metrics['gene_variant_types'] = gene_variant_counts.to_dict()
                
            return metrics
            
        except Exception as e:
            self.logger.error(f"Error calculating gene metrics: {str(e)}")
            return {}

class MutationSignatureAnalyzer:
    """Analyze mutation signatures from MAF data"""
    def __init__(self):
        self.setup_substitutions()
        self.logger = logging.getLogger(__name__)

    def setup_substitutions(self):
        """Setup mutation substitution types"""
        self.CONTEXTS = ["A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T",
                        "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T",
                        "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T",
                        "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T"]
        
        self.COLORS = {
            'C>A': 'skyblue',
            'C>G': 'black',
            'C>T': 'red',
            'T>A': 'grey',
            'T>C': 'green',
            'T>G': 'pink'
        }

    def calculate_signature(self, df: pd.DataFrame) -> Dict:
        """Calculate mutation signature - matches the pipeline's expected interface"""
        try:
            snp_df = df[df['Variant_Type'] == 'SNP'].copy()
            
            # Process reference and alternate alleles
            ref_alleles = snp_df['Reference_Allele'].fillna('-')
            alt_alleles = snp_df['Tumor_Seq_Allele2'].fillna('-')
            
            # Calculate substitutions
            substitutions = []
            for ref, alt in zip(ref_alleles, alt_alleles):
                if len(ref) == 1 and len(alt) == 1:
                    substitutions.append(f"{ref}>{alt}")
            
            # Calculate percentages
            total = len(substitutions)
            percentages = []
            for context in self.CONTEXTS:
                sub_type = context.split('[')[1].split(']')[0]
                count = sum(1 for s in substitutions if s == sub_type)
                percentages.append((count/total * 100) if total > 0 else 0)
            
            return {
                'contexts': self.CONTEXTS,
                'percentages': percentages,
                'colors': [self.COLORS[ctx.split('[')[1].split(']')[0]] 
                          for ctx in self.CONTEXTS],
                'total_mutations': total
            }
            
        except Exception as e:
            self.logger.error(f"Error calculating signature: {str(e)}")
            return None

    def calculate_mutation_context(self, df: pd.DataFrame) -> Dict:
        """Calculate mutation context percentages"""
        try:
            # Filter for SNPs
            snp_df = df[df['Variant_Type'] == 'SNP'].copy()
            
            # Process reference and alternate alleles
            ref_alleles = snp_df['Reference_Allele'].fillna('-')
            alt_alleles = snp_df['Tumor_Seq_Allele2'].fillna('-')
            
            # Calculate substitutions
            substitutions = []
            for ref, alt in zip(ref_alleles, alt_alleles):
                if len(ref) == 1 and len(alt) == 1:
                    substitutions.append(f"{ref}>{alt}")
            
            # Calculate percentages
            total = len(substitutions)
            percentages = []
            for context in self.CONTEXTS:
                sub_type = context.split('[')[1].split(']')[0]
                count = sum(1 for s in substitutions if s == sub_type)
                percentages.append((count/total * 100) if total > 0 else 0)
            
            return {
                'contexts': self.CONTEXTS,
                'percentages': percentages,
                'colors': [self.COLORS[ctx.split('[')[1].split(']')[0]] * 16 
                          for ctx in self.CONTEXTS]
            }
            
        except Exception as e:
            self.logger.error(f"Error calculating mutation context: {str(e)}")
            return None

    def plot_mutation_signature(self, signature_data: Dict, output_file: str):
        """Create mutation signature plot"""
        try:
            plt.figure(figsize=(10, 5))
            
            plt.bar(range(len(signature_data['contexts'])), 
                   signature_data['percentages'],
                   color=signature_data['colors'])
            
            plt.xticks(range(len(signature_data['contexts'])), 
                      signature_data['contexts'], 
                      rotation=90, fontsize=6)
            
            plt.ylabel("% of base substitutions")
            plt.title("Mutation Signature")
            
            plt.tight_layout()
            plt.savefig(output_file)
            plt.close()
            
        except Exception as e:
            self.logger.error(f"Error plotting mutation signature: {str(e)}")

class MAFSummaryGenerator:
    """Generate MAF summaries for clinical reports"""
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        
    def generate_summary_plots(self, analysis_data: Dict, sample_name: str):
        """Generate summary plots for clinical report"""
        try:
            if not analysis_data.get('basic_metrics'):
                return
                
            with PdfPages(f'data/results/reports/{sample_name}_summary.pdf') as pdf:
                # Plot 1: Variant Types
                plt.figure(figsize=(8, 6))
                variant_types = analysis_data['basic_metrics']['variant_types']
                plt.bar(list(variant_types.keys()), list(variant_types.values()))
                plt.title('Variant Type Distribution')
                plt.xticks(rotation=45)
                plt.tight_layout()
                pdf.savefig()
                plt.close()
                
                # Plot 2: Gene Metrics
                if 'gene_metrics' in analysis_data:
                    plt.figure(figsize=(10, 6))
                    top_genes = analysis_data['gene_metrics']['top_mutated_genes']
                    top_10 = {k: v for k, v in list(top_genes.items())[:10]}
                    plt.barh(list(top_10.keys()), list(top_10.values()))
                    plt.title('Top 10 Mutated Genes')
                    plt.tight_layout()
                    pdf.savefig()
                    plt.close()
                
        except Exception as e:
            self.logger.error(f"Error generating summary plots: {str(e)}")

    def create_summary_report(self, analysis_data: Dict, sample_name: str) -> bool:
        """Create summary report"""
        try:
            output_file = f'data/results/reports/{sample_name}_report.json'
            
            summary = {
                'sample_name': sample_name,
                'total_variants': analysis_data.get('basic_metrics', {}).get('total_variants', 0),
                'variant_types': len(analysis_data.get('basic_metrics', {}).get('variant_types', {})),
                'total_genes': analysis_data.get('gene_metrics', {}).get('genes_mutated', 0)
            }
            
            with open(output_file, 'w') as f:
                json.dump(summary, f, indent=4)
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error creating summary report: {str(e)}")
            return False

class MAFVisualizer:
    """Generate enhanced visualizations from MAF analysis results"""
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        # Set style parameters
        # plt.style.use('seaborn-v0_8-darkgrid')
        sns.set_theme(style="whitegrid")
        self.colors = ['#3498db', '#2ecc71', '#e74c3c', '#f1c40f', '#9b59b6']
        
    def _setup_figure_style(self, fig, ax):
        """Apply consistent styling to figures"""
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', linestyle='--', alpha=0.2, color='gray')
        fig.patch.set_facecolor('white')
        ax.set_facecolor('#f8f9fa')
        for spine in ax.spines.values():
            spine.set_color('#cccccc')

    def _convert_to_png(self, output_file: str) -> str:
        """Convert file path to PNG"""
        return str(Path(output_file).with_suffix('.png'))

class MAFVisualizer:
    """Generate enhanced visualizations from MAF analysis results"""
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        # Set style parameters
        sns.set_theme(style="whitegrid")
        self.colors = ['#3498db', '#2ecc71', '#e74c3c', '#f1c40f', '#9b59b6']
        
    def _setup_figure_style(self, fig, ax):
        """Apply consistent styling to figures"""
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', linestyle='--', alpha=0.2, color='gray')
        fig.patch.set_facecolor('white')
        ax.set_facecolor('#f8f9fa')
        for spine in ax.spines.values():
            spine.set_color('#cccccc')
            
    def _convert_to_png(self, output_file: str) -> str:
        """Convert file path to PNG"""
        return str(Path(output_file).with_suffix('.png'))

    def plot_mutation_signature(self, signature_data: Dict, output_file: str):
        """Create enhanced mutation signature visualization with base substitutions"""
        try:
            if not signature_data:
                return
                
            # Convert output_file to string and handle extension
            output_file = str(output_file)
            png_output = output_file.replace('.pdf', '.png') if '.pdf' in output_file else output_file
                
            # Extract data and prepare for plotting
            percentages = np.array(signature_data['percentages'])
            contexts = signature_data['contexts']
            
            # Create mapping of substitution types to colors
            color_mapping = {
                'C>A': '#87CEEB',  # Sky blue
                'C>G': '#2F4F4F',  # Dark slate gray  
                'C>T': '#CD5C5C',  # Indian red
                'T>A': '#808080',  # Gray
                'T>C': '#228B22',  # Forest green
                'T>G': '#FF69B4'   # Hot pink
            }

            # Group data by context
            data = {}
            for i, ctx in enumerate(contexts):
                full_context = ctx  # Keep full trinucleotide context
                sub_type = ctx.split('[')[1].split(']')[0]
                if sub_type not in data:
                    data[sub_type] = {'contexts': [], 'percentages': [], 'color': color_mapping[sub_type]}
                data[sub_type]['contexts'].append(full_context)
                data[sub_type]['percentages'].append(percentages[i])

            # Create figure
            plt.figure(figsize=(15, 8), dpi=300)
            
            # Plot bars for each substitution type
            x_pos = 0
            xticks_pos = []
            xticks_labels = []
            
            for sub_type, sub_data in data.items():
                positions = range(x_pos, x_pos + len(sub_data['contexts']))
                plt.bar(positions, sub_data['percentages'], 
                    color=sub_data['color'], 
                    alpha=0.8,
                    width=0.8)
                
                # Store positions for xticks
                xticks_pos.extend(positions)
                xticks_labels.extend(sub_data['contexts'])
                
                # Add substitution type label
                mid_pos = x_pos + len(sub_data['contexts'])/2 - 0.5
                plt.text(mid_pos, -1, sub_type, 
                        ha='center', va='top', 
                        fontsize=12, fontweight='bold')
                
                x_pos += len(sub_data['contexts'])
                
                # Add separator
                if sub_type != list(data.keys())[-1]:
                    plt.axvline(x_pos - 0.2, color='gray', linestyle='--', alpha=0.3)

            # Customize plot
            plt.xticks(xticks_pos, xticks_labels, rotation=45, ha='right', fontsize=8)
            plt.ylabel('% of Base Substitutions', fontsize=12)
            plt.title('Mutation Signature Profile', fontsize=14, pad=20)
            
            # Add grid and adjust layout
            plt.grid(axis='y', linestyle='--', alpha=0.3)
            plt.gca().spines['top'].set_visible(False)
            plt.gca().spines['right'].set_visible(False)
            
            # Add mean line for each substitution type
            for sub_type, sub_data in data.items():
                mean_val = np.mean(sub_data['percentages'])
                start_pos = xticks_pos[xticks_labels.index(sub_data['contexts'][0])]
                end_pos = xticks_pos[xticks_labels.index(sub_data['contexts'][-1])]
                plt.hlines(mean_val, start_pos - 0.4, end_pos + 0.4, 
                        colors='black', linestyles='dashed', linewidth=1)
                plt.text(end_pos + 0.5, mean_val, f'{mean_val:.1f}%', 
                        va='center', fontsize=8)

            # Add legend explaining the plot
            legend_text = (
                'Bar: Individual trinucleotide context frequency\n'
                'Dashed line: Mean frequency per substitution type\n'
                'Substitution types shown below contexts'
            )
            plt.figtext(0.99, 0.98, legend_text, 
                    ha='right', va='top', 
                    fontsize=8, 
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.8))

            plt.tight_layout()
            
            # Adjust bottom margin for substitution type labels
            plt.subplots_adjust(bottom=0.15)
            
            # Save as PNG
            plt.savefig(png_output, dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()
                
        except Exception as e:
            self.logger.error(f"Error plotting mutation signature: {str(e)}")   
            
    def plot_summary_dashboard(self, df: pd.DataFrame, analysis_results: Dict, output_file: str):
        """Create enhanced summary dashboard with multiple plots"""
        try:
            fig = plt.figure(figsize=(15, 12), dpi=300)
            gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
            
            # Plot 1: Variant Type Distribution
            ax1 = fig.add_subplot(gs[0, 0])
            variant_counts = df['Variant_Type'].value_counts()
            df_counts = pd.DataFrame({'Variant_Type': variant_counts.index, 'count': variant_counts.values})
            sns.barplot(data=df_counts, x='count', y='Variant_Type', 
                       hue='Variant_Type', palette=self.colors[:len(variant_counts)],
                       legend=False, ax=ax1)
            ax1.set_title('Variant Type Distribution', pad=10, fontsize=12, fontweight='bold')
            self._setup_figure_style(fig, ax1)
            
            # Plot 2: Variant Classification
            ax2 = fig.add_subplot(gs[0, 1])
            if 'Variant_Classification' in df.columns:
                class_counts = df['Variant_Classification'].value_counts()
                df_class = pd.DataFrame({'Classification': class_counts.index, 'count': class_counts.values})
                sns.barplot(data=df_class, x='count', y='Classification',
                          hue='Classification', 
                          palette=sns.color_palette("husl", len(class_counts)),
                          legend=False, ax=ax2)
                ax2.set_title('Variant Classification', pad=10, fontsize=12, fontweight='bold')
            self._setup_figure_style(fig, ax2)
            
            # Plot 3: Top Mutated Genes
            ax3 = fig.add_subplot(gs[1, 0])
            gene_counts = df['Hugo_Symbol'].value_counts().head(10)
            df_genes = pd.DataFrame({'Gene': gene_counts.index, 'count': gene_counts.values})
            sns.barplot(data=df_genes, x='count', y='Gene',
                       hue='Gene', palette=sns.color_palette("viridis", len(gene_counts)),
                       legend=False, ax=ax3)
            ax3.set_title('Top 10 Mutated Genes', pad=10, fontsize=12, fontweight='bold')
            self._setup_figure_style(fig, ax3)
            
            # Plot 4: Chromosome Distribution
            ax4 = fig.add_subplot(gs[1, 1])
            chrom_counts = df['Chromosome'].value_counts()
            df_chrom = pd.DataFrame({'Chromosome': chrom_counts.index, 'count': chrom_counts.values})
            sns.barplot(data=df_chrom, x='count', y='Chromosome',
                       hue='Chromosome', palette=sns.color_palette("mako", len(chrom_counts)),
                       legend=False, ax=ax4)
            ax4.set_title('Variants by Chromosome', pad=10, fontsize=12, fontweight='bold')
            self._setup_figure_style(fig, ax4)
            
            plt.savefig(self._convert_to_png(output_file), dpi=300, bbox_inches='tight', facecolor='white')
            plt.close(fig)
            
        except Exception as e:
            self.logger.error(f"Error creating summary dashboard: {str(e)}")

    def plot_comparison_summary(self, comparison_results: Dict, name1: str, name2: str, output_file: str):
        """Plot enhanced comparison summary"""
        try:
            # Close any existing figures
            plt.close('all')
            
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7), dpi=300)
            fig.patch.set_facecolor('white')
            
            colors = ['#2ecc71', '#3498db', '#9b59b6']
            
            # Plot 1: Variant Distribution
            variant_data = {
                f'Shared': comparison_results['shared_variants']['count'],
                f'Unique to {name1}': comparison_results[f'unique_to_{name1}']['count'],
                f'Unique to {name2}': comparison_results[f'unique_to_{name2}']['count']
            }
            
            wedges, texts, autotexts = ax1.pie(
                variant_data.values(),
                labels=variant_data.keys(),
                colors=colors,
                autopct='%1.1f%%',
                pctdistance=0.85,
                explode=(0.05, 0, 0),
                shadow=True
            )
            
            plt.setp(autotexts, size=9, weight="bold")
            plt.setp(texts, size=10)
            ax1.set_title(f'Variant Distribution\n{name1} vs {name2}', 
                         pad=20, size=12, weight='bold')
            
            # Plot 2: Gene Comparison
            categories = [f'Shared\nGenes', f'Unique to\n{name1}', f'Unique to\n{name2}']
            values = [
                comparison_results['gene_comparison']['shared_genes'],
                comparison_results['gene_comparison'][f'unique_to_{name1}'],
                comparison_results['gene_comparison'][f'unique_to_{name2}']
            ]
            
            bars = ax2.bar(categories, values, color=colors)
            ax2.set_title('Gene Distribution Comparison', 
                         pad=20, size=12, weight='bold')
            ax2.set_ylabel('Number of Genes', size=10)
            
            # Enhanced bar labels
            for bar in bars:
                height = bar.get_height()
                ax2.text(bar.get_x() + bar.get_width()/2., height,
                        f'{int(height):,}',
                        ha='center', va='bottom', size=10,
                        bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
            
            ax2.spines['top'].set_visible(False)
            ax2.spines['right'].set_visible(False)
            ax2.tick_params(labelsize=9)
            ax2.grid(axis='y', linestyle='--', alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(self._convert_to_png(output_file), dpi=300, bbox_inches='tight', facecolor='white')
            plt.close(fig)
            
        except Exception as e:
            self.logger.error(f"Error plotting comparison summary: {str(e)}")

    def plot_variant_distribution(self, df: pd.DataFrame, output_file: str):
        """Plot enhanced variant type distribution"""
        try:
            plt.close('all')  # Close any existing figures
            
            fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
            fig.patch.set_facecolor('white')
            ax.set_facecolor('#f8f9fa')
            
            variant_counts = df['Variant_Type'].value_counts()
            df_variants = pd.DataFrame({'Variant_Type': variant_counts.index, 'count': variant_counts.values})
            
            # Use seaborn with proper hue parameter
            sns.barplot(data=df_variants, x='Variant_Type', y='count',
                       hue='Variant_Type', palette=self.colors[:len(variant_counts)],
                       legend=False, ax=ax)
            
            # Add value labels
            for i, v in enumerate(variant_counts.values):
                ax.text(i, v, str(v), ha='center', va='bottom')
            
            ax.set_title('Variant Type Distribution', pad=20, size=12, weight='bold')
            ax.set_ylabel('Number of Variants', size=10)
            ax.set_xlabel('Variant Type', size=10)
            plt.xticks(rotation=45)
            
            self._setup_figure_style(fig, ax)
            
            plt.tight_layout()
            plt.savefig(self._convert_to_png(output_file), dpi=300, bbox_inches='tight', facecolor='white')
            plt.close(fig)
            
        except Exception as e:
            self.logger.error(f"Error plotting variant distribution: {str(e)}")

class MAFAnalysisPipeline:
    """Main analysis pipeline integrating all components"""
    def __init__(self):
        self.processor = MAFProcessor()
        self.signature_analyzer = MutationSignatureAnalyzer()
        self.visualizer = MAFVisualizer()
        self.summary_generator = MAFSummaryGenerator()
        self.logger = logging.getLogger(__name__)
        
        # Setup directory structure
        self.base_dir = Path('data')
        self.maf_dir = self.base_dir / 'maf'
        self.results_dir = self.base_dir / 'results'
        self.viz_dir = self.results_dir / 'visualizations'
        self.reports_dir = self.results_dir / 'reports'
        
        for directory in [self.maf_dir, self.results_dir, self.viz_dir, self.reports_dir]:
            directory.mkdir(exist_ok=True, parents=True)

    def process_sample(self, maf_file: Path) -> Dict:
        """Process a single MAF file with comprehensive analysis"""
        try:
            sample_name = maf_file.stem
            self.logger.info(f"Processing {sample_name}")
            
            # Read and validate MAF file
            df = self.processor.read_maf(str(maf_file))
            if df is None:
                return {}
            
            # Calculate basic metrics
            basic_metrics = {
                'total_variants': len(df),
                'variant_types': df['Variant_Type'].value_counts().to_dict(),
                'variant_classification': df['Variant_Classification'].value_counts().to_dict() 
                                       if 'Variant_Classification' in df.columns else {},
                'by_chromosome': df['Chromosome'].value_counts().to_dict()
            }
            
            # Calculate VAF statistics if available
            if all(col in df.columns for col in ['t_depth', 't_alt_count']):
                df['VAF'] = df['t_alt_count'] / df['t_depth']
                basic_metrics['vaf_stats'] = {
                    'mean': df['VAF'].mean(),
                    'median': df['VAF'].median(),
                    'std': df['VAF'].std()
                }
            
            # Calculate gene metrics
            gene_metrics = {
                'genes_mutated': df['Hugo_Symbol'].nunique(),
                'top_mutated_genes': df['Hugo_Symbol'].value_counts().head(20).to_dict()
            }
            
            # Calculate mutation signature
            mutation_signature = self.signature_analyzer.calculate_signature(df)
            
            # Generate additional mutation context analysis
            mutation_context = self.signature_analyzer.calculate_mutation_context(df)
            
            # Compile complete analysis results
            analysis_results = {
                'sample_name': sample_name,
                'basic_metrics': basic_metrics,
                'gene_metrics': gene_metrics,
                'mutation_signature': mutation_signature,
                'mutation_context': mutation_context
            }
            
            # Generate visualizations
            # if mutation_signature:
            #     self.visualizer.plot_mutation_signature(
            #         mutation_signature,
            #         self.viz_dir / f"{sample_name}_mutation_signature.pdf"
            #     )
            
            self.visualizer.plot_summary_dashboard(
                df, analysis_results,
                self.viz_dir / f"{sample_name}_summary_dashboard.pdf"
            )
            
            # Generate summary plots for clinical report
            self.summary_generator.generate_summary_plots(
                analysis_results, 
                sample_name
            )
            
            # Create summary report
            self.summary_generator.create_summary_report(
                analysis_results,
                sample_name
            )
            
            return analysis_results
            
        except Exception as e:
            self.logger.error(f"Error processing {maf_file}: {str(e)}")
            return {}

    def compare_samples(self, df1: pd.DataFrame, df2: pd.DataFrame, 
                    name1: str, name2: str) -> Dict:
        """Compare two samples with enhanced analysis"""
        try:
            def create_variant_id(row):
                return f"{row['Chromosome']}_{row['Start_Position']}_" \
                    f"{row['Reference_Allele']}_{row['Tumor_Seq_Allele2']}"
            
            # First generate variant IDs
            df1['variant_id'] = df1.apply(create_variant_id, axis=1)
            df2['variant_id'] = df2.apply(create_variant_id, axis=1)
            
            # Get variant sets
            shared = set(df1['variant_id']) & set(df2['variant_id'])
            unique_1 = set(df1['variant_id']) - set(df2['variant_id'])
            unique_2 = set(df2['variant_id']) - set(df1['variant_id'])
            
            # Get gene sets
            genes1 = set(df1['Hugo_Symbol'].unique())
            genes2 = set(df2['Hugo_Symbol'].unique())
            shared_genes = genes1 & genes2
            unique_genes1 = genes1 - genes2
            unique_genes2 = genes2 - genes1
            
            # Create the comparison results
            comparison_results = {
                'shared_variants': {
                    'count': len(shared),
                    'variants': list(shared)
                },
                f'unique_to_{name1}': {
                    'count': len(unique_1),
                    'variants': list(unique_1)
                },
                f'unique_to_{name2}': {
                    'count': len(unique_2),
                    'variants': list(unique_2)
                },
                'gene_comparison': {
                    'shared_genes': len(shared_genes),
                    f'unique_to_{name1}': len(unique_genes1),
                    f'unique_to_{name2}': len(unique_genes2)
                }
            }

            # Generate individual variant distributions
            self.visualizer.plot_variant_distribution(
                df1, self.viz_dir / f"{name1}_variant_distribution.pdf"
            )
            self.visualizer.plot_variant_distribution(
                df2, self.viz_dir / f"{name2}_variant_distribution.pdf"
            )
            
            # Generate comparison visualization
            self.visualizer.plot_comparison_summary(
                comparison_results,
                name1,
                name2,
                str(self.viz_dir / f"comparison_{name1}_{name2}.pdf")
            )
            
            return comparison_results
            
        except Exception as e:
            self.logger.error(f"Error comparing samples: {str(e)}")
            return {}

    def run_pipeline(self):
        """Run the complete analysis pipeline"""
        try:
            # Find all MAF files
            maf_files = list(self.maf_dir.glob('*.maf'))
            if not maf_files:
                self.logger.warning("No MAF files found")
                return
            
            results = {
                'total_samples': len(maf_files),
                'sample_analyses': {},
                'comparisons': {}
            }
            
            # Process each file
            samples = {}
            for maf_file in maf_files:
                df = self.processor.read_maf(str(maf_file))
                if df is not None:
                    samples[maf_file.stem] = df
                    results['sample_analyses'][maf_file.stem] = self.process_sample(maf_file)
            
            # Perform pairwise comparisons
            sample_names = list(samples.keys())
            for i in range(len(sample_names)):
                for j in range(i+1, len(sample_names)):
                    name1, name2 = sample_names[i], sample_names[j]
                    comparison = self.compare_samples(
                        samples[name1], samples[name2],
                        name1, name2
                    )
                    results['comparisons'][f"{name1}_vs_{name2}"] = comparison
            
            # Save complete results
            with open(self.results_dir / 'complete_analysis.json', 'w') as f:
                json.dump(results, f, indent=4, default=str)
            
            # Generate final summary report
            self.logger.info("\nAnalysis Summary:")
            for sample, analysis in results['sample_analyses'].items():
                self.logger.info(f"\n{sample}:")
                if 'basic_metrics' in analysis:
                    metrics = analysis['basic_metrics']
                    self.logger.info(f"  Total variants: {metrics['total_variants']}")
                    self.logger.info(f"  Variant types: {len(metrics['variant_types'])}")
                    if 'vaf_stats' in metrics:
                        self.logger.info(f"  Mean VAF: {metrics['vaf_stats']['mean']:.3f}")
                
                if 'gene_metrics' in analysis:
                    self.logger.info(f"  Total genes affected: {analysis['gene_metrics']['genes_mutated']}")
            
            # Log comparison results
            self.logger.info("\nSample Comparisons:")
            for comparison_name, comparison in results['comparisons'].items():
                self.logger.info(f"\n{comparison_name}:")
                self.logger.info(f"  Shared variants: {comparison['shared_variants']['count']}")
                self.logger.info(f"  Shared genes: {comparison['gene_comparison']['shared_genes']}")
                
        except Exception as e:
            self.logger.error(f"Pipeline error: {str(e)}")
            return None
        
        return results

def main():
    try:
        pipeline = MAFAnalysisPipeline()
        results = pipeline.run_pipeline()
        
    except Exception as e:
        logging.error(f"Error in main execution: {str(e)}")
        return 1
        
    return 0

if __name__ == "__main__":
    sys.exit(main())
