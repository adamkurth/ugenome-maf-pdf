#!/usr/bin/env python3
import pandas as pd
import numpy as np
from pathlib import Path
import logging
import json
from typing import Dict, List, Tuple, Optional
import glob
import os

class MAFAnalyzer:
    def __init__(self):
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)
        
        # Define comprehensive column sets
        self.variant_info_columns = [
            'Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build',
            'Chromosome', 'Start_Position', 'End_Position', 'Strand',
            'Variant_Classification', 'Variant_Type', 'Reference_Allele',
            'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS',
            'dbSNP_Val_Status', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode'
        ]
        
        self.sequencing_info_columns = [
            'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2',
            'Tumor_Validation_Allele1', 'Tumor_Validation_Allele2',
            'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2',
            'Verification_Status', 'Validation_Status', 'Mutation_Status',
            'Sequencing_Phase', 'Sequence_Source', 'Validation_Method',
            'Score', 'BAM_File', 'Sequencer'
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

    def setup_directories(self) -> Tuple[Path, Path]:
        """Setup necessary directories for data and results"""
        base_dir = Path('data')
        maf_dir = base_dir / 'maf'
        results_dir = base_dir / 'results'
        
        for directory in [base_dir, maf_dir, results_dir]:
            directory.mkdir(exist_ok=True)
            
        return maf_dir, results_dir

    def read_maf(self, filepath: str) -> Optional[pd.DataFrame]:
        """Read MAF file with comprehensive error handling"""
        try:
            # Read the file and skip metadata lines starting with #
            df = pd.read_csv(filepath, sep='\t', comment='#', low_memory=False)
            
            # Basic validation
            required_cols = ['Hugo_Symbol', 'Chromosome', 'Start_Position', 'Variant_Type']
            missing = [col for col in required_cols if col not in df.columns]
            if missing:
                self.logger.warning(f"Missing required columns in {filepath}: {missing}")
            
            # Convert chromosome to string and ensure proper formatting
            df['Chromosome'] = df['Chromosome'].astype(str)
            
            # Handle missing values
            df = df.replace({'': np.nan, '.': np.nan})
            
            return df
            
        except Exception as e:
            self.logger.error(f"Error reading MAF file {filepath}: {str(e)}")
            return None

    def calculate_variant_metrics(self, df: pd.DataFrame) -> Dict:
        """Calculate comprehensive variant metrics"""
        metrics = {
            'variant_counts': {
                'total': len(df),
                'by_type': df['Variant_Type'].value_counts().to_dict(),
                'by_classification': df['Variant_Classification'].value_counts().to_dict(),
                'by_chromosome': df['Chromosome'].value_counts().to_dict()
            },
            'impact_distribution': df['IMPACT'].value_counts().to_dict() if 'IMPACT' in df.columns else {},
            'consequence_distribution': df['Consequence'].value_counts().to_dict() if 'Consequence' in df.columns else {},
            'clinical_significance': df['CLIN_SIG'].value_counts().to_dict() if 'CLIN_SIG' in df.columns else {}
        }
        
        # Calculate variant allele frequencies if coverage data is available
        if all(col in df.columns for col in ['t_depth', 't_alt_count']):
            df['VAF'] = df['t_alt_count'] / df['t_depth']
            metrics['vaf_statistics'] = {
                'mean': df['VAF'].mean(),
                'median': df['VAF'].median(),
                'std': df['VAF'].std()
            }
            
        return metrics

    def extract_gene_level_summary(self, df: pd.DataFrame) -> Dict:
        """Extract gene-level mutation summary"""
        gene_summary = {
            'most_mutated_genes': df['Hugo_Symbol'].value_counts().head(20).to_dict(),
            'genes_by_impact': {}
        }
        
        if 'IMPACT' in df.columns:
            impact_groups = df.groupby(['Hugo_Symbol', 'IMPACT']).size().unstack(fill_value=0)
            gene_summary['genes_by_impact'] = impact_groups.to_dict()
            
        return gene_summary

    def analyze_sample(self, filepath: str, results_dir: Path) -> Dict:
        """Analyze a single MAF file with comprehensive metrics"""
        df = self.read_maf(filepath)
        if df is None:
            return {}
            
        sample_name = Path(filepath).stem
        
        # Generate comprehensive analysis
        analysis = {
            'sample_name': sample_name,
            'file_stats': {
                'total_variants': len(df),
                'file_path': str(filepath),
                'creation_date': os.path.getctime(filepath)
            },
            'variant_metrics': self.calculate_variant_metrics(df),
            'gene_summary': self.extract_gene_level_summary(df)
        }
        
        # Save processed data
        processed_dir = results_dir / 'processed_data'
        processed_dir.mkdir(exist_ok=True)
        
        # Save different variant types separately
        for variant_type in df['Variant_Type'].unique():
            subset = df[df['Variant_Type'] == variant_type]
            subset.to_csv(processed_dir / f'{sample_name}_{variant_type.lower()}.csv', index=False)
        
        # Save complete processed data
        df.to_csv(processed_dir / f'{sample_name}_complete.csv', index=False)
        
        return analysis

    def process_all_samples(self) -> Dict:
        """Process all MAF files in the data directory"""
        maf_dir, results_dir = self.setup_directories()
        
        # Find all MAF files
        maf_files = glob.glob(str(maf_dir / '*.maf'))
        if not maf_files:
            self.logger.warning("No MAF files found in data/maf directory")
            return {}
            
        # Process each file
        results = {
            'sample_analyses': {},
            'comparative_analysis': {
                'total_samples': len(maf_files),
                'shared_variants': {},
                'unique_variants': {}
            }
        }
        
        # Analyze individual samples
        for maf_file in maf_files:
            sample_analysis = self.analyze_sample(maf_file, results_dir)
            results['sample_analyses'][Path(maf_file).stem] = sample_analysis
            
        # Save results
        with open(results_dir / 'complete_analysis.json', 'w') as f:
            json.dump(results, f, indent=4)
            
        self.logger.info(f"Analysis complete. Processed {len(maf_files)} MAF files.")
        return results

def main():
    analyzer = MAFAnalyzer()
    results = analyzer.process_all_samples()
    
    # Log summary statistics
    logging.info("Analysis Summary:")
    for sample, analysis in results['sample_analyses'].items():
        if analysis:  # Check if analysis exists
            variant_count = analysis.get('file_stats', {}).get('total_variants', 0)
            logging.info(f"{sample}: {variant_count} variants processed")

if __name__ == "__main__":
    main()


"""
RUN IN TERMINAL:
> python scripts/process-maf-data.py
2024-11-20 11:10:42,709 - INFO - Analysis complete. Processed 4 MAF files.
2024-11-20 11:10:42,709 - INFO - Analysis Summary:
2024-11-20 11:10:42,709 - INFO - MMCID-26B: 535 variants processed
2024-11-20 11:10:42,709 - INFO - MMCID-30B: 275 variants processed
2024-11-20 11:10:42,709 - INFO - TCMK1-14B: 327 variants processed
2024-11-20 11:10:42,709 - INFO - TMCK1-23B: 322 variants processed
"""