import numpy as np
import pandas as pd
from typing import List, Tuple, Union, Dict
from pathlib import Path
from glob import glob
import logging
import json
import os, re
from pprint import pprint

root = Path(os.getcwd()).absolute()
data_path = root / 'data'
maf_path = data_path / 'maf'
# MAFParse
class MAFParse: 
    def __init__(self) -> None:
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)

        self.crucial_columns = [
            'Hugo_Symbol',            # Gene name
            'Chromosome',             # Chromosome location
            'Start_Position',         # Start position
            'End_Position',           # End position
            'Variant_Classification', # Type of variant
            'Variant_Type',          # SNP, INS, DEL, etc.
            'Reference_Allele',      # Reference allele
            'Tumor_Seq_Allele1',     # Tumor allele 1
            'Tumor_Seq_Allele2',     # Tumor allele 2
            'dbSNP_RS',              # Known variant ID
            'HGVSp_Short',           # Protein change
            'SIFT',                  # SIFT prediction
            'PolyPhen',              # PolyPhen prediction
            'IMPACT',                # Variant impact
            'Consequence'            # Variant consequence
        ]

    def read_maf(self, filepath: str) -> pd.DataFrame:
        """Read and validate MAF file"""
        try:
            df = pd.read_csv(filepath, sep='\t', comment='#', low_memory=False)
            
            # Validate minimum required columns
            min_required = ['Hugo_Symbol', 'Chromosome', 'Start_Position', 'Variant_Classification']
            missing_cols = [col for col in min_required if col not in df.columns]
            if missing_cols:
                raise ValueError(f"Missing required columns: {missing_cols}")
                
            self.logger.info(f"Successfully read MAF file: {filepath}")
            return df
            
        except Exception as e:
            self.logger.error(f"Error reading MAF file {filepath}: {str(e)}")
            raise
            
    def extract_crucial_cols(self, df:pd.DataFrame) -> pd.DataFrame:
        try:
            available_columns = [col for col in self.crucial_columns if col in df.columns]
            df_crucial = df[available_columns].copy()

            # Add derived columns
            df_crucial['Is_Novel'] = df['dbSNP_RS'].isna() # is/is not present in dbSNP database
            
            # clean specific columns 
            # check if 'HGVSp_Short' col else fill with 'Unknown'
            if 'HGVSp_Short' in df_crucial.columns:
                df_crucial['HGVSp_Short'] = df_crucial['HGVSp_Short'].fillna('Unknown')
                
            return df_crucial
        
        except Exception as e: # what is the best way to handle this?
            self.logger.error(f"Error extracting crucial info: {str(e)}")
            raise
            
    def process_file(self, filepath:str) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Process single MAF file and retunr both full/crucial dataframes"""
        try: 
            df_full = self.read_maf(filepath=filepath)
            df_crucial = self.extract_crucial_cols(df=df_full)
            return df_full, df_crucial
        except Exception as e:
            self.logger.error(f"Error processing MAF file: {str(e)}")
            raise

# MAFAnalyze
class MAFAnalyze:
    def __init__(self) -> None:
        self.logger = logging.getLogger(__name__)

    def _classify_clinical_significance(self, variant_info:pd.Series) -> str:
        """Classify variant clinical significance"""
        try:
            # check for high impact variants
            if variant_info['Variant_Classification'] in ['Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'Splice_Site']:
                return 'High'
            
            # check low impact variants
            if variant_info['Variant_Classification'] in ['Silent', 'Intron', "3'UTR", "5'UTR"]:
                return 'Low'

        except Exception as e:
            return 'Unknown'
    
    def get_variant_summary(self, df:pd.DataFrame, filename:str='') -> Dict:
        """get comprehensive variant summary"""
        try: 
            df['Clinical_Significance'] = df.apply(self._classify_clinical_significance, axis=1)
            
            summary = {
                'filename': filename,
                'total_variants': len(df),
                'variant_types': df['Variant_Type'].value_counts().to_dict(),
                'variant_classifications': df['Variant_Classification'].value_counts().to_dict(),
                'chromosomal_distribution': df['Chromosome'].value_counts().to_dict(),
                'novel_variants': int(df['Is_Novel'].sum()),
                'known_variants': int((~df['Is_Novel']).sum()),
                'clinical_significance': df['Clinical_Significance'].value_counts().to_dict()
            }
            
            # Add high impact variants detail
            high_impact = df[df['Clinical_Significance'] == 'High']
            summary['high_impact_variants'] = {
                'count': len(high_impact),
                'genes': high_impact['Hugo_Symbol'].unique().tolist()
            }
            
            return summary
            
        except Exception as e:
            self.logger.error(f"Error generating variant summary: {str(e)}")
            raise

    def compare_samples(self, df1:pd.DataFrame, df2:pd.DataFrame,
                        name1:str = 'Sample1', name2:str = 'Sample2') -> Dict:
        """Compare two MAF samples"""
        try:
            def create_variant_id(df):
                # create id like: chr1_12345_A_T
                return df.apply(lambda x: f"{x['Chromosome']}_{x['Start_Position']}_"f"{x['Reference_Allele']}_{x['Tumor_Seq_Allele2']}", axis=1)
            
            df1, df2 = df1.copy(), df2.copy()
            df1['Variant_ID'] = create_variant_id(df1)
            df2['Variant_ID'] = create_variant_id(df2)
        
            shared_variants = set(df1['Variant_ID']).intersection(set(df2['Variant_ID']))
            unique_to_1 = set(df1['Variant_ID']).difference(set(df2['Variant_ID']))
            unique_to_2 = set(df2['Variant_ID']).difference(set(df1['Variant_ID']))

            # clinical significance if not present
            if 'Clinical_Significance' not in df1.columns:
                df1['Clinical_Significance'] = df1.apply(self._classify_clinical_significance, axis=1)
            if 'Clinical_Significance' not in df2.columns:
                df2['Clinical_Significance'] = df2.apply(self._classify_clinical_significance, axis=1)
            
            comparison = {
                'shared_variants': {
                    'count': len(shared_variants),
                    'variants': list(shared_variants)
                },
                f'unique_to_{name1}': {
                    'count': len(unique_to_1),
                    'high_impact_genes': df1[
                        (df1['Variant_ID'].isin(unique_to_1)) & 
                        (df1['Clinical_Significance'] == 'High')
                    ]['Hugo_Symbol'].unique().tolist()
                },
                f'unique_to_{name2}': {
                    'count': len(unique_to_2),
                    'high_impact_genes': df2[
                        (df2['Variant_ID'].isin(unique_to_2)) & 
                        (df2['Clinical_Significance'] == 'High')
                    ]['Hugo_Symbol'].unique().tolist()
                }
            }
            
            return comparison

        except Exception as e:
            self.logger.error(f"Error comparing samples: {str(e)}")
            raise

def find_maf_pairs(dir: Path) -> List[Tuple[str, str]]:
    """
    Find paired MAF files based on common prefixes (MMCID or TCMK1)
    """
    try:
        # Get all MAF files
        files = list(dir.glob('*.maf'))
        
        # Group files by their prefix
        mmcid_files = [f for f in files if 'MMCID' in f.name]
        tcmk1_files = [f for f in files if 'TCMK1' in f.name]
        
        pairs = []
        
        # Pair MMCID files
        if len(mmcid_files) == 2:
            pairs.append((str(mmcid_files[0]), str(mmcid_files[1])))
            
        # Pair TCMK1 files
        if len(tcmk1_files) == 2:
            pairs.append((str(tcmk1_files[0]), str(tcmk1_files[1])))
            
        logging.info(f"Found {len(pairs)} pairs: {pairs}")
        return pairs
        
    except Exception as e:
        logging.error(f"Error finding MAF pairs: {str(e)}")
        return []


def analyze_maf_dir(dir:str) -> Dict:
    """analyze all maf files"""
    parse = MAFParse()
    analyzer = MAFAnalyze()
    results = {}

    # find paired MAF files
    pairs = find_maf_pairs(dir=dir)
    for f1, f2 in pairs:
        try:
            # parse files
            _, df_crucial1 = parse.process_file(filepath=f1)
            _, df_crucial2 = parse.process_file(filepath=f2)
            
            # generate summaries
            summary1 = analyzer.get_variant_summary(df=df_crucial1, filename=f1)
            summary2 = analyzer.get_variant_summary(df=df_crucial2, filename=f2)

            # comparison
            comparison = analyzer.compare_samples(
                df1=df_crucial1, df2=df_crucial2,
                name1=Path(f1).stem, name2=Path(f2).stem
            )

            pair_name = f'{Path(f1).stem}_{Path(f2).stem}'
            results[pair_name] = {
                'sample_summary1': summary1,
                'sample_summary2': summary2,
                'comparison': comparison
            }
        except Exception as e:
            logging.error(f"Error processing pair {f1}, {f2}: {str(e)}")
            continue
        return results

def save_json(data: Dict, filename: str):
    """
    Save a dictionary to a JSON file with indentation.

    Parameters:
    data (Dict): The dictionary to save.
    filename (str): The name of the file to save the dictionary to.
    """
    with open(filename, 'w') as file:
        json.dump(data, file, indent=4)



parse = MAFParse()
analyze = MAFAnalyze()

sample1 = maf_path/'MMCID-26B.maf'
sample2 = maf_path/'MMCID-30B.maf'
parse = MAFParse()

df1_full, df1_crucial = parse.process_file(filepath=sample1)
summary = analyze.get_variant_summary(df=df1_crucial, filename=sample1)

df2_full, df2_crucial = parse.process_file(filepath=sample2)
summary = analyze.get_variant_summary(df=df2_crucial, filename=sample2)

# display(df1_full, df1_crucial, summary)
compare = analyze.compare_samples(df1=df1_crucial, df2=df2_crucial, name1='MMCID-26B', name2='MMCID-30B')
pprint(compare)

results = analyze_maf_dir(dir=maf_path)
pprint(results)

save_path = data_path / 'results'
save_json(data=results, filename=save_path/'results.json')
save_json(data=compare, filename=save_path/'compare.json')