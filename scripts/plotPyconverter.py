# -*- coding: utf-8 -*-

"""
Created on Wed Nov  6 16:06:11 2024

@author: Joseph Rosen - joemrosen@gmail.com
"""
#Mouse Exome Seq MAF Plot Generation -> Py version

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
from matplotlib.backends.backend_pdf import PdfPages
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

os.chdir('mouseExome')

# Calculate percentages
MMCID_TNM_percentages = [ 2.45775729646697, 0.307219662058372, 0.460829493087558, 1.38248847926267, 2.30414746543779, 1.0752688172043, 0.921658986175115, 2.91858678955453, 2.1505376344086, 1.22887864823349, 0.153609831029186, 1.38248847926267, 0.614439324116744, 0.921658986175115, 0, 2.45775729646697, 0.768049155145929, 0.153609831029186, 0.307219662058372, 0.153609831029186, 0.307219662058372, 0.153609831029186, 0.460829493087558, 0.307219662058372, 0.460829493087558, 0.460829493087558, 0.614439324116744, 0.307219662058372, 0.307219662058372, 0, 0, 0.460829493087558, 1.68970814132104, 1.68970814132104, 3.68663594470046, 0.921658986175115, 2.30414746543779, 2.1505376344086, 5.22273425499232, 3.68663594470046, 1.99692780337942, 2.76497695852535, 5.06912442396313, 4.91551459293395, 0.921658986175115, 1.0752688172043, 2.76497695852535, 2.76497695852535, 0.307219662058372, 0, 0, 0.460829493087558, 0.614439324116744, 0.768049155145929, 1.22887864823349, 0, 0, 0.307219662058372, 0.921658986175115, 0, 1.68970814132104, 0, 0, 0.153609831029186, 0.921658986175115, 0.153609831029186, 0.921658986175115, 0.307219662058372, 0.921658986175115, 1.68970814132104, 1.22887864823349, 2.30414746543779, 0.768049155145929, 0.614439324116744, 2.30414746543779, 1.99692780337942, 0.768049155145929, 0.614439324116744, 1.38248847926267, 1.53609831029186, 0.153609831029186, 0.614439324116744, 0.153609831029186, 0.307219662058372, 0, 0.614439324116744, 0.307219662058372, 0.307219662058372, 0.153609831029186, 0.307219662058372, 1.22887864823349, 0.614439324116744, 0, 0.307219662058372, 0.460829493087558, 0 ]

MMCID_TNM_percentages_Names = ["A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T", "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T", "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T", "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T", "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T", "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T"]

TMCK1_TNM_percentages = [ 0.922509225092251, 1.1070110701107, 0.18450184501845, 0.3690036900369, 0.738007380073801, 0.922509225092251, 0.738007380073801, 1.4760147601476, 1.1070110701107, 0.3690036900369, 0.3690036900369, 1.8450184501845, 1.8450184501845, 1.66051660516605, 0.3690036900369, 3.690036900369, 1.29151291512915, 0.55350553505535, 0, 0.3690036900369, 0.738007380073801, 0.18450184501845, 0, 0.738007380073801, 0.738007380073801, 2.02952029520295, 0.55350553505535, 0.738007380073801, 0.55350553505535, 0, 0.18450184501845, 0.55350553505535, 2.2140221402214, 1.1070110701107, 2.02952029520295, 0.55350553505535, 1.29151291512915, 2.39852398523985, 0.18450184501845, 2.39852398523985, 2.76752767527675, 0.738007380073801, 5.53505535055351, 2.02952029520295, 0.738007380073801, 2.5830258302583, 1.29151291512915, 0.922509225092251, 0.738007380073801, 0, 0, 2.02952029520295, 0, 1.66051660516605, 0.738007380073801, 0, 0, 0, 0.18450184501845, 0, 1.1070110701107, 0, 0, 0.3690036900369, 1.1070110701107, 0, 1.29151291512915, 0.3690036900369, 0.922509225092251, 1.8450184501845, 2.02952029520295, 1.1070110701107, 1.1070110701107, 3.13653136531365, 2.39852398523985, 0.738007380073801, 0, 1.8450184501845, 1.66051660516605, 0.3690036900369, 0.3690036900369, 0.3690036900369, 2.02952029520295, 0.922509225092251, 0.18450184501845, 0, 1.66051660516605, 3.13653136531365, 0.3690036900369, 1.8450184501845, 2.2140221402214, 1.66051660516605, 0.55350553505535, 0.55350553505535, 0.738007380073801, 0.922509225092251 ]

TMCK1_TNM_percentages_Names = ["A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T", "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T", "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T", "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T", "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T", "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T"]
# Define colors
colors = ["skyblue"] * 16 + ["black"] * 16 + ["purple"] * 16 + ["grey"] * 16 + ["green"] * 16 + ["pink"] * 16

# Plot 1: MMCID
def make_MM_plot() :
    plt.figure(figsize=(10, 5))
    plt.bar(range(len(MMCID_TNM_percentages)), MMCID_TNM_percentages, color=colors)
    plt.xticks(range(len(MMCID_TNM_percentages)), MMCID_TNM_percentages_Names, rotation=90, fontsize=2)  # Use names for X-axis
    plt.ylabel("% of base substitutions")
    plt.title("MMCID Mutational Signature")
    plt.tight_layout()
    plt.savefig("MMCID_Mutational_Signature_Python.pdf")
    plt.show()
    plt.close()

    

#make_MM_plot()

# Plot 2: TMCK1
def make_TM_plot() : 
    plt.figure(figsize=(10, 5))
    plt.bar(range(len(TMCK1_TNM_percentages)), TMCK1_TNM_percentages, color=colors)
    plt.xticks(range(len(TMCK1_TNM_percentages)), TMCK1_TNM_percentages_Names, rotation=90, fontsize=2)  # Use names for X-axis
    plt.ylabel("% of base substitutions")
    plt.title("TMCK1 Mutational Signature")
    plt.tight_layout()
    plt.savefig("TMCK1_Mutational_Signature_Python.pdf")
    plt.show()
    plt.close()

#make_TM_plot()


#Creating MAF summary statistics plots


#reading .txt MAF data in
MMCID_maf_df = pd.read_table("MMCID_maf_data.txt", sep="\t")

TMCK1_maf_df = pd.read_table("TMCK1_maf_data.txt", sep="\t")


def get_MAF_summary(outputPDFname: str, maf_file_path: str) -> str:
    # Activate the automatic conversion between R data frames and Pandas DataFrames
    pandas2ri.activate()
    
    # Import the maftools package from R
    maftools = importr('maftools')
    base =importr('base')
    
    # Load the MAF file as an MAF object in R
    laml = maftools.read_maf(maf=maf_file_path)
    
    #Convert each summary to a plain R data frame before returning to Python
    
    sample_summary_r = base.as_data_frame(maftools.getSampleSummary(laml))
    gene_summary_r = base.as_data_frame(maftools.getGeneSummary(laml))
    variant_classification_r = base.as_data_frame(maftools.getFields(laml))
    
    # Convert R data frames to Pandas DataFrames
    sample_summary = pandas2ri.rpy2py(sample_summary_r)
    gene_summary = pandas2ri.rpy2py(gene_summary_r)
    variant_classification = pandas2ri.rpy2py(variant_classification_r)
    
        
    def plot_maf_summary(sample_summary, gene_summary, variant_classification, file_name=outputPDFname):
        # Set up a 2x2 grid of subplots
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # Plot 1: Variants per Sample (Stacked Barplot)
        sns.barplot(x='Tumor_Sample_Barcode', y='Mutation_Count', data=sample_summary, ax=axes[0, 0], palette='viridis')
        axes[0, 0].set_title('Number of Variants per Sample')
        axes[0, 0].set_xlabel('Sample')
        axes[0, 0].set_ylabel('Mutation Count')
        axes[0, 0].tick_params(axis='x', rotation=90)
    
        # Plot 2: Variant Classification Summary (Barplot)
        sns.barplot(x='Variant_Classification', y='Count', data=variant_classification, ax=axes[0, 1], palette='coolwarm')
        axes[0, 1].set_title('Variant Classification Summary')
        axes[0, 1].set_xlabel('Variant Classification')
        axes[0, 1].set_ylabel('Count')
        axes[0, 1].tick_params(axis='x', rotation=45)
    
        # Plot 3: Top Mutated Genes (Top 10 Mutated Genes)
        top_genes = gene_summary.nlargest(10, 'Mutation_Count')  # Top 10 genes by mutation count
        sns.barplot(y='Hugo_Symbol', x='Mutation_Count', data=top_genes, ax=axes[1, 0], palette='plasma')
        axes[1, 0].set_title('Top 10 Mutated Genes')
        axes[1, 0].set_xlabel('Mutation Count')
        axes[1, 0].set_ylabel('Gene')
        
        # Plot 4: Boxplot of Mutations by Variant Classification
        sns.boxplot(x='Variant_Classification', y='Mutation_Count', data=sample_summary, ax=axes[1, 1], palette='coolwarm')
        axes[1, 1].set_title('Boxplot of Mutation Counts by Variant Classification')
        axes[1, 1].set_xlabel('Variant Classification')
        axes[1, 1].set_ylabel('Mutation Count')
        axes[1, 1].tick_params(axis='x', rotation=45)
    
        # Adjust layout and save to a single PDF page
        plt.tight_layout()
        with PdfPages(file_name) as pdf:
            pdf.savefig(fig)
            plt.close(fig)
    
        print(f"Plots have been saved to {file_name}")

    # Call the plotting function with the extracted summaries
    plot_maf_summary(sample_summary, gene_summary, variant_classification, file_name=outputPDFname)
    
    # Return the name of the saved PDF
    return outputPDFname

get_MAF_summary('TMCK1_summary_plots_pyTSST.pdf', 'TMCK1_maf_data.txt')

# Example usage:
# output_pdf = get_MAF_summary("MAF_Summary_Plots.pdf", "MMCID_maf_data.txt")
        

        
        
