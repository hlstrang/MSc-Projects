import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

#reading in input file
def read_input(csv_file):
    try:
        df = pd.read_csv(csv_file)
        print(f'{"-" * 70}\nFile successfully read\n{"-" * 70}')
        return df
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)
        return None #debugging statement to ensure the file is a csv file

#before QC stats - mean, sd, median, min, max values for both confidence scores alongside NofPmids and NofSnps
def compute_stats(df, label="Before quality control"):
    print(f'{label}:')
    selected_cols = [
        'rna2locus_conf_score',
        'gene2disease_conf_score',
        'NofPmids',
        'NofSnps'
    ] #selects columns of interest
    subset_df = df[selected_cols] #creates table with only these columns
    stats_table = pd.DataFrame({
        'Mean': subset_df.mean(),
        'Standard Deviation': subset_df.std(),
        'Median': subset_df.median(),
        'Minimum Value': subset_df.min(),
        'Maximum Value': subset_df.max()
    }) #creates data frame with these columns and their associated statistics
    stats_table = stats_table.T.round(2) #transposes the table so stats are rows and variables are columns and rounds the numbers to 2 significant figures
    stats_table = stats_table.rename(columns={
        'rna2locus_conf_score': 'RNA to Locus Confidence Score',
        'gene2disease_conf_score': 'Gene to Disease Confidence Score',
        'NofPmids': 'No. of PMIDs',
        'NofSnps': 'No. of SNPs'
    }) #renames the columns in the table for better readability
    print(stats_table)
    return stats_table

#filtering data based on confidence score
def filter_data(df, rna2locus, gene2disease):
    filtered_df = df[
        (df['rna2locus_conf_score'] > rna2locus) &
        (df['gene2disease_conf_score'] > gene2disease)
    ]
    print(f'{"-" * 70}\nData successfully filtered.\nMinimum RNA to Locus score = {rna2locus}\nMinimum Gene to Disease score = {gene2disease}\n{"-" * 70}')
    return filtered_df #deletes rows that are below the threshold values set by user

#stats after QC
def compute_QC_stats(filtered_df, label="After quality control"):
    print(f'{label}:')
    selected_cols = [
        'rna2locus_conf_score',
        'gene2disease_conf_score',
        'NofPmids',
        'NofSnps'
    ]
    subset_QCdf = filtered_df[selected_cols]
    stats_tableQC = pd.DataFrame({
        'Mean': subset_QCdf.mean(),
        'Standard Deviation': subset_QCdf.std(),
        'Median': subset_QCdf.median(),
        'Minimum Value': subset_QCdf.min(),
        'Maximum Value': subset_QCdf.max()
    })
    stats_tableQC = stats_tableQC.T.round(2)
    stats_tableQC = stats_tableQC.rename(columns={
        'rna2locus_conf_score': 'RNA to Locus Confidence Score',
        'gene2disease_conf_score': 'Gene to Disease Confidence Score',
        'NofPmids': 'No. of PMIDs',
        'NofSnps': 'No. of SNPs'
    })
    print(stats_tableQC)
    return stats_tableQC

#top 3 organisms
def top_organisms(filtered_df, n=3):
    organisms = filtered_df['Organism'].value_counts().head(n)
    print(f'{"-" * 70}\nTop {n} Organisms by count:\n{organisms}')
    return organisms.index.tolist() #prints the top 3 organisms by count and returns it to a list

#plotting before QC
def plot_before_QC(stats_table, save_plots=False):
    stats_to_plot = ['Mean', 'Median', 'Minimum Value', 'Maximum Value']
    variables = ['RNA to Locus Confidence Score', 'Gene to Disease Confidence Score', 'No. of PMIDs','No. of SNPs']
    comparison = stats_table.loc[stats_to_plot, variables]
    ax = comparison.plot(kind='bar', figsize=(10, 6))
    ax.set_title("Log-Scaled Summary Statistics Before QC")
    ax.set_ylabel("Value")
    ax.set_xlabel("Statistic")
    ax.set_yscale('log')
    ax.set_xticks(range(len(stats_to_plot)))
    ax.set_xticklabels(stats_to_plot, rotation=0)
    plt.tight_layout()
    if save_plots:
        plt.savefig("before_QC.png")
        print(f"{'-'*70}\nPlot saved to before_QC.png\n{'-'*70}")
    else:
        plt.show()

#plotting after QC
def plot_after_QC(stats_tableQC, save_plots=False):
    stats_to_plot = ['Mean', 'Median', 'Minimum Value', 'Maximum Value']
    variables = ['RNA to Locus Confidence Score', 'Gene to Disease Confidence Score','No. of PMIDs','No. of SNPs']
    comparison = stats_tableQC.loc[stats_to_plot, variables]
    ax = comparison.plot(kind='bar', figsize=(10, 6))
    ax.set_title("Log-Scaled Summary Statistics After QC")
    ax.set_ylabel("Value")
    ax.set_xlabel("Statistic")
    ax.set_yscale('log')
    ax.set_xticks(range(len(stats_to_plot)))
    ax.set_xticklabels(stats_to_plot, rotation=0)
    plt.tight_layout()
    if save_plots:
        plt.savefig("after_QC.png")
        print(f"{'-'*70}\nPlot saved to after_QC.png\n{'-'*70}")
    else:
        plt.show()

#output filtered data as tsv file
def output_file(filtered_df, tsv_output):
    filtered_df.to_csv(tsv_output, sep='\t', index=False)
    print(f"{'-'*70}\nFiltered data saved to {tsv_output}\n{'-'*70}")

def statistics_mirna(csv_file, tsv_output, rna2locus, gene2disease, save_plots=False): #main function for use in terminal
    df = read_input(csv_file)
    stats_before = compute_stats(df, "Before quality control")
    plot_before_QC(stats_before, save_plots)
    filtered_df = filter_data(df, rna2locus, gene2disease)
    stats_after = compute_QC_stats(filtered_df, "After quality control")
    plot_after_QC(stats_after, save_plots)
    orgs = top_organisms(filtered_df)
    output_file(filtered_df, tsv_output)

#allows the file to be used in command line
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Analyze miRNA statistics and generate graphs using quality controlled data."
    )
    parser.add_argument(
        "-i", "--input-csv",
        required=True,
        metavar="FILE",
        help="Input CSV file containing miRNA data."
    )
    parser.add_argument(
        "-o", "--output-tsv",
        required=True,
        metavar="FILE",
        help="Output TSV file for results after quality control."
    )
    parser.add_argument(
        "-r2l", "--rna2locus",
        type=float, required=True,
        metavar="FLOAT",
        help="Minimum RNA-to-locus confidence score (50–100)."
    )
    parser.add_argument(
        "-g2d", "--gene2disease",
        type=float, required=True,
        metavar="FLOAT",
        help="Minimum gene-to-disease confidence score (0.3–1.0)."
    )
    parser.add_argument(
        "-saveplots", "--save-plots",
        action="store_true",
        help="If present, save plots to PNG files instead of showing them."
    )
    args = parser.parse_args()
    statistics_mirna(
        args.input_csv,
        args.output_tsv,
        args.rna2locus,
        args.gene2disease,
        args.save_plots
    )

