from scipy.stats import linregress
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Paths to the uploaded files
m3_ko_m6a_only_file_path = 'datasets/de_Mettl3_KO_vs_WT_m6A_m6a_only.csv'
ythdf_ko_m6a_only_file_path = 'datasets/de_YTHDF_Tri_KO_vs_WT_m6A_m6a_only.csv'

path_to_name = {
    m3_ko_m6a_only_file_path: "Mettle3 (m6A only)",
    ythdf_ko_m6a_only_file_path: "YTHDF (m6A only)",
}

experiments = [
    [m3_ko_m6a_only_file_path, ythdf_ko_m6a_only_file_path],
]

def log2_scatterplot(filepath_1, filepath_2):
    # Read the datasets
    first_ko = pd.read_csv(filepath_1)
    second_ko = pd.read_csv(filepath_2)

    # Get the names
    name_1 = path_to_name[filepath_1]
    name_2 = path_to_name[filepath_2]

    # Merge the datasets on gene_id
    merged_data = pd.merge(first_ko, second_ko, on='gene_id', suffixes=('_first', '_second'))

    # Split the dataset based on gene_type
    lncRNA_data = merged_data[
        (merged_data['gene_type_first'] == 'lncRNA') & (merged_data['gene_type_second'] == 'lncRNA')]
    protein_coding_data = merged_data[
        (merged_data['gene_type_first'] == 'protein_coding') & (merged_data['gene_type_second'] == 'protein_coding')]

    print("lncRNA Count: ", len(lncRNA_data))
    print("Protein Coding Count: ", len(protein_coding_data))

    for data, gene_type in zip([lncRNA_data, protein_coding_data], ['lncRNA', 'protein_coding']):
        # Drop rows with pvalue > 0.05
        data = data[(data['padj_first'] < 0.05) & (data['padj_second'] < 0.05)]


        # Perform linear regression
        slope, intercept, r_value, p_value, std_err = linregress(data['log2FoldChange_first'],
                                                                 data['log2FoldChange_second'])

        # Calculate the line of best fit
        line = slope * data['log2FoldChange_first'] + intercept

        # Calculate R squared value
        r_squared = r_value ** 2

        # Calculate sample size
        sample_size = len(data)

        # Create a figure with gridspec for custom layout
        fig = plt.figure(figsize=(12, 10))
        gs = fig.add_gridspec(2, 2, width_ratios=(5, 1), height_ratios=(1, 5), hspace=0.05, wspace=0.05)

        # Main scatter plot
        ax_scatter = fig.add_subplot(gs[1, 0])
        ax_scatter.scatter(data['log2FoldChange_first'], data['log2FoldChange_second'], alpha=0.5, color='plum')
        ax_scatter.plot(data['log2FoldChange_first'], line, color='red')  # line of best fit
        ax_scatter.set_xlabel(f'Log Fold Change ({name_1} KO)')
        ax_scatter.set_ylabel(f'Log Fold Change ({name_2} KO)')
        ax_scatter.annotate(f'Slope: {slope:.3f}', xy=(0.05, 0.95), xycoords='axes fraction',
                            fontsize=10, backgroundcolor='white',
                            horizontalalignment='left', verticalalignment='top')
        plt.grid()

        # Histogram for first KO
        ax_histx = fig.add_subplot(gs[0, 0], sharex=ax_scatter)
        ax_histx.hist(data['log2FoldChange_first'], bins=150, density=True, color='salmon')
        ax_histx.set_ylabel('Density')

        # Histogram for second KO (rotated 90 degrees)
        ax_histy = fig.add_subplot(gs[1, 1], sharey=ax_scatter)
        ax_histy.hist(data['log2FoldChange_second'], bins=150, density=True, orientation='horizontal', color='skyblue')
        ax_histy.set_xlabel('Density')

        # Set the title for the main scatter plot
        ax_histx.set_title(f'{gene_type} - Scatter Plot of Log Fold Change\n$R^2 = {r_squared:.3f}, N={sample_size}$')

        # Hide x labels and tick labels for top plots and y ticks for right plots
        plt.setp(ax_histx.get_xticklabels(), visible=False)
        plt.setp(ax_histy.get_yticklabels(), visible=False)

        # Save the plot
        plt.savefig(f"{name_1}_vs_{name_2}_{gene_type}_ko_scatterplot.png", dpi=300)
        plt.close()


for e in experiments:
    fpath1, fpath2 = e
    log2_scatterplot(fpath1, fpath2)
