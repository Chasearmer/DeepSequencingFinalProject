import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_m6a_histogram(names, filepaths):

    # Choose the column for the histogram
    m6a_column = 'm6A_ESC_mm_num'
    gene_type_column = 'gene_type'
    colors = ['salmon', 'skyblue']

    # Prepare plot
    plt.figure(figsize=(10, 6))

    # Loop through each name for overlaying plots
    for i, name in enumerate(names):

        # Read in the data
        data = pd.read_csv(filepaths[i])

        # Modify values greater than 6
        data.loc[data[m6a_column] > 5, m6a_column] = 6  # Replace values greater than 6 with 7

        # Create the counts for each bin
        counts = data[m6a_column].value_counts().sort_index()
        percents = counts / counts.sum()

        # Adjusting the position of the bars for clarity
        bar_positions = np.arange(len(counts)) + 0.2 * (i - 1)

        # Create a bar plot
        plt.bar(bar_positions, percents, width=0.2, label=name, edgecolor='black', color=colors[i])

    # Customize the x-axis to show '5+' label and proper bin labels
    plt.xticks(np.arange(len(counts)), labels=[0, 1, 2, 3, 4, 5, '5+'])

    # Add labels, title, and legend
    plt.xlabel('m6A Count')
    plt.ylabel('Counts')
    plt.title('Comparison of m6A Counts')
    plt.legend()

    # Show the plot
    plt.savefig('comparison_histogram_housekeeping_frequency.png')

plot_m6a_histogram(['Housekeepers', 'Non-housekeepers'], ['datasets/housekeepers.csv', 'datasets/non_housekeepers.csv'])


