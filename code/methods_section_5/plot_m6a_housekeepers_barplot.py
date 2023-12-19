import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def plot_m6a_stacked_histogram(names, filepaths):
    # Choose the column for the histogram
    m6a_column = 'm6A_ESC_mm_num'
    colors = ['salmon', 'skyblue']

    # Prepare plot
    plt.figure(figsize=(10, 6))

    # Read in the data for both housekeepers and non-housekeepers
    housekeepers = pd.read_csv(filepaths[0])
    non_housekeepers = pd.read_csv(filepaths[1])

    # Modify values greater than 6
    housekeepers.loc[housekeepers[m6a_column] > 5, m6a_column] = 6
    non_housekeepers.loc[non_housekeepers[m6a_column] > 5, m6a_column] = 6

    num_bins = 7

    # Initialize a dictionary to hold the count data
    count_data = {i: {'housekeepers': 0, 'non_housekeepers': 0} for i in range(num_bins)}  # 0 to 6

    # Create the counts for each bin for housekeepers
    counts_housekeepers = housekeepers[m6a_column].value_counts().sort_index()
    for m6a_count in range(num_bins):  # 0 to 6
        count_data[m6a_count]['housekeepers'] = counts_housekeepers.get(m6a_count, 0)

    # Create the counts for each bin for non-housekeepers
    counts_non_housekeepers = non_housekeepers[m6a_column].value_counts().sort_index()
    for m6a_count in range(num_bins):  # 0 to 6
        count_data[m6a_count]['non_housekeepers'] = counts_non_housekeepers.get(m6a_count, 0)

    # Plot the bars
    bar_positions = np.arange(num_bins)
    for m6a_count in range(num_bins):
        # Calculate the total count for the current m6A bin
        total_count = count_data[m6a_count]['housekeepers'] + count_data[m6a_count]['non_housekeepers']
        # Calculate the proportion for housekeepers and non-housekeepers
        prop_housekeepers = count_data[m6a_count]['housekeepers'] / total_count if total_count else 0
        prop_non_housekeepers = count_data[m6a_count]['non_housekeepers'] / total_count if total_count else 0

        # Stack the proportions in the bar plot
        plt.bar(bar_positions[m6a_count], prop_housekeepers, width=0.4, label='Housekeepers' if m6a_count == 0 else "",
                color=colors[0])
        plt.bar(bar_positions[m6a_count], prop_non_housekeepers, width=0.4, bottom=prop_housekeepers,
                label='Non-housekeepers' if m6a_count == 0 else "", color=colors[1])

    # Customize the x-axis to show '6+' label and proper bin labels
    plt.xticks(bar_positions, labels=[0, 1, 2, 3, 4, 5, '5+'])

    # Add labels, title, and legend
    plt.xlabel('Number of m6A genes')
    plt.ylabel('Proportion within m6A category')
    plt.title('Proportional Comparison of m6A Gene Counts (Stacked)')
    plt.legend()

    # Show the plot
    plt.savefig('comparison_stacked_histogram.png')


plot_m6a_stacked_histogram(['Housekeepers', 'Non-housekeepers'],
                           ['datasets/housekeepers.csv', 'datasets/non_housekeepers.csv'])
