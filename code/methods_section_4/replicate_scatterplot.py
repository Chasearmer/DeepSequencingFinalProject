import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Path to the uploaded file
file_path = 'datasets/normalized_count_matrix_2.csv'  # Update this with the actual file path

# Dictionary mapping column names to sample names
column_to_sample = {
    "SRR11459699": "WT",
    "SRR11459700": "WT",
}

# Pairs of replicates for plotting
replicate_pairs = [
    ("SRR11459699", "SRR11459700"),  # WT
]

def remove_outliers(data, column_1, column_2):
    removal_factor = 50

    Q1 = data[column_1].quantile(0.25)
    Q3 = data[column_1].quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - removal_factor * IQR
    upper_bound = Q3 + removal_factor * IQR

    filtered_data = data[(data[column_1] >= lower_bound) & (data[column_1] <= upper_bound)]

    Q1 = filtered_data[column_2].quantile(0.25)
    Q3 = filtered_data[column_2].quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - removal_factor * IQR
    upper_bound = Q3 + removal_factor * IQR

    filtered_data = filtered_data[(filtered_data[column_2] >= lower_bound) & (filtered_data[column_2] <= upper_bound)]

    print("Number of outliers removed: ", len(data) - len(filtered_data))

    return filtered_data

def plot_replicate_correlation(file_path, replicate_pair):
    # Read the dataset
    df = pd.read_csv(file_path)

    # Extract the replicates
    replicate_1, replicate_2 = replicate_pair
    data = df[[replicate_1, replicate_2]]

    # Remove outliers
    data = remove_outliers(data, replicate_1, replicate_2)

    # Calculate sample size
    sample_size = len(data)

    # Perform linear regression
    slope, intercept, r_value, p_value, std_err = linregress(data[replicate_1], data[replicate_2])

    # Calculate the line of best fit
    line = slope * data[replicate_1] + intercept

    # Calculate R squared value
    r_squared = r_value ** 2

    # Create a figure with gridspec for custom layout
    fig = plt.figure(figsize=(12, 10))
    gs = fig.add_gridspec(2, 2, width_ratios=(5, 1), height_ratios=(1, 5), hspace=0.05, wspace=0.05)

    # Main scatter plot
    ax_scatter = fig.add_subplot(gs[1, 0])
    ax_scatter.scatter(data[replicate_1], data[replicate_2], alpha=0.5)
    ax_scatter.plot(data[replicate_1], line, color='red')  # line of best fit
    ax_scatter.set_xlabel(f'Normalized Counts ({column_to_sample[replicate_1]})')
    ax_scatter.set_ylabel(f'Normalized Counts ({column_to_sample[replicate_2]})')
    plt.grid()

    # Histogram for first replicate
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax_scatter)
    ax_histx.hist(data[replicate_1], bins=150, density=True, color='blue')
    ax_histx.set_ylabel('Density')

    # Histogram for second replicate (rotated 90 degrees)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax_scatter)
    ax_histy.hist(data[replicate_2], bins=150, density=True, orientation='horizontal', color='blue')
    ax_histy.set_xlabel('Density')

    # Set the title for the main scatter plot
    ax_histx.set_title(f'Scatter Plot of Normalized Counts\n{column_to_sample[replicate_1]} vs {column_to_sample[replicate_2]}\n$R^2 = {r_squared:.3f}, N={sample_size}$')

    # Hide x labels and tick labels for top plots and y ticks for right plots
    plt.setp(ax_histx.get_xticklabels(), visible=False)
    plt.setp(ax_histy.get_yticklabels(), visible=False)

    # Save the plot
    plt.savefig(f"{column_to_sample[replicate_1]}_vs_{column_to_sample[replicate_2]}_scatterplot.png", dpi=300)
    plt.close()

# Generate scatterplots for each pair of replicates
for pair in replicate_pairs:
    plot_replicate_correlation(file_path, pair)
