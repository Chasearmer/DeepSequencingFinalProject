import pandas as pd

def filter_m6a_genes(file_path):
    # Read the dataset
    data = pd.read_csv(file_path)

    # Filter out rows where 'm6A_ESC_mm' is not 'm6a'
    m6a_data = data[data['m6A_ESC_mm'] == 'm6A']

    # Create new file name
    new_file_name = file_path.rsplit('.', 1)[0] + '_m6a_only.csv'

    # Save the filtered dataset
    m6a_data.to_csv(new_file_name, index=False)
    print(f"Filtered data saved to {new_file_name}")

# Example usage
file_path = 'datasets/de_YTHDF_Tri_KO_vs_WT_m6A.csv'  # Replace with your actual file path
filter_m6a_genes(file_path)
