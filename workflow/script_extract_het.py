import os
import pandas as pd

# Get the list of samples to analyze
df_sple = pd.read_csv("mapped/sample_list", header=None, index_col=False)
# Make the list of samples
sample_list = df_sple[0].to_list()

print("Sample List:", sample_list)

# Initialize lists to store sample names and results
samples = []
final_results = []

for f in sample_list:
    file_path = 'heterozygosity/' + f + '.ml'
    print(f"Processing {f} from file: {file_path}")

    if os.path.isfile(file_path) and os.path.getsize(file_path) > 0:
        try:
            # Read the file
            df = pd.read_csv(file_path, sep='\\s+', header=None, index_col=False)

            # Debug print to see the data read from the file
            print("Data from file:", df.head())

            # Check if the DataFrame has at least one row and two columns
            if not df.empty and df.shape[0] > 0 and df.shape[1] > 1:
                # Calculate heterozygosity for each site
                het_by_site = df.iloc[:,1] / (df[0] + df[1] + df[2])

                
                # Ensure no division by zero
                if het_by_site.sum() != 0:
                    # Make an average of heterozygosity on all sites
                    result = het_by_site.mean(axis=0)
                else:
                    result = 0
                    print(f"Warning: Sum of elements in row for {file_path} is zero. Result set to 0.")
            else:
                result = 0
                print(f"Warning: Data in {file_path} is not in the expected format. Result set to 0.")
        except pd.errors.EmptyDataError:
            result = 0
            print(f"Warning: {file_path} is empty or has no columns to parse. Result set to 0.")
        except Exception as e:
            result = 0
            print(f"Error processing {file_path}: {e}. Result set to 0.")
    else:
        result = 0
        print(f"Warning: {file_path} does not exist or is empty. Result set to 0.")

    # Append the sample name and result
    samples.append(f)
    final_results.append(result)

# Create a DataFrame for the results
results_df = pd.DataFrame({
    'Sample': samples,
    'Heterozygosity': final_results
})

print(results_df)
# Save the DataFrame to a CSV file
results_df.to_csv('heterozygosity/heterozygosity_results_angsd.csv', index=False)
