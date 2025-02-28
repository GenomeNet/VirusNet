#!/usr/bin/env python

import os
import subprocess
import argparse
import json
import urllib.request
import hashlib
import sys
import pandas as pd
import tempfile as tempfile_module
import shutil
from Bio import SeqIO

# List of known bacteriophage genera
BACTERIOPHAGE_GENERA = [
    "Cystovirus", "Inovirus", "Nonagvirus", "Ceduovirus", "Biseptimavirus", 
    "Moineauvirus", "Felixounavirus", "Claudivirus", "Pbunavirus", "Phikmvvirus", 
    "Tequatrovirus", "Lambdavirus", "Brussowvirus", "Kayfunavirus", "Bruynoghevirus", 
    "Lederbergvirus", "Peduovirus", "Kuttervirus", "Kayvirus", "Dhillonvirus", 
    "Friunavirus", "Drulisvirus", "Casadabanvirus", "Benedictvirus", "Skunavirus", 
    "Jiaodavirus", "Phietavirus", "Pakpunavirus", "Aureusvirus", "Fletchervirus"
]

def download_models(download_path, verify=False):
    models_path = os.path.join(sys.prefix, 'bin', 'models.json')
    print(f"Looking for models.json at: {models_path}")
    if not os.path.exists(models_path):
        raise FileNotFoundError(f"models.json not found at {models_path}")
    
    with open(models_path, 'r') as file:
        models = json.load(file)
    os.makedirs(download_path, exist_ok=True)
    print(f"Models will be downloaded to: {download_path}")
    print("You can change the download location by using the --path argument.")
    
    for key, model_info in models.items():
        url = model_info['url']
        expected_hash = model_info['hash']
        file_path = os.path.join(download_path, os.path.basename(url))
        
        if os.path.exists(file_path):
            if verify and not verify_file_hash(file_path, expected_hash):
                print(f"Hash mismatch for {key}, will re-download.")
            else:
                print(f"{key} already downloaded" + (" and verified." if verify else "."))
                user_input = input(f"File already exists. Do you want to re-download it? (yes/no): ")
                if user_input.lower() != 'yes':
                    continue
        
        print(f"Downloading {key} to {file_path}...")
        urllib.request.urlretrieve(url, file_path, reporthook=download_progress)
        if verify:
            if verify_file_hash(file_path, expected_hash):
                print(f"Successfully verified the hash for {key}.")
            else:
                print(f"Warning: Hash mismatch for {key}, download might be corrupted.")
                user_input = input("Continue anyway? (yes/no): ")
                if user_input.lower() != 'yes':
                    raise ValueError(f"Hash mismatch for {key}, download might be corrupted.")
        os.environ[f"VIRUSNET_{key.upper()}"] = file_path

def download_progress(block_num, block_size, total_size):
    downloaded = block_num * block_size
    if total_size > 0:
        progress_percentage = min(downloaded * 100 / total_size, 100)  # Ensure percentage does not exceed 100%
        progress_bar = f"[{'=' * int(progress_percentage // 2)}{' ' * (50 - int(progress_percentage // 2))}]"
        print(f"\rDownloading: {progress_bar} {progress_percentage:.2f}%", end='')
        if downloaded >= total_size:
            print()

def verify_file_hash(file_path, expected_hash):
    sha256 = hashlib.sha256()
    with open(file_path, 'rb') as f:
        for chunk in iter(lambda: f.read(4096), b""):
            sha256.update(chunk)
    calculated_hash = sha256.hexdigest()
    return calculated_hash == expected_hash

def check_files(download_path, verify=False):
    models_path = os.path.join(sys.prefix, 'bin', 'models.json')
    if not os.path.exists(models_path):
        raise FileNotFoundError(f"models.json not found at {models_path}")
    
    with open(models_path, 'r') as file:
        models = json.load(file)
    for key, model_info in models.items():
        file_path = os.path.join(download_path, os.path.basename(model_info['url']))
        if not os.path.exists(file_path):
            return False
        if verify and not verify_file_hash(file_path, model_info['hash']):
            print(f"Warning: Hash mismatch for {key}, file might be corrupted.")
            return False
    return True

def generate_simplified_genus_output(genus_summary_path, output_path):
    """
    Generate a simplified genus output file with just the contig name, most likely genus, its probability,
    and a flag indicating if it's a bacteriophage.
    """
    try:
        # Read the genus summary CSV
        try:
            genus_df = pd.read_csv(genus_summary_path)
        except Exception:
            # Try with different CSV reading options if the default fails
            genus_df = pd.read_csv(genus_summary_path, sep=None, engine='python')
        
        # Create a new simplified dataframe
        simplified_df = pd.DataFrame()
        simplified_df['contig_name'] = genus_df['contig_name']
        
        # Determine which column to use for the most likely genus
        if 'clean_taxon' in genus_df.columns:
            simplified_df['most_likely_genus'] = genus_df['clean_taxon']
        elif 'max_virus' in genus_df.columns:
            # Strip 'mean_' prefix if present
            simplified_df['most_likely_genus'] = genus_df['max_virus'].apply(
                lambda x: x[5:] if isinstance(x, str) and x.startswith('mean_') else x
            )
        else:
            # Calculate the most likely genus from mean columns
            raise ValueError("Could not find clean_taxon or max_virus columns in genus summary")
        
        # Add the probability column
        if 'max_probability' in genus_df.columns:
            simplified_df['probability'] = genus_df['max_probability']
        else:
            # Try to find the maximum probability among mean columns
            mean_cols = [col for col in genus_df.columns if col.startswith('mean_')]
            if mean_cols:
                simplified_df['probability'] = genus_df.apply(lambda row: max([row[col] for col in mean_cols]), axis=1)
            else:
                raise ValueError("Could not find probability information in genus summary")
        
        # Ensure probability is formatted as a percentage
        simplified_df['probability'] = simplified_df['probability'].apply(lambda x: round(float(x) * 100, 2))
        
        # Add the bacteriophage flag
        simplified_df['is_bacteriophage'] = simplified_df['most_likely_genus'].apply(
            lambda x: 'TRUE' if any(phage.lower() == str(x).lower() for phage in BACTERIOPHAGE_GENERA) else 'FALSE'
        )
        
        # Save to CSV
        simplified_df.to_csv(output_path, index=False)
        print(f"Simplified genus output saved to {output_path}")
        return True
    except Exception as e:
        print(f"Error generating simplified genus output: {str(e)}")
        return False

def generate_combined_output(binary_results_path, genus_summary_path, output_path, virus_threshold=0.5, include_std_dev=True):
    """
    Generate a combined output file that includes all contigs with their binary classification
    and genus information (for viral contigs).
    
    Parameters:
    -----------
    binary_results_path : str
        Path to the binary classification results CSV file
    genus_summary_path : str
        Path to the genus summary CSV file
    output_path : str
        Path to save the combined output file
    virus_threshold : float
        Threshold for classifying a contig as viral
    include_std_dev : bool
        Whether to include the standard deviation column in the output
    """
    try:
        # Read the binary classification results
        binary_df = pd.read_csv(binary_results_path)
        
        # Check if we have a summarized version with standard deviation
        binary_summarized_path = binary_results_path.replace("binary_results.csv", "binary_results_summarized.csv")
        if os.path.exists(binary_summarized_path) and include_std_dev:
            try:
                binary_summarized_df = pd.read_csv(binary_summarized_path)
                
                # Check if this is from metagenome mode (simplified format) or regular mode
                if 'sd_is_virus' in binary_summarized_df.columns:
                    # Regular mode with multiple predictions per contig
                    # Create a mapping of contig names to their standard deviation
                    sd_mapping = {}
                    for _, row in binary_summarized_df.iterrows():
                        sd_mapping[row['contig_name']] = row.get('sd_is_virus', 0)
                    
                    # Add standard deviation column to the main dataframe
                    binary_df['virus_std_dev'] = binary_df['contig_name'].apply(
                        lambda x: round(sd_mapping.get(x, 0) * 100, 2) if x in sd_mapping else 0
                    )
                else:
                    # Metagenome mode with one prediction per contig
                    print("Metagenome mode detected - using simplified format with no standard deviation")
                    # Add a default std_dev column of zeros
                    binary_df['virus_std_dev'] = 0
            except Exception as e:
                print(f"Warning: Could not process binary summarized results: {str(e)}")
                # Add a default std_dev column
                binary_df['virus_std_dev'] = 0
        elif include_std_dev:
            # If no summarized file exists but std_dev is requested, add a default std_dev column
            binary_df['virus_std_dev'] = 0
        
        # Initialize genus-related columns
        binary_df['most_likely_genus'] = 'NA'
        binary_df['genus_probability'] = 0.0
        binary_df['is_bacteriophage'] = 'FALSE'
        
        # Read genus summary if it exists
        if os.path.exists(genus_summary_path):
            try:
                genus_df = pd.read_csv(genus_summary_path)
                
                # Create a mapping of contig names to genus information
                genus_info = {}
                
                # Determine which columns to use
                genus_col = 'clean_taxon' if 'clean_taxon' in genus_df.columns else 'max_virus'
                prob_col = 'max_probability' if 'max_probability' in genus_df.columns else None
                
                if genus_col and prob_col:
                    for _, row in genus_df.iterrows():
                        contig_name = row['contig_name']
                        # Handle possible different naming conventions
                        contig_id = contig_name.split()[0]
                        
                        # Get the genus name, removing 'mean_' prefix if present
                        genus = row[genus_col]
                        if isinstance(genus, str) and genus.startswith('mean_'):
                            genus = genus[5:]
                            
                        # Store both the full name and the ID for matching
                        genus_info[contig_name] = {
                            'genus': genus,
                            'probability': row[prob_col]
                        }
                        genus_info[contig_id] = {
                            'genus': genus,
                            'probability': row[prob_col]
                        }
                    
                    # Update the binary dataframe with genus information for viral contigs
                    for idx, row in binary_df.iterrows():
                        if row['virus'] >= virus_threshold:
                            contig_name = row['contig_name']
                            contig_id = contig_name.split()[0]
                            
                            if contig_name in genus_info:
                                binary_df.at[idx, 'most_likely_genus'] = genus_info[contig_name]['genus']
                                binary_df.at[idx, 'genus_probability'] = round(genus_info[contig_name]['probability'] * 100, 2)
                            elif contig_id in genus_info:
                                binary_df.at[idx, 'most_likely_genus'] = genus_info[contig_id]['genus']
                                binary_df.at[idx, 'genus_probability'] = round(genus_info[contig_id]['probability'] * 100, 2)
                            else:
                                binary_df.at[idx, 'most_likely_genus'] = 'Unknown_virus'
                
                # Add the bacteriophage flag
                binary_df['is_bacteriophage'] = binary_df['most_likely_genus'].apply(
                    lambda x: 'TRUE' if any(phage.lower() == str(x).lower() for phage in BACTERIOPHAGE_GENERA) else 'FALSE'
                )
                
            except Exception as e:
                print(f"Warning: Could not process genus summary: {str(e)}")
                # Mark viral contigs with unknown genus
                for idx, row in binary_df.iterrows():
                    if row['virus'] >= virus_threshold:
                        binary_df.at[idx, 'most_likely_genus'] = 'Unknown_virus'
        
        # Format probabilities as percentages
        binary_df['virus'] = binary_df['virus'].apply(lambda x: round(float(x) * 100, 2))
        if 'non_virus' in binary_df.columns:
            binary_df['non_virus'] = binary_df['non_virus'].apply(lambda x: round(float(x) * 100, 2))
        
        # Rename columns for clarity
        binary_df = binary_df.rename(columns={
            'virus': 'virus_probability',
            'non_virus': 'non_virus_probability'
        })
        
        # Add a clear virus classification column
        binary_df['is_virus'] = binary_df['virus_probability'].apply(
            lambda x: 'TRUE' if x >= virus_threshold * 100 else 'FALSE'
        )
        
        # Save to CSV
        binary_df.to_csv(output_path, index=False)
        print(f"Combined output saved to {output_path}")
        return True
    except Exception as e:
        print(f"Error generating combined output: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

def run_prediction(input, output, model_paths, step_size=1000, batch_size=100, mode='binary', metagenome=False, virus_threshold=0.5, genus_threshold=0.5):
    """
    Function to run the R script for virus prediction using the specified arguments.
    """
    if mode == 'binary':
        if metagenome:
            r_script_path = os.path.join(os.path.dirname(__file__), "predict_binary_metagenome.r")
        else:
            r_script_path = os.path.join(os.path.dirname(__file__), "predict_binary.r")
    elif mode == 'genus':
        r_script_path = os.path.join(os.path.dirname(__file__), "predict_genus.r")
    else:
        raise ValueError(f"Invalid mode: {mode}")

    command = [
        "Rscript", r_script_path,
        '--input', input,
        '--output', output,
        '--model_binary', model_paths['binary_model'], 
        '--model_genus', model_paths['genus_model'], 
        '--labels_genus', model_paths['genus_labels'],
        '--step_size', str(step_size),
        '--batch_size', str(batch_size),
        '--virus_threshold', str(virus_threshold),
        '--genus_threshold', str(genus_threshold)
    ]
    subprocess.run(command)
    
    # Generate the simplified genus output if in genus mode
    if mode == 'genus':
        genus_summary_path = os.path.join(output, "genus_summary_output.csv")
        if os.path.exists(genus_summary_path):
            simplified_output_path = os.path.join(output, "simplified_genus_output.csv")
            generate_simplified_genus_output(genus_summary_path, simplified_output_path)

def run_auto_mode(input, output, model_paths, step_size=1000, batch_size=100, virus_threshold=0.5, genus_threshold=0.5):
    """
    Function to run the auto mode workflow:
    1. Run binary metagenome prediction to identify viral contigs
    2. Run genus prediction only on viral contigs
    3. Merge results into a combined output with viral taxonomy
    """
    print("Running auto mode...")
    
    # Create temporary directories for intermediate outputs
    temp_dir = tempfile_module.mkdtemp()
    binary_output_dir = os.path.join(temp_dir, "binary_output")
    genus_output_dir = os.path.join(temp_dir, "genus_output")
    os.makedirs(binary_output_dir, exist_ok=True)
    os.makedirs(genus_output_dir, exist_ok=True)
    
    try:
        # First run metagenome mode to classify viral vs non-viral
        print("Step 1: Running binary metagenome classification...")
        run_prediction(input, binary_output_dir, model_paths, step_size, batch_size, 
                      mode='binary', metagenome=True, virus_threshold=virus_threshold, genus_threshold=genus_threshold)
        
        viral_fasta_path = os.path.join(binary_output_dir, "viral_contigs.fasta")
        binary_results_path = os.path.join(binary_output_dir, "binary_results.csv")
        
        # Check if binary results file exists
        if not os.path.exists(binary_results_path):
            print("Error: Binary classification results not found. Auto mode could not be completed.")
            return
            
        # Try to read the binary results
        try:
            # Directly examine the CSV file content to debug the format
            print("Debug - CSV file content (first few lines):")
            with open(binary_results_path, 'r') as f:
                for i, line in enumerate(f):
                    if i < 5:  # Print first 5 lines
                        print(line.strip())
            
            # Read with delimiter explicitly specified and quoting enabled
            binary_df = pd.read_csv(binary_results_path, low_memory=False, delimiter=',', quotechar='"')
            
            print(f"Debug - Binary results head: {binary_df.head()}")
            print(f"Debug - Binary results shape: {binary_df.shape}")
            print(f"Debug - Binary results columns: {binary_df.columns}")
            
            if binary_df.empty:
                print("Error: Binary classification results are empty. Auto mode could not be completed.")
                return
            
            # Check if virus column exists and has values
            if 'virus' not in binary_df.columns:
                print("Error: 'virus' column not found in binary results.")
                return
                
            # Convert columns to numeric, handling missing or non-numeric values
            print(f"Debug - Before conversion virus column values: {binary_df['virus'].head().values}")
            
            # Use a direct numeric conversion approach
            binary_df['virus'] = pd.to_numeric(binary_df['virus'], errors='coerce').fillna(0)
            binary_df['non_virus'] = pd.to_numeric(binary_df['non_virus'], errors='coerce').fillna(0)
            
            print(f"Debug - After conversion virus column values: {binary_df['virus'].head().values}")
            
        except Exception as e:
            print(f"Error reading binary classification results: {str(e)}")
            import traceback
            traceback.print_exc()
            return
        
        # Check if there are any viral contigs based on the virus column
        viral_contigs_count = (binary_df['virus'] >= virus_threshold).sum()
        print(f"Debug - Found {viral_contigs_count} viral contigs based on virus column")
        
        # Check if viral contigs were found in the CSV
        if viral_contigs_count == 0:
            print("No viral contigs found in classification results. Skipping genus prediction.")
            # Copy binary results to final output and add NA for taxon
            binary_df['most_probable_taxon'] = 'NA'
            os.makedirs(output, exist_ok=True)
            binary_df.to_csv(os.path.join(output, "auto_results.csv"), index=False)
            # Also copy binary_results.csv to output
            shutil.copy(binary_results_path, os.path.join(output, "binary_results.csv"))
            return
        
        # Check if viral FASTA file exists and is not empty
        fasta_missing = not os.path.exists(viral_fasta_path) or os.path.getsize(viral_fasta_path) == 0
        if fasta_missing:
            print("Warning: Viral contigs identified but FASTA file is missing or empty.")
            print("Skipping genus prediction and marking viral contigs as 'Unknown_virus'.")
            
            # Skip genus prediction and mark viral contigs
            binary_df['most_probable_taxon'] = binary_df.apply(
                lambda row: 'Unknown_virus' if row['virus'] >= virus_threshold else 'NA', axis=1)
            os.makedirs(output, exist_ok=True)
            binary_df.to_csv(os.path.join(output, "auto_results.csv"), index=False)
            shutil.copy(binary_results_path, os.path.join(output, "binary_results.csv"))
            return
        
        # Run genus prediction on viral contigs
        print("Step 2: Running genus classification on viral contigs...")
        try:
            run_prediction(viral_fasta_path, genus_output_dir, model_paths, step_size, batch_size, 
                       mode='genus', metagenome=False, virus_threshold=virus_threshold, genus_threshold=genus_threshold)
        except Exception as e:
            print(f"Warning: Error during genus prediction: {str(e)}")
            print("Will attempt to use any available genus results.")
        
        # Merge the results
        print("Step 3: Merging results...")
        
        # Initialize most_probable_taxon column with 'NA'
        binary_df['most_probable_taxon'] = 'NA'
        
        # Check if genus summary output exists
        genus_summary_path = os.path.join(genus_output_dir, "genus_summary_output.csv")
        
        if os.path.exists(genus_summary_path):
            try:
                # Directly examine the genus summary CSV file content
                print("Debug - Genus CSV file content (first few lines):")
                with open(genus_summary_path, 'r') as f:
                    for i, line in enumerate(f):
                        if i < 5:  # Print first 5 lines
                            print(line.strip())
            
                # Try to read the CSV file with different options if needed
                try:
                    genus_df = pd.read_csv(genus_summary_path)
                except Exception as csv_error:
                    print(f"Warning: Error reading genus CSV with default settings: {str(csv_error)}")
                    print("Trying with different CSV reading options...")
                    try:
                        genus_df = pd.read_csv(genus_summary_path, sep=None, engine='python')
                    except Exception as e:
                        print(f"Error: Could not read genus CSV file: {str(e)}")
                        raise
                
                print(f"Debug - Genus summary head: {genus_df.head()}")
                print(f"Debug - Genus summary columns: {genus_df.columns}")
                
                # Look for the clean_taxon column first, then fall back to max_virus
                if 'clean_taxon' in genus_df.columns:
                    print("Using clean_taxon column for mapping")
                    # Fix the contig name format issue by extracting the actual contig ID
                    taxon_mapping = {}
                    for _, row in genus_df.iterrows():
                        contig_name = row['contig_name']
                        # Extract the contig ID from the full name - get the first part before any spaces or semicolons
                        # This should be the actual contig ID like NZ_CAQWZH010000007.1
                        contig_id = contig_name.split()[0]
                        
                        # Print some debug info for the first few entries
                        if len(taxon_mapping) < 5:
                            print(f"Debug - Mapping contig: {contig_name} -> ID: {contig_id} -> Taxon: {row['clean_taxon']}")
                        
                        # Store mappings for both the ID and full name
                        taxon_mapping[contig_id] = row['clean_taxon']
                        taxon_mapping[contig_name] = row['clean_taxon']
                    
                    print(f"Debug - Fixed taxon mapping (first 10 entries): {dict(list(taxon_mapping.items())[:10])}")
                elif 'max_virus' in genus_df.columns:
                    print("Using max_virus column for mapping")
                    taxon_mapping = {}
                    for _, row in genus_df.iterrows():
                        contig_name = row['contig_name']
                        # Extract the contig ID from the full name
                        contig_id = contig_name.split()[0]
                        
                        taxon = row['max_virus']
                        if isinstance(taxon, str) and taxon.startswith('mean_'):
                            taxon = taxon[5:]  # Remove 'mean_' prefix
                        
                        # Store mappings for both the ID and full name
                        taxon_mapping[contig_id] = taxon
                        taxon_mapping[contig_name] = taxon
                else:
                    print("Using column with highest mean value for mapping")
                    taxon_mapping = {}
                    for _, row in genus_df.iterrows():
                        contig_name = row['contig_name']
                        # Extract the contig ID from the full name
                        contig_id = contig_name.split()[0]
                        
                        taxa_columns = [col for col in genus_df.columns if col.startswith('mean_')]
                        if taxa_columns:
                            max_col = max(taxa_columns, key=lambda x: row[x])
                            taxon = max_col[5:]  # Remove 'mean_' prefix
                            
                            # Store mappings for both the ID and full name
                            taxon_mapping[contig_id] = taxon
                            taxon_mapping[contig_name] = taxon
                
                print(f"Debug - Taxon mapping: {taxon_mapping}")
                
                # Create list of viral contigs
                viral_contigs = binary_df[binary_df['virus'] >= virus_threshold]['contig_name'].tolist()
                print(f"Debug - Number of viral contigs identified: {len(viral_contigs)}")
                print(f"Debug - First few viral contigs: {viral_contigs[:5]}")
                
                # Update the binary_df with taxon information for viral contigs
                updated_count = 0
                for idx, row in binary_df.iterrows():
                    if row['virus'] >= virus_threshold:
                        # Try different matching approaches for contig names
                        matched = False
                        contig_name = row['contig_name']
                        
                        # Extract the contig ID - get the first part before any spaces
                        # This is the actual contig ID like NZ_CAQWZH010000007.1
                        contig_id = contig_name.split()[0]
                        
                        # Debug the first few entries
                        if updated_count < 5:
                            print(f"Debug - Looking for match for: {contig_name} -> ID: {contig_id}")
                        
                        # Try matching with just the contig ID
                        if contig_id in taxon_mapping:
                            binary_df.at[idx, 'most_probable_taxon'] = taxon_mapping[contig_id]
                            if updated_count < 5:
                                print(f"Debug - Found match by ID: {contig_id} -> {taxon_mapping[contig_id]}")
                            updated_count += 1
                            matched = True
                        # Direct match with full name
                        elif contig_name in taxon_mapping:
                            binary_df.at[idx, 'most_probable_taxon'] = taxon_mapping[contig_name]
                            if updated_count < 5:
                                print(f"Debug - Found match by full name: {contig_name} -> {taxon_mapping[contig_name]}")
                            updated_count += 1
                            matched = True
                        else:
                            # Try partial matching - but be more careful to avoid incorrect matches
                            for genus_contig in taxon_mapping.keys():
                                # Only match if the genus contig ID is the same as the binary contig ID
                                # or if the full genus contig name is contained in the binary contig name
                                genus_contig_id = genus_contig.split()[0] if ' ' in genus_contig else genus_contig
                                
                                if (contig_id == genus_contig_id or 
                                    (len(genus_contig) > 10 and genus_contig in contig_name)):
                                    binary_df.at[idx, 'most_probable_taxon'] = taxon_mapping[genus_contig]
                                    if updated_count < 5:
                                        print(f"Debug - Found partial match: {contig_id} with {genus_contig_id} -> {taxon_mapping[genus_contig]}")
                                    updated_count += 1
                                    matched = True
                                    break
                        
                        # If no match found, mark as Unknown_virus
                        if not matched:
                            binary_df.at[idx, 'most_probable_taxon'] = 'Unknown_virus'
                            if updated_count < 5:
                                print(f"Debug - No match found for: {contig_id}, marking as Unknown_virus")
                
                print(f"Debug - Updated {updated_count} rows with taxon information")
                
            except Exception as e:
                print(f"Warning: Could not process genus summary results: {str(e)}")
                import traceback
                traceback.print_exc()
                print("Continuing with binary results only.")
                
                # Even if there was an error in genus classification, mark viral contigs
                for idx, row in binary_df.iterrows():
                    if row['virus'] >= virus_threshold:
                        binary_df.at[idx, 'most_probable_taxon'] = 'Unknown_virus'
        else:
            print("Warning: Genus summary results not found. Continuing with binary results only.")
            # Even without genus results, mark viral contigs
            for idx, row in binary_df.iterrows():
                if row['virus'] >= virus_threshold:
                    binary_df.at[idx, 'most_probable_taxon'] = 'Unknown_virus'
        
        # Create the output directory if it doesn't exist
        os.makedirs(output, exist_ok=True)
        
        # Print the final dataframe before saving
        print(f"Debug - Final dataframe structure: {binary_df.shape}")
        print(f"Debug - Final dataframe columns: {binary_df.columns}")
        print(f"Debug - Number of rows with non-NA taxon: {(binary_df['most_probable_taxon'] != 'NA').sum()}")
        
        # Save the combined results
        binary_df.to_csv(os.path.join(output, "auto_results.csv"), index=False)
        
        # Also save a copy of the original binary results
        shutil.copy(binary_results_path, os.path.join(output, "binary_results.csv"))
        
        # Also copy the summarized binary results if they exist
        binary_summarized_path = os.path.join(binary_output_dir, "binary_results_summarized.csv")
        if os.path.exists(binary_summarized_path):
            shutil.copy(binary_summarized_path, os.path.join(output, "binary_results_summarized.csv"))
        
        # Generate a simplified genus output from genus summary results
        simplified_output_path = os.path.join(output, "simplified_genus_output.csv")
        if os.path.exists(genus_summary_path):
            generate_simplified_genus_output(genus_summary_path, simplified_output_path)
        
        # Also create a simplified output that only includes viral contigs with their genus
        viral_contigs_df = binary_df[binary_df['virus'] >= virus_threshold].copy()
        if not viral_contigs_df.empty:
            viral_contigs_df = viral_contigs_df[['contig_name', 'virus', 'most_probable_taxon']]
            viral_contigs_df.columns = ['contig_name', 'virus_probability', 'most_likely_genus']
            viral_contigs_df['virus_probability'] = viral_contigs_df['virus_probability'].apply(lambda x: round(float(x) * 100, 2))
            
            # Add bacteriophage flag
            viral_contigs_df['is_bacteriophage'] = viral_contigs_df['most_likely_genus'].apply(
                lambda x: 'TRUE' if any(phage.lower() == str(x).lower() for phage in BACTERIOPHAGE_GENERA) else 'FALSE'
            )
            
            viral_contigs_df.to_csv(os.path.join(output, "viral_contigs_genus.csv"), index=False)
        
        # Generate standardized combined output
        combined_output_path = os.path.join(output, "combined_results.csv")
        generate_combined_output(binary_results_path, genus_summary_path, combined_output_path, virus_threshold, include_std_dev=False)
        
        # Copy other relevant files
        if os.path.exists(viral_fasta_path):
            viral_output_path = os.path.join(output, "viral_contigs.fasta")
            shutil.copy(viral_fasta_path, viral_output_path)
            
            # Update the viral contigs FASTA with genus information and virus probability
            binary_summarized_path = os.path.join(binary_output_dir, "binary_results_summarized.csv")
            if os.path.exists(genus_summary_path) and os.path.exists(binary_summarized_path):
                update_fasta_with_genus(
                    viral_output_path, 
                    genus_summary_path, 
                    binary_results_path=binary_summarized_path,
                    output_path=viral_output_path
                )
        
        non_viral_fasta_path = os.path.join(binary_output_dir, "non_viral_contigs.fasta")
        if os.path.exists(non_viral_fasta_path):
            shutil.copy(non_viral_fasta_path, os.path.join(output, "non_viral_contigs.fasta"))
        
        print(f"Auto mode completed. Results saved to {output}")
    
    finally:
        # Clean up temporary directories
        shutil.rmtree(temp_dir, ignore_errors=True)

def update_fasta_with_genus(fasta_path, genus_summary_path, binary_results_path=None, output_path=None):
    """
    Update FASTA headers with genus information from genus summary and virus probability.
    
    Parameters:
    -----------
    fasta_path : str
        Path to the FASTA file to update
    genus_summary_path : str
        Path to the genus summary CSV file
    binary_results_path : str, optional
        Path to the binary results CSV file for virus probability
    output_path : str, optional
        Path to save the updated FASTA file. If None, overwrites the input file.
    """
    if not os.path.exists(fasta_path) or not os.path.exists(genus_summary_path):
        print(f"Warning: Could not update FASTA with genus information. Missing files.")
        return False
    
    try:
        # Read the genus summary CSV
        try:
            genus_df = pd.read_csv(genus_summary_path)
        except Exception:
            # Try with different CSV reading options if the default fails
            genus_df = pd.read_csv(genus_summary_path, sep=None, engine='python')
        
        # Create a mapping of contig names to genus information
        genus_info = {}
        
        # Determine which columns to use
        genus_col = 'clean_taxon' if 'clean_taxon' in genus_df.columns else 'max_virus'
        
        if genus_col:
            for _, row in genus_df.iterrows():
                contig_name = row['contig_name']
                # Handle possible different naming conventions
                contig_id = contig_name.split()[0]
                
                # Get the genus name, removing 'mean_' prefix if present
                genus = row[genus_col]
                if isinstance(genus, str) and genus.startswith('mean_'):
                    genus = genus[5:]
                    
                # Store both the full name and the ID for matching
                genus_info[contig_name] = genus
                genus_info[contig_id] = genus
        
        # Read virus probability information if available
        virus_prob_info = {}
        if binary_results_path and os.path.exists(binary_results_path):
            try:
                binary_df = pd.read_csv(binary_results_path)
                
                # Check if we have the summarized results with mean_is_virus
                if 'mean_is_virus' in binary_df.columns:
                    for _, row in binary_df.iterrows():
                        contig_name = row['contig_name']
                        contig_id = contig_name.split()[0]
                        virus_prob = row['mean_is_virus']
                        virus_prob_info[contig_name] = virus_prob
                        virus_prob_info[contig_id] = virus_prob
            except Exception as e:
                print(f"Warning: Could not read virus probability information: {str(e)}")
        
        # If no output path specified, overwrite the input file
        if output_path is None:
            output_path = fasta_path
        
        # Read the FASTA file and update headers
        records = []
        for record in SeqIO.parse(fasta_path, "fasta"):
            # Get the original ID and description
            original_id = record.id
            original_desc = record.description
            
            # Check if we have genus information for this contig
            genus = None
            if original_id in genus_info:
                genus = genus_info[original_id]
            elif original_desc in genus_info:
                genus = genus_info[original_desc]
            else:
                # Try to match by the first part of the ID (before any spaces)
                contig_id = original_id.split()[0]
                if contig_id in genus_info:
                    genus = genus_info[contig_id]
            
            # Check if we have virus probability information
            virus_prob = None
            if original_id in virus_prob_info:
                virus_prob = virus_prob_info[original_id]
            elif original_desc in virus_prob_info:
                virus_prob = virus_prob_info[original_desc]
            else:
                contig_id = original_id.split()[0]
                if contig_id in virus_prob_info:
                    virus_prob = virus_prob_info[contig_id]
            
            # Build the new description
            new_desc = original_desc
            
            # Add virus probability if available
            if virus_prob is not None:
                if " | " in new_desc:
                    new_desc = f"{new_desc} | virus_prop={virus_prob:.4f}"
                else:
                    new_desc = f"{new_desc} | virus_prop={virus_prob:.4f}"
            
            # Add genus information if found
            if genus:
                if " | " in new_desc:
                    new_desc = f"{new_desc} | predicted_genus={genus}"
                else:
                    new_desc = f"{new_desc} | predicted_genus={genus}"
            
            record.description = new_desc
            records.append(record)
        
        # Write the updated FASTA file
        with open(output_path, "w") as output_handle:
            SeqIO.write(records, output_handle, "fasta")
        
        print(f"Updated FASTA headers with additional information: {output_path}")
        return True
    except Exception as e:
        print(f"Error updating FASTA headers: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='VirusNet Tool')
    subparsers = parser.add_subparsers(dest='command')

    download_parser = subparsers.add_parser('download')
    download_parser.add_argument('--path', type=str, default=os.path.expanduser('~/.virusnet'))
    download_parser.add_argument('--verify', action='store_true', help='Enable hash verification of downloaded models')

    predict_parser = subparsers.add_parser('predict')
    predict_parser.add_argument('--input', type=str, required=True)
    predict_parser.add_argument('--output', type=str, required=True)
    predict_parser.add_argument('--path', type=str, default=os.path.expanduser('~/.virusnet'))
    predict_parser.add_argument('--step_size', type=int, default=1000, 
                               help='Step size for prediction (Note: ignored in metagenome mode which uses one prediction per contig)')
    predict_parser.add_argument('--batch_size', type=int, default=100, help='Batch size for prediction')
    predict_parser.add_argument('--mode', type=str, choices=['binary', 'genus', 'both'], default='binary', 
                               help='Prediction mode: binary (viral vs. non-viral), genus (taxonomic classification), or both (run binary first, then genus on viral contigs)')
    predict_parser.add_argument('--metagenome', action='store_true', help='Enable metagenome mode (only applicable for binary mode)')
    predict_parser.add_argument('--auto', action='store_true', 
                              help='Enable auto mode: First runs binary metagenome classification to identify viral contigs, '
                                   'then performs genus prediction only on those viral contigs. '
                                   'The final output includes the binary classification results with an additional column '
                                   'showing the most probable virus taxon for viral contigs (NA for non-viral contigs).')
    predict_parser.add_argument('--verify', action='store_true', help='Enable hash verification of model files')
    predict_parser.add_argument('--virus_threshold', type=float, default=0.5, 
                              help='Threshold for classifying a contig as viral (default: 0.5)')
    predict_parser.add_argument('--genus_threshold', type=float, default=0.5,
                              help='Threshold for assigning a genus classification (default: 0.5)')
    
    args = parser.parse_args()

    if args.command == 'download':
        download_models(args.path, args.verify)
    elif args.command == 'predict':
        if check_files(args.path, args.verify):
            model_paths = {key: os.path.join(args.path, os.path.basename(info['url'])) for key, info in json.load(open(os.path.join(sys.prefix, 'bin', 'models.json'))).items()}
            
            # Check for mode conflicts
            if args.auto and args.mode != 'binary':
                print("Warning: When using --auto mode, the --mode argument is ignored.")
                print("Auto mode will run binary metagenome prediction followed by genus prediction on viral contigs.")
                
            if args.mode == 'both':
                # Run binary first, then genus on viral contigs
                print("Running in 'both' mode: binary classification followed by genus classification on viral contigs...")
                
                # Create temporary directory for binary output
                temp_dir = tempfile_module.mkdtemp()
                binary_output_dir = os.path.join(temp_dir, "binary_output")
                os.makedirs(binary_output_dir, exist_ok=True)
                os.makedirs(args.output, exist_ok=True)
                
                try:
                    # Step 1: Run binary classification in metagenome mode
                    print("Step 1: Running binary metagenome classification...")
                    run_prediction(args.input, binary_output_dir, model_paths, args.step_size, args.batch_size, 
                                 'binary', True, args.virus_threshold, args.genus_threshold)
                    
                    # Check for the viral contigs file
                    viral_fasta_path = os.path.join(binary_output_dir, "viral_contigs.fasta")
                    if not os.path.exists(viral_fasta_path) or os.path.getsize(viral_fasta_path) == 0:
                        print("No viral contigs found or viral contigs file is empty. Skipping genus classification.")
                        
                        # Copy binary results to the output directory
                        binary_results_path = os.path.join(binary_output_dir, "binary_results.csv")
                        if os.path.exists(binary_results_path):
                            shutil.copy(binary_results_path, os.path.join(args.output, "binary_results.csv"))
                            
                        # Also copy the summarized binary results if they exist
                        binary_summarized_path = os.path.join(binary_output_dir, "binary_results_summarized.csv")
                        if os.path.exists(binary_summarized_path):
                            shutil.copy(binary_summarized_path, os.path.join(args.output, "binary_results_summarized.csv"))
                            
                        # Copy viral and non-viral fasta files if they exist
                        for file_name in ["viral_contigs.fasta", "non_viral_contigs.fasta"]:
                            source_path = os.path.join(binary_output_dir, file_name)
                            if os.path.exists(source_path):
                                dest_path = os.path.join(args.output, file_name)
                                shutil.copy(source_path, dest_path)
                                
                                # Update viral contigs FASTA with genus information and virus probability
                                if file_name == "viral_contigs.fasta":
                                    genus_summary_path = os.path.join(args.output, "genus_summary_output.csv")
                                    if os.path.exists(genus_summary_path) and os.path.exists(binary_summarized_path):
                                        update_fasta_with_genus(
                                            dest_path,
                                            genus_summary_path,
                                            binary_results_path=binary_summarized_path,
                                            output_path=dest_path
                                        )
                    else:
                        # Step 2: Run genus classification on viral contigs
                        print("Step 2: Running genus classification on viral contigs...")
                        run_prediction(viral_fasta_path, args.output, model_paths, args.step_size, args.batch_size, 
                                     'genus', False, args.virus_threshold, args.genus_threshold)
                        
                        # Copy binary results to the output directory
                        binary_results_path = os.path.join(binary_output_dir, "binary_results.csv")
                        if os.path.exists(binary_results_path):
                            shutil.copy(binary_results_path, os.path.join(args.output, "binary_results.csv"))
                            
                        # Also copy the summarized binary results if they exist
                        binary_summarized_path = os.path.join(binary_output_dir, "binary_results_summarized.csv")
                        if os.path.exists(binary_summarized_path):
                            shutil.copy(binary_summarized_path, os.path.join(args.output, "binary_results_summarized.csv"))
                            
                        # Copy viral and non-viral fasta files if they exist
                        for file_name in ["viral_contigs.fasta", "non_viral_contigs.fasta"]:
                            source_path = os.path.join(binary_output_dir, file_name)
                            if os.path.exists(source_path):
                                dest_path = os.path.join(args.output, file_name)
                                shutil.copy(source_path, dest_path)
                                
                                # Update viral contigs FASTA with genus information and virus probability
                                if file_name == "viral_contigs.fasta":
                                    genus_summary_path = os.path.join(args.output, "genus_summary_output.csv")
                                    if os.path.exists(genus_summary_path) and os.path.exists(binary_summarized_path):
                                        update_fasta_with_genus(
                                            dest_path,
                                            genus_summary_path,
                                            binary_results_path=binary_summarized_path,
                                            output_path=dest_path
                                        )
                        
                        # Generate combined output with both binary and genus results
                        genus_summary_path = os.path.join(args.output, "genus_summary_output.csv")
                        if os.path.exists(binary_results_path):
                            combined_output_path = os.path.join(args.output, "combined_results.csv")
                            print("Generating combined output with binary and genus classification results...")
                            generate_combined_output(
                                binary_results_path, 
                                genus_summary_path, 
                                combined_output_path,
                                args.virus_threshold,
                                include_std_dev=False  # Exclude standard deviation in 'both' mode
                            )
                    
                    print(f"'Both' mode completed. Results saved to {args.output}")
                    
                finally:
                    # Clean up temporary directories
                    shutil.rmtree(temp_dir, ignore_errors=True)
                    
            # Run in auto mode
            elif args.auto:
                run_auto_mode(args.input, args.output, model_paths, args.step_size, args.batch_size, args.virus_threshold, args.genus_threshold)
            else:
                # Run in regular mode
                run_prediction(args.input, args.output, model_paths, args.step_size, args.batch_size, args.mode, args.metagenome, args.virus_threshold, args.genus_threshold)
        else:
            print("Model files are missing or corrupted. Please download them again.")