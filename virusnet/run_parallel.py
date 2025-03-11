import os
import subprocess
import multiprocessing
import pandas as pd

def run_r_script(file, gpu_id, r_script_path, output_dir, model_binary, model_genus, genus_labels, window_size, step):
    """
    Run the R script for a single FASTA file on the specified GPU.
    
    Args:
        file (str): Path to the FASTA file.
        gpu_id (int): GPU index to use (e.g., 0 for GPU0, 1 for GPU1).
        r_script_path (str): Path to the R script.
        output_dir (str): Output directory for results.
        model_binary (str): Path to the binary model.
        model_genus (str): Path to the genus model.
        genus_labels (str): Path to the genus labels RDS file.
        window_size (int): Window size for processing.
        step (int): Step size for processing.
    """
    # Set the environment variable to assign the GPU
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    
    # Construct the Rscript command with your existing arguments
    cmd = [
        "Rscript", r_script_path,
        "--fasta_file", file,
        "--output_dir", output_dir,
        "--model_binary", model_binary,
        "--model_genus", model_genus,
        "--genus_labels", genus_labels,
        "--window_size", str(window_size),
        "--step", str(step)
    ]
    
    # Execute the command
    try:
        subprocess.run(cmd, check=True, env=env)
        print(f"Completed processing {file} on GPU {gpu_id}")
    except subprocess.CalledProcessError as e:
        print(f"Error processing {file} on GPU {gpu_id}: {e}")

def main():
    # Configuration (adjust these to match your setup)
    input_folder = "../test"  # Directory with FASTA files
    output_dir = "output_folder_3"               # Directory for output files
    r_script_path = "./process_fasta.r"        # Path to your R script
    model_binary = "~/.virusnet/transfer_learning_virus_bert_5_85.h5"   # Binary model path
    model_genus = "~/.virusnet/virus_genus_2023-01-23.hdf5"   # Genus model path
    genus_labels = "~/.virusnet/genus_labels.rds"  # Genus labels path
    window_size = 1000                         # Window size
    step = 5000                                # Step size
    num_gpus = 2                               # Number of GPUs available (e.g., 2 for GPU0 and GPU1)

    # Create output directory if it doesn’t exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all FASTA files
    fasta_files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith(".fasta")]
    if not fasta_files:
        raise ValueError(f"No FASTA files found in {input_folder}")
    
    print(f"Found {len(fasta_files)} FASTA files to process.")
    
    # Use a process pool with the number of GPUs
    with multiprocessing.Pool(processes=num_gpus) as pool:
        # Assign each FASTA file to a GPU (cycling through GPU0, GPU1, etc.)
        tasks = [
            (file, i % num_gpus, r_script_path, output_dir, model_binary, model_genus, genus_labels, window_size, step)
            for i, file in enumerate(fasta_files)
        ]
        
        # Run tasks in parallel
        pool.starmap(run_r_script, tasks)
    
    # Combine summary CSV files into one (optional)
    summary_files = [os.path.join(output_dir, f"{os.path.splitext(os.path.basename(f))[0]}_summary.csv") 
                     for f in fasta_files]
    summary_dfs = [pd.read_csv(sf, sep=";") for sf in summary_files if os.path.exists(sf)]
    if summary_dfs:
        combined_summary = pd.concat(summary_dfs, ignore_index=True)
        combined_path = os.path.join(output_dir, "combined_summary.csv")
        combined_summary.to_csv(combined_path, sep=";", index=False)
        print(f"Combined summary saved to {combined_path}")
    else:
        print("No summary files generated.")

if __name__ == "__main__":
    main()