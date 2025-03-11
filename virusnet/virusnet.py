#!/usr/bin/env python

import os
import subprocess
import multiprocessing
import argparse
import pandas as pd

# Hardcoded R script path (assumes process_fasta.r is in the same directory as this script)
R_SCRIPT_PATH = "./process_fasta.r"

def check_files(model_binary, model_genus, genus_labels):
    """
    Check if the specified model files exist.
    Raises FileNotFoundError if any are missing.
    """
    missing = []
    for path in (model_binary, model_genus, genus_labels):
        if not os.path.exists(path):
            missing.append(path)
    if missing:
        raise FileNotFoundError(
            f"Missing model files: {', '.join(missing)}. "
            "Please ensure they exist or run 'virusnet download' to fetch them."
        )

def run_r_script(file, gpu_id, output_dir, model_binary, model_genus, genus_labels, window_size, step):
    """
    Run the R script for a single FASTA file on the specified GPU or CPU.

    Args:
        file (str): Path to the FASTA file.
        gpu_id (int or None): GPU index to use (e.g., 0 for GPU0), or None for CPU.
        output_dir (str): Output directory for results.
        model_binary (str): Path to the binary model.
        model_genus (str): Path to the genus model.
        genus_labels (str): Path to the genus labels RDS file.
        window_size (int): Window size for processing.
        step (int): Step size for processing.
    """
    # Copy the current environment
    env = os.environ.copy()

    # Set CUDA_VISIBLE_DEVICES based on gpu_id
    if gpu_id is not None:
        env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)  # Use specified GPU
    else:
        env["CUDA_VISIBLE_DEVICES"] = ""  # Hide all GPUs to force CPU usage

    # Command to run the R script
    cmd = [
        "Rscript", "./process_fasta.r",
        "--fasta_file", file,
        "--output_dir", output_dir,
        "--model_binary", model_binary,
        "--model_genus", model_genus,
        "--genus_labels", genus_labels,
        "--window_size", str(window_size),
        "--step", str(step)
    ]
    
    try:
        subprocess.run(cmd, check=True, env=env)
        print(f"Completed processing {file}" + (f" on GPU {gpu_id}" if gpu_id is not None else " on CPU"))
    except subprocess.CalledProcessError as e:
        print(f"Error processing {file}" + (f" on GPU {gpu_id}" if gpu_id is not None else " on CPU") + f": {e}")
        
def download_models(download_path, verify=False):
    """
    Placeholder for downloading model files.
    In a full implementation, this would download models to the specified path and verify hashes if requested.
    """
    # For this example, it’s a placeholder. Replace with actual download logic if available.
    os.makedirs(download_path, exist_ok=True)
    print(f"Downloading models to {download_path} with verify={verify}")
    # Example: Add actual urllib.request.urlretrieve calls and hash verification here

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="VirusNet Tool")
    subparsers = parser.add_subparsers(dest="command", required=True, help="Available commands")

    # Download subcommand
    download_parser = subparsers.add_parser("download", help="Download model files")
    download_parser.add_argument(
        "--path", 
        type=str, 
        default=os.path.expanduser("~/.virusnet"), 
        help="Path to save model files"
    )
    download_parser.add_argument(
        "--verify", 
        action="store_true", 
        help="Verify downloaded files with hash"
    )

    # Predict subcommand
    predict_parser = subparsers.add_parser("predict", help="Process FASTA files for viral classification")
    predict_parser.add_argument(
        "--input_folder", 
        type=str, 
        required=True, 
        help="Input folder containing FASTA files"
    )
    predict_parser.add_argument(
        "--output_dir", 
        type=str, 
        required=True, 
        help="Output directory for results"
    )
    predict_parser.add_argument(
        "--model_binary", 
        type=str, 
        default=os.path.expanduser("~/.virusnet/transfer_learning_virus_bert_5_85.h5"), 
        help="Path to the binary model"
    )
    predict_parser.add_argument(
        "--model_genus", 
        type=str, 
        default=os.path.expanduser("~/.virusnet/virus_genus_2023-01-23.hdf5"), 
        help="Path to the genus model"
    )
    predict_parser.add_argument(
        "--genus_labels", 
        type=str, 
        default=os.path.expanduser("~/.virusnet/genus_labels.rds"), 
        help="Path to the genus labels RDS file"
    )
    predict_parser.add_argument(
        "--window_size", 
        type=int, 
        default=1000, 
        help="Window size for processing"
    )
    predict_parser.add_argument(
        "--step", 
        type=int, 
        default=10000, 
        help="Step size for processing"
    )
    predict_parser.add_argument(
        "--num_gpus", 
        type=int, 
        default=0, 
        help="Number of GPUs to use (0 for CPU)"
    )
    predict_parser.add_argument(
        "--num_cpu_processes", 
        type=int, 
        default=4, 
        help="Number of CPU processes if num_gpus=0"
    )

    args = parser.parse_args()

    if args.command == "download":
        download_models(args.path, args.verify)

    elif args.command == "predict":
        # Verify model files exist
        check_files(args.model_binary, args.model_genus, args.genus_labels)

        # Ensure output directory exists
        os.makedirs(args.output_dir, exist_ok=True)

        # Collect FASTA files from input folder
        fasta_files = [
            os.path.join(args.input_folder, f) 
            for f in os.listdir(args.input_folder) 
            if f.endswith(".fasta")
        ]
        if not fasta_files:
            raise ValueError(f"No FASTA files found in {args.input_folder}")
        print(f"Found {len(fasta_files)} FASTA files to process.")

        # Set up parallel processing based on GPU or CPU
        num_processes = args.num_gpus if args.num_gpus > 0 else args.num_cpu_processes
        use_gpu = args.num_gpus > 0

        with multiprocessing.Pool(processes=num_processes) as pool:
            if use_gpu:
                tasks = [
                    (file, i % args.num_gpus, args.output_dir, args.model_binary, 
                     args.model_genus, args.genus_labels, args.window_size, args.step)
                    for i, file in enumerate(fasta_files)
                ]
            else:
                tasks = [
                    (file, None, args.output_dir, args.model_binary, 
                     args.model_genus, args.genus_labels, args.window_size, args.step)
                    for file in fasta_files
                ]
            pool.starmap(run_r_script, tasks)

        # Combine summary CSV files
        summary_files = [
            os.path.join(args.output_dir, f"{os.path.splitext(os.path.basename(f))[0]}_summary.csv")
            for f in fasta_files
        ]
        summary_dfs = [pd.read_csv(sf, sep=";") for sf in summary_files if os.path.exists(sf)]
        if summary_dfs:
            combined_summary = pd.concat(summary_dfs, ignore_index=True)
            combined_path = os.path.join(args.output_dir, "combined_summary.csv")
            combined_summary.to_csv(combined_path, sep=";", index=False)
            print(f"Combined summary saved to {combined_path}")
        else:
            print("No summary files generated.")