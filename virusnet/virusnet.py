#!/usr/bin/env python

import os
import subprocess
import multiprocessing as mp
import argparse
import pandas as pd
import warnings
import time
import logging
import queue
import threading
import tempfile

# Suppress h5py UserWarning
warnings.filterwarnings("ignore", category=UserWarning, module="h5py")

# Use a relative path based on the location of this script
R_SCRIPT_PATH = os.path.join(os.path.dirname(__file__), "process_fasta.r")

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

def init_worker():
    """
    Initialize logging for each worker process with a simpler format.
    """
    logging.basicConfig(level=logging.INFO, format='%(message)s')

def run_batch_r_script(gpu_id, fasta_list_file, output_dir, model_binary, model_genus, genus_labels, window_size, step, binary_batch_size, genus_batch_size):
    """
    Run the R script in batch mode for multiple files with a single model load.
    
    Args:
        gpu_id: GPU ID to use (set to None for CPU)
        fasta_list_file: Path to file containing list of FASTA files
        Other args are passed to the R script
    """
    device_type = f"GPU {gpu_id}" if gpu_id is not None else "CPU"
    logging.info(f"Started batch processing on {device_type}")
    
    # Set up the environment for this process
    env = os.environ.copy()
    if gpu_id is not None:
        env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    else:
        env["CUDA_VISIBLE_DEVICES"] = ""
    
    # Command to run the R script in batch mode
    cmd = [
        "Rscript", R_SCRIPT_PATH,
        "--fasta_list", fasta_list_file,
        "--output_dir", output_dir,
        "--model_binary", model_binary,
        "--model_genus", model_genus,
        "--genus_labels", genus_labels,
        "--window_size", str(window_size),
        "--step", str(step),
        "--binary_batch_size", str(binary_batch_size),
        "--genus_batch_size", str(genus_batch_size)
    ]
    
    try:
        subprocess.run(cmd, check=True, env=env)
        logging.info(f"Completed batch processing on {device_type}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in batch processing on {device_type}: {e}")

def download_models(download_path, verify=False):
    """
    Placeholder for downloading model files.
    """
    os.makedirs(download_path, exist_ok=True)
    logging.info(f"Downloading models to {download_path} with verify={verify}")

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
        "--binary_batch_size", 
        type=int, 
        default=1, 
        help="Batch size for binary prediction (default: 1)"
    )
    predict_parser.add_argument(
        "--genus_batch_size", 
        type=int, 
        default=10, 
        help="Batch size for genus prediction (default: 10)"
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

    # Set up logging for the main process
    logging.basicConfig(level=logging.INFO, format='%(message)s')

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
        logging.info(f"Found {len(fasta_files)} FASTA files to process.")

        # Record the start time
        start_time = time.time()

        # Distribute files based on GPU or CPU processing
        if args.num_gpus > 0:
            # Create a list of files for each GPU
            gpu_file_lists = [[] for _ in range(args.num_gpus)]
            
            # Distribute files among GPUs (round-robin)
            for i, file in enumerate(fasta_files):
                gpu_id = i % args.num_gpus
                gpu_file_lists[gpu_id].append(file)
            
            # Create temporary list files for each GPU
            temp_list_files = []
            for gpu_id, file_list in enumerate(gpu_file_lists):
                # Create a temporary file with the list of FASTA files for this GPU
                temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=f'_gpu{gpu_id}.txt')
                for file_path in file_list:
                    temp_file.write(f"{file_path}\n")
                temp_file.close()
                temp_list_files.append(temp_file.name)
                logging.info(f"GPU {gpu_id} will process {len(file_list)} files")
            
            # Process each batch in parallel
            processes = []
            for gpu_id, fasta_list_file in enumerate(temp_list_files):
                process = mp.Process(
                    target=run_batch_r_script,
                    args=(
                        gpu_id, fasta_list_file, args.output_dir, args.model_binary,
                        args.model_genus, args.genus_labels, args.window_size, args.step,
                        args.binary_batch_size, args.genus_batch_size
                    )
                )
                process.start()
                processes.append(process)
            
            # Wait for all processes to complete
            for process in processes:
                process.join()
                
            # Clean up temporary files
            for temp_file in temp_list_files:
                os.unlink(temp_file)
                
            logging.info("All GPU processing completed")
        else:
            # For CPU mode, split files among CPU processes
            cpu_file_lists = [[] for _ in range(args.num_cpu_processes)]
            
            # Distribute files among CPU processes (round-robin)
            for i, file in enumerate(fasta_files):
                process_id = i % args.num_cpu_processes
                cpu_file_lists[process_id].append(file)
            
            # Create temporary list files for each CPU process
            temp_list_files = []
            for process_id, file_list in enumerate(cpu_file_lists):
                # Create a temporary file with the list of FASTA files for this process
                temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=f'_cpu{process_id}.txt')
                for file_path in file_list:
                    temp_file.write(f"{file_path}\n")
                temp_file.close()
                temp_list_files.append(temp_file.name)
                logging.info(f"CPU process {process_id} will process {len(file_list)} files")
            
            # Process each batch in parallel
            processes = []
            for process_id, fasta_list_file in enumerate(temp_list_files):
                process = mp.Process(
                    target=run_batch_r_script,
                    args=(
                        None, fasta_list_file, args.output_dir, args.model_binary,
                        args.model_genus, args.genus_labels, args.window_size, args.step,
                        args.binary_batch_size, args.genus_batch_size
                    )
                )
                process.start()
                processes.append(process)
            
            # Wait for all processes to complete
            for process in processes:
                process.join()
                
            # Clean up temporary files
            for temp_file in temp_list_files:
                os.unlink(temp_file)
                
            logging.info("All CPU processing completed")

        logging.info("All files have been processed")

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
            logging.info(f"Combined summary saved to {combined_path}")
        else:
            logging.warning("No summary files generated.")

        # Calculate and display the total processing time
        end_time = time.time()
        total_time = end_time - start_time
        hours = int(total_time // 3600)
        minutes = int((total_time % 3600) // 60)
        seconds = int(total_time % 60)
        logging.info(f"Total processing time: {hours} hours, {minutes} minutes, and {seconds} seconds.")