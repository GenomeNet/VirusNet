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

def run_r_script(file, gpu_id, output_dir, model_binary, model_genus, genus_labels, window_size, step, binary_batch_size, genus_batch_size):
    # Get base filename for cleaner logs
    base_filename = os.path.basename(file)
    pid = os.getpid()
    
    # Set up the environment for this process
    env = os.environ.copy()
    if gpu_id is not None:
        env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
        device = f"GPU {gpu_id}"
    else:
        env["CUDA_VISIBLE_DEVICES"] = ""
        device = "CPU"
    
    # Simplified logging - just the essential information
    logging.info(f"Started processing {base_filename} on {device}")
    
    # Command to run the R script
    cmd = [
        "Rscript", R_SCRIPT_PATH,
        "--fasta_file", file,
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
        logging.info(f"Completed processing {base_filename} on {device}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error processing {base_filename} on {device}: {e}")

def download_models(download_path, verify=False):
    """
    Placeholder for downloading model files.
    """
    os.makedirs(download_path, exist_ok=True)
    logging.info(f"Downloading models to {download_path} with verify={verify}")

def gpu_worker(gpu_id, task_queue, output_dir, model_binary, model_genus, genus_labels, window_size, step, binary_batch_size, genus_batch_size):
    """
    Worker function that processes tasks for a specific GPU.
    Pulls tasks from the queue until it receives None (indicating no more tasks).
    """
    logging.info(f"GPU {gpu_id} worker started")
    while True:
        file = task_queue.get()
        if file is None:  # Sentinel value to indicate end of queue
            task_queue.task_done()
            break
            
        run_r_script(file, gpu_id, output_dir, model_binary, model_genus, genus_labels, 
                     window_size, step, binary_batch_size, genus_batch_size)
        task_queue.task_done()

def cpu_worker_pool(fasta_files, num_cpu_processes, output_dir, model_binary, model_genus, genus_labels, window_size, step, binary_batch_size, genus_batch_size):
    """
    Process files using CPU workers in a pool.
    """
    logging.info(f"Starting processing with {num_cpu_processes} CPU processes")
    tasks = [
        (file, None, output_dir, model_binary, model_genus, genus_labels, 
         window_size, step, binary_batch_size, genus_batch_size)
        for file in fasta_files
    ]
    
    with mp.get_context("spawn").Pool(processes=num_cpu_processes, initializer=init_worker) as pool:
        pool.starmap(run_r_script, tasks)

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

        # Use GPU or CPU processing
        if args.num_gpus > 0:
            # Create a task queue for each GPU
            gpu_queues = [queue.Queue() for _ in range(args.num_gpus)]
            
            # Distribute files among GPU queues (round-robin)
            for i, file in enumerate(fasta_files):
                gpu_id = i % args.num_gpus
                gpu_queues[gpu_id].put(file)
            
            # Add sentinel values to signal the end of tasks
            for q in gpu_queues:
                q.put(None)
            
            # Create and start a worker thread for each GPU
            threads = []
            for gpu_id, task_queue in enumerate(gpu_queues):
                worker_thread = threading.Thread(
                    target=gpu_worker,
                    args=(gpu_id, task_queue, args.output_dir, args.model_binary, 
                          args.model_genus, args.genus_labels, args.window_size, 
                          args.step, args.binary_batch_size, args.genus_batch_size)
                )
                worker_thread.start()
                threads.append(worker_thread)
            
            # Wait for all worker threads to complete
            for thread in threads:
                thread.join()
                
            logging.info("All GPU processing completed")
        else:
            # Use CPU processing with a pool
            cpu_worker_pool(
                fasta_files, args.num_cpu_processes, args.output_dir, args.model_binary,
                args.model_genus, args.genus_labels, args.window_size, args.step,
                args.binary_batch_size, args.genus_batch_size
            )

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