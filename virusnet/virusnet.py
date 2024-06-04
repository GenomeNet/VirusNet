#!/usr/bin/env python

import os
import subprocess
import argparse
import json
import urllib.request
import hashlib

def download_models(download_path):
    with open('models.json', 'r') as file:
        models = json.load(file)
    os.makedirs(download_path, exist_ok=True)
    print(f"Models will be downloaded to: {download_path}")
    print("You can change the download location by using the --path argument.")
    
    for key, model_info in models.items():
        url = model_info['url']
        expected_hash = model_info['hash']
        file_path = os.path.join(download_path, os.path.basename(url))
        
        if os.path.exists(file_path) and verify_file_hash(file_path, expected_hash):
            print(f"{key} already downloaded and verified.")
            user_input = input("File already exists and is verified. Do you want to re-download it? (yes/no): ")
            if user_input.lower() != 'yes':
                continue
        
        print(f"Downloading {key} to {file_path}...")
        urllib.request.urlretrieve(url, file_path, reporthook=download_progress)
        if verify_file_hash(file_path, expected_hash):
            print(f"Successfully verified the hash for {key}.")
        else:
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

def check_files(download_path):
    with open('models.json', 'r') as file:
        models = json.load(file)
    for key, model_info in models.items():
        file_path = os.path.join(download_path, os.path.basename(model_info['url']))
        if not os.path.exists(file_path) or not verify_file_hash(file_path, model_info['hash']):
            return False
    return True

def run_prediction(input, output, model_paths):
    """
    Function to run the R script for virus prediction using the specified arguments.
    """
    r_script_path = os.path.join(os.path.dirname(__file__), "predict.r")
    command = [
        "Rscript", r_script_path, 
        '--input', input, 
        '--output', output, 
        '--model_binary', model_paths['binary_model'], 
        '--model_genus', model_paths['genus_model'], 
        '--labels_genus', model_paths['genus_labels']
    ]
    subprocess.run(command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='VirusNet Tool')
    subparsers = parser.add_subparsers(dest='command')

    download_parser = subparsers.add_parser('download')
    download_parser.add_argument('--path', type=str, default=os.path.expanduser('~/.virusnet'))

    predict_parser = subparsers.add_parser('predict')
    predict_parser.add_argument('--input', type=str, required=True)
    predict_parser.add_argument('--output', type=str, required=True)
    predict_parser.add_argument('--path', type=str, default=os.path.expanduser('~/.virusnet'))

    args = parser.parse_args()

    if args.command == 'download':
        download_models(args.path)
    elif args.command == 'predict':
        if check_files(args.path):
            model_paths = {key: os.path.join(args.path, os.path.basename(info['url'])) for key, info in json.load(open('models.json')).items()}
            run_prediction(args.input, args.output, model_paths)
        else:
            print("Model files are missing or corrupted. Please download them again.")