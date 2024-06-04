#!/usr/bin/env python

import os
import subprocess
import argparse

def run_prediction(model='genus', input='test.fasta', output='prediction.csv', step=1000, batch_size=32):
    """
    Function to run the R script for virus prediction using the specified arguments.
    """

    # Define the path to your R script
    r_script_path = os.path.join(os.path.dirname(__file__), "predict.r")

    # Define the command that you would use to run the R script from the command line
    command = ["Rscript", r_script_path, 
               '--input', input,
               '--output', output,
               '--step', str(step),
               '--batch_size', str(batch_size)]

    # Run the command
    subprocess.run(command)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run virus predictions.')
    parser.add_argument('-i', '--input', type=str, default='test.fasta',
                        help='Input fasta file.')
    parser.add_argument('-o', '--output', type=str, default='prediction.csv',
                        help='Output CSV file.')
    parser.add_argument('-s', '--step', type=int, default=1000,
                        help='Step size to iterate though sequences.')
    parser.add_argument('-b', '--batch_size', type=int, default=32,
                        help='Number of samples processed in one batch.')

    args = parser.parse_args()

    run_prediction(input=args.input, output=args.output, step=args.step, batch_size=args.batch_size)