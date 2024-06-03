# cli.py
import subprocess

def main():
    subprocess.run(["Rscript", "virusnet/test.R"])

if __name__ == "__main__":
    main()