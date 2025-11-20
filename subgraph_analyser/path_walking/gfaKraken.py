#!/usr/bin/env python3

# Script to run kraken on each segment of a GFA file.
import sys
import subprocess

def gfa2fa(gfa_input):
  """
  Using subprocess, runs gfatools gfa2fa to convert a GFA file segments to a FASTA file.
  """
  fasta_file = gfa_input.split(".")[0] + ".fasta"
  with open(fasta_file, "w") as f:
      subprocess.run(["gfatools", "gfa2fa", gfa_input], stdout=f)

def run_kraken(gfa_input, kraken_db):
    # create a fasta file of the node sequences in the graph
    gfa2fa(gfa_input)

    # define the fasta file name
    fasta_file = gfa_input.split(".")[0] + ".fasta"
    
    #run kraken on the fasta file
    kraken_out=subprocess.check_output(["kraken2", 
                    "--db", kraken_db, 
                    "--threads", "1",
                    fasta_file])
    # process kraken output using decode
    kraken_out = kraken_out.decode("utf-8")

    return kraken_out

if __name__ == '__main__':
    gfa_input = sys.argv[1]
    kraken_db = sys.argv[2]
    kraken_out = run_kraken(gfa_input, kraken_db)
    #write kraken output to file
    with open(gfa_input.split(".")[0] + "_kraken_out.txt", "w") as f:
        f.write(kraken_out)
