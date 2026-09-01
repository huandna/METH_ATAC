import argparse
from Bio import SeqIO
def process_fasta(input_fasta, output_file):
	with open(output_file, 'w') as out:
	    for record in SeqIO.parse(input_fasta, 'fasta'):
	        out.write(f"{record.id}\t0\t{len(record.seq)}\n")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process FASTA file to extract sequence ID and length.")
    parser.add_argument("--input", required=True, help="Input FASTA file path.")
    parser.add_argument("--output", required=True, help="Output text file path.")
    args = parser.parse_args()
    process_fasta(args.input, args.output)