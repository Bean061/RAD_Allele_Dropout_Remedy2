from Bio import SeqIO
import re
import argparse
import os

def split_to_separate_files(input_dir, output_r1_dir, output_r2_dir):
    if not os.path.exists(output_r1_dir):
        os.makedirs(output_r1_dir)
    if not os.path.exists(output_r2_dir):
        os.makedirs(output_r2_dir)

    fasta_files = [f for f in os.listdir(input_dir) if f.endswith(".fasta") or f.endswith(".fa")]

    for fasta_file in fasta_files:
        input_path = os.path.join(input_dir, fasta_file)
        output_r1_path = os.path.join(output_r1_dir, fasta_file)
        output_r2_path = os.path.join(output_r2_dir, fasta_file)

        with open(output_r1_path, 'w') as r1_handle, open(output_r2_path, 'w') as r2_handle:
            for record in SeqIO.parse(input_path, 'fasta'):
                separator = re.search(r'N{3,}[N\-]*', str(record.seq))
                if separator:
                    pos_start = separator.start()
                    pos_end = separator.end()

                    r1 = record[:pos_start]
                    r1.id = record.id
                    r1.description = ""
                    SeqIO.write(r1, r1_handle, 'fasta')

                    r2 = record[pos_end:]
                    r2.id = record.id
                    r2.description = ""
                    SeqIO.write(r2, r2_handle, 'fasta')
                else:
                    print(f"Warning: No separator found in {record.id} from {fasta_file}")

    print("Splitting complete.")

def main():
    parser = argparse.ArgumentParser(description='Batch split FASTA files by NNNN separator into R1 and R2 directories.')
    parser.add_argument("-i", "--inputdir", dest='input_dir', help='Input directory with FASTA files', required=True)
    parser.add_argument("-o", "--outputR1dir", dest='output_r1_dir', help='Output directory for R1 FASTA files', required=True)
    parser.add_argument("-p", "--outputR2dir", dest='output_r2_dir', help='Output directory for R2 FASTA files', required=True)

    args = parser.parse_args()

    split_to_separate_files(args.input_dir, args.output_r1_dir, args.output_r2_dir)

if __name__ == '__main__':
    main()
