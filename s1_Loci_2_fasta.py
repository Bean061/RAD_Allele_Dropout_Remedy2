import argparse
import os

def Loci2partition(ipyrad_locifile, output_dir):
    os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists

    with open(ipyrad_locifile, 'r') as f1:
        lines = f1.readlines()

    sum_header = "#nexus\nbegin sets;\n"
    locus_sequence = ""
    accession_n = 0
    a = 0

    for text in lines:
        text = text.strip()

        if not text.startswith("/"):
            name = text.split(" ")[0].split(".trimmed")[0]
            sequence = text.split(" ")[-1]
            a = len(sequence)
            locus_sequence += f">{name}\n{sequence}\n"
            accession_n += 1

        else:
            number = text.find("|")
            number2 = text.rfind("|")
            loci_number = text[number + 1:number2]
            output_locus = os.path.join(output_dir, f"locus_{loci_number}.fasta")

            with open(output_locus, 'a') as fw:
                fw.write(locus_sequence)

            locus_sequence = ""
            a = 0
            accession_n = 0

    print(f"Output fasta files saved to: {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert ipyrad .loci file to per-locus fasta files.",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-i", "--inputlocifile", dest='input_file', type=str, required=True,
                        help='Full path to the .loci file from ipyrad.')

    parser.add_argument("-o", "--outputdir", dest='output_dir', type=str, default="R1R2_output",
                        help='Directory to save per-locus fasta files.')

    args = parser.parse_args()

    input_file = os.path.realpath(args.input_file)
    output_dir = os.path.realpath(args.output_dir)

    Loci2partition(input_file, output_dir)


if __name__ == '__main__':
    main()
