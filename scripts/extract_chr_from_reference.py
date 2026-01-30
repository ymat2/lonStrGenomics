import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta")
    parser.add_argument("-g", "--gff")
    parser.add_argument("-c", "--chrom")
    parser.add_argument("-o", "--outdir")
    args = parser.parse_args()

    fa = extract_chrom_from_fasta(args.fasta, args.chrom)
    outpath4fasta = Path(args.outdir, args.chrom + ".fa")
    with open(outpath4fasta, "w") as f:
        for line in fa:
            f.write(line)

    gff = extract_chrom_from_gff(args.gff, args.chrom)
    outpath4gff = Path(args.outdir, args.chrom + ".gff")
    with open(outpath4gff, "w") as f:
        for line in gff:
            f.write(line)


def extract_chrom_from_fasta(fasta: Path, chrom: str) -> list:
    extracted_fasta = []
    contig_name = ""
    with open(fasta) as f:
        for line in f:
            if line[0] == ">":
                contig_name = line.split()[0][1:]
                if contig_name == chrom:
                    extracted_fasta.append(line)
                else:
                    continue
            elif contig_name == chrom:
                extracted_fasta.append(line)
            else:
                continue
    return extracted_fasta


def extract_chrom_from_gff(gff: Path, chrom: str) -> list:
    extracted_gff = []
    with open(gff) as f:
        for line in f:
            if line[0] == "#":
                extracted_gff.append(line)
            elif line.split("\t")[0] == chrom:
                extracted_gff.append(line)
            else:
                continue
    return extracted_gff


if __name__ == "__main__":
    main()
