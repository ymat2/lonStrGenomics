
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    args = parser.parse_args()

    gtf = clean_gtf(args.infile)
    with open(args.outfile, "w") as f:
        for l in gtf:
            f.write("\t".join(l)+"\n")


def clean_gtf(gtf):
    new_lines = []
    with open(gtf) as f:
        for l in f:
            if l[0] == "#":
                new_lines.append([l.rstrip("\n")])
            else:
                l = l.rstrip("\n").split("\t")
                l[-1] = extract_gene_from_gtf(l[-1])
                new_lines.append(l)
    return new_lines


def extract_gene_from_gtf(desc):
    descs = desc.split("; ")
    descs = {i.split(" ")[0]: i.split(" ")[1].strip('"') for i in descs if " " in i}
    gene = descs.get("gene_id", "NA")
    return gene


if __name__ == "__main__":
    main()
