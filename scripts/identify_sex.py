import argparse
import glob
import pandas as pd


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--indir", help="Directory contains bam files")
  parser.add_argument("-c", "--sex_chrom", help="Sex chromosome")
  parser.add_argument("-o", "--output", help="Output file path")
  args = parser.parse_args()

  samples = glob.glob(args.indir+"/*")
  with open(args.output, "w") as f:

    _header = "sample\tmeandp\tmeandp_sex_chrom\n"
    f.write(_header)

    for sample in samples:
      sample = sample.split("/")[-1]
      cov_file = args.indir+"/"+sample+"/"+sample+".cov"
      print(cov_file)
      stats = summarize_stat(cov_file, args.sex_chrom)
      f.write(sample+"\t"+"\t".join(stats)+"\n")
      #print(cov_file)

  print("Finish identifying sex.")


def summarize_stat(path, chrom):

  df = pd.read_csv(path, sep = "\t")
  df = df[df["numreads"] != 0]

  total_depth = (df["covbases"]*df["meandepth"]).sum()
  total_coverage = (df["covbases"]*df["coverage"]).sum()
  total_covbase = df["covbases"].sum()
  mean_depth = total_depth/total_covbase

  mean_depth_sex = df[df["#rname"] == chrom]["meandepth"].tolist()[0]

  stats = [mean_depth, mean_depth_sex]
  stats = list(map(str, stats))

  return stats


if __name__ == "__main__":
  main()
