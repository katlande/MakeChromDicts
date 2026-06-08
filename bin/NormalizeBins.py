import pandas as pd
import numpy as np
import argparse
import sys, getopt

# checks if chrnames are messed up:
def is_numeric_string(val):
    try:
        float(str(val))
        return True
    except ValueError:
        return False

def main(argv):
   parser = argparse.ArgumentParser(description="Check the transcript denisty around a slice.")
   
   # Input
   parser.add_argument("-b", "--bed", required=True, help="multisamsummary output")
   parser.add_argument("-s", "--samples", type=str, required=True, help="sample names from meta")
   parser.add_argument("-o", "--output", required=True, help="path to output file")
   
   args = parser.parse_args()
   inbed = (args.bed)
   samps = args.samples.split(",")
   outFile = args.output
   
   # read file in:
   bedfile = pd.read_csv(inbed, sep="\t",  low_memory=False)
   # rename columns:
   bedfile.columns = ["chr", "start", "end"] + samps
   # replace NAs with 0 counts:
   bedfile = bedfile.fillna(0)
   # normalize:
   bedfile.iloc[:, 3:] = (bedfile.iloc[:, 3:] / bedfile.iloc[:, 3:].sum())*1000000
   
   # get bin centers and remove unnecessary columns:
   binsize = bedfile.iat[0, 2]-bedfile.iat[0, 1]
   bedfile["start"] = bedfile["start"] + int(binsize/2)
   bedfile = bedfile.drop(columns=['end'])
   bedfile = bedfile.rename(columns={"start": "pos"})
   
   # fix chromosome names if necessary
   # add 'chr' character to the beginning of the chrom names if any numeric-only names are present;
   # this will prevent downstream errors
   chrnames = bedfile['chr'].unique()
   
   if any(is_numeric_string(x) for x in chrnames):
       bedfile['chr'] = 'chr' + bedfile['chr']
   
   
   # write output file:
   bedfile.to_csv(outFile, sep="\t", index=False) 
   
if __name__ == "__main__":
   main(sys.argv[1:])