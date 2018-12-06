#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument("infile", nargs="?", help="input bam file",  type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument("outfile",nargs="?", help="output bam file", type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument('-d', action="store_true", default=False)
args = parser.parse_args()

import glob
import os
import sys
import re
import itertools 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


args.infile = "overlap.by.asm.bed"

def pairNum():
	f = open("asmbeds.txt")
	conv = {}
	for idx, line in enumerate(f):
		if(line[0] == "#"):
			continue
		line = line.strip().split("/")
		#print(line[7])
		conv[str(idx+1)] = line[7].lower()
	return(conv)


df = pd.read_table(args.infile)
df.columns = ["chr", "start", "end", "golga_id", "strand", "genome", "aln_chr", "aln_start", "aln_end", "contig", "contig_start", "contig_end", "contig_len", "aln_ID", "overlap"]
conv = pairNum()

df["golga_id"]=df["golga_id"].astype(str)
df["strand"]=df["strand"].astype(str)

df["genome"] = df["genome"].replace(conv)

df["%_covered"] = df["overlap"]/(df["end"]-df["start"])*100


# require the majority of the region to be covered
df = df[df["%_covered"] > 90.0]

# require the majority of the cotig to be aligned 
df["contig_end"] = df["contig_end"].astype(int)
df["contig_len"] = df["contig_len"].astype(int)
df["contig_start"] = df["contig_start"].astype(int)
df["frac_of_contig_in_aln"] = (df["contig_end"]-df["contig_start"])/df["contig_len"]


#df.sort_values(by=["frac_of_contig_in_aln"], inplace=True)
df['frac_max'] = df.groupby(['contig'])['frac_of_contig_in_aln'].transform(max)
df = df[ df["frac_of_contig_in_aln"] >= df["frac_max"] ]
df.drop(columns = ["frac_max"], inplace=True)

df.sort_values(by=["chr", "start", "end", "genome"], inplace=True)
df.to_csv("tbl.tab", index=False, sep="\t")


print(df)
if(True):
	df2 = df.loc[:, "chr":"genome"]
	df2 = pd.DataFrame(df2.groupby(["chr", "start", "end", "golga_id", "strand"]).aggregate(set))
	df2["genome"] = df2["genome"].apply(sorted).str.join(", ")
	df2.to_csv("simple.tab",  sep="\t")


