import pandas as pd
import numpy as np
import argparse


def calc(s, l):
    rate = np.log(s) - np.log(l)
    denom = np.log(np.sum(np.exp(rate)))
    return np.exp(rate - denom + np.log(1e6))



def process_count(file, ref, output):
    df = pd.read_csv(file, sep="\t", index_col=0)
    refdf = pd.read_csv(ref, sep="\t", index_col="gene_id")
    refdf = refdf[["gene_name", "exon_length"]]
    if not np.array_equal(df.index, refdf.index):
        print("Gene id is not valid.")
        return
    np.seterr(divide='ignore')
    print("Start calculating...")
    tpmdf = df.apply(calc, axis=0, args=(refdf['exon_length'],))
    print("Saving...")
    tpmdf.to_csv(output, sep="\t")



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--counts', '-c', help='counts file', required=True)
    parser.add_argument('--ref', '-r', help='reference file', required=True)
    parser.add_argument('--output', '-o', help='Output filename', default="./tpm.csv")
    args = parser.parse_args()

    try:
        process_count(args.counts, args.ref, args.output)
    except Exception as e:
        print(e)

    