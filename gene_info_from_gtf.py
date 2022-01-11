import pandas as pd
import argparse
from tqdm import tqdm

tqdm.pandas()


def parseinfo(s):
    items = s['additional_info'].replace(" ", "").split(";")
    gene_id = None
    gene_type = None
    gene_name = None
    for i in items:
        if 'gene_id' in i:
            gene_id = i.split('"')[1]
        elif 'gene_type' in i:
            gene_type = i.split('"')[1]
        elif 'gene_name' in i:
            gene_name = i.split('"')[1]

    return gene_id, gene_type, gene_name


def unioninterval(l):    
    def get_union(a, b):
        a0, a1 = a
        b0, b1 = b
        if a1 + 1 == b0:
            return [1, a0, b1]
        if a1 < b0:
            return [2, a, b]
        if a0 <= b0:
            return [1, a0, max(a1, b1)]
    
    l.sort(key=lambda x: x[0])
    i, j = 0, len(l) - 1
    while True:
        if i == j:
            break
        rangea, rangeb = l.pop(i), l.pop(i)
        union_i = get_union(rangea, rangeb)
        if union_i[0] == 1:
            l.insert(i, [union_i[1], union_i[2]])
        elif union_i[0] == 2:
            l.insert(i, union_i[1])
            l.insert(i + 1, union_i[2])
            i += 1
        j = len(l) - 1
    return l


def addup(sdf):
    interval = unioninterval(
        sdf[['genomic_start_location', 'genomic_end_location']].values.tolist()
    )
    length = 0
    for x in interval:
        length += x[1] - x[0] + 1        
    return pd.Series([sdf.iloc[0, :]['gene_name'], length])


def process_gtf(file, output):
    print("Loading GTF file...")
    df = pd.read_csv(
        file, sep="\t", comment="#", compression='gzip', 
        usecols=[0, 1, 2, 3, 4, 8],
        header=None, names=(
            "chromosome_name", "annotation_source", "feature_type",
            "genomic_start_location", "genomic_end_location", 
            # "score", "genomic_strand", "genomic_phase", 
            "additional_info"
            ),
        low_memory=False
    )

    print("GTF file loaded. Parsing additional info...")
    df[['ensembl_id', 'gene_type', 'gene_name']] = df[['additional_info']].progress_apply(parseinfo, axis=1, result_type='expand')

    new_df = df[["chromosome_name", "annotation_source", "feature_type",
            "genomic_start_location", "genomic_end_location", 
            # "score", "genomic_strand", "genomic_phase", 
            "ensembl_id", "gene_type", "gene_name"]]

    print("Additional info parsing completed. Calculating exon length...")
    exon_df = new_df.loc[new_df['feature_type'] == 'exon', ].groupby('ensembl_id').progress_apply(addup)
    exon_df.columns = ['symbol', 'length']

    print("Exon length calculation completed. Saving results.")
    exon_df.to_csv(output)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gtf', '-f', help='gtf file', required=True)
    parser.add_argument('--output', '-o', help='Output filename', default="./gtf.info.csv")
    args = parser.parse_args()

    try:
        process_gtf(args.gtf, args.output)
    except Exception as e:
        print(e)

    