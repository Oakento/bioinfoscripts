import pandas as pd
import json
import argparse
import multiprocessing
import os


def load_metadata(metadata_path) -> pd.DataFrame:
    with open(metadata_path, 'r') as f:
        metadata = json.load(f)

    return pd.DataFrame(
        [(
            item['file_id'],
            item['file_name'],
            item['associated_entities'][0]['entity_submitter_id'],
            item['associated_entities'][0]['case_id']
        ) for item in metadata],
        columns=['id', 'filename', 'TCGA_ID', 'case_id']
    )


def read_file(data_dir, file_dir_id, filename, TCGA_ID) -> pd.DataFrame:
    filepath = os.path.join(data_dir, file_dir_id, filename)
    return pd.read_csv(
        filepath,
        compression='gzip',
        sep='\t', 
        index_col=0, 
        names=['gene_id', TCGA_ID],
        # skipfooter=5, 
        # engine='python'
    )


def start_task(data_dir, sub_metadata_df: pd.DataFrame) -> pd.DataFrame:
    result_df = pd.DataFrame()
    temp = []
    for row in sub_metadata_df.iterrows():
        temp.append(read_file(data_dir, row[1]['id'], row[1]['filename'], row[1]['TCGA_ID']))
        if len(temp) > len(sub_metadata_df) / 3:
            result_df = pd.concat([result_df] + temp, axis=1)
            temp = []
    result_df = pd.concat([result_df] + temp, axis=1)

    return result_df


def mergefile(data_dir, metadata_path, cpu_num=0) -> pd.DataFrame:
    metadata_df = load_metadata(metadata_path)

    if cpu_num == 0:
        cpu_num = multiprocessing.cpu_count()

    slice_length = metadata_df.shape[0] // cpu_num

    task_list = [
        (data_dir, metadata_df[i * slice_length: (i + 1) * slice_length]) for i in range(cpu_num - 1)
    ] + [
        (data_dir, metadata_df[(cpu_num - 1) * slice_length:])
    ]

    result = []
    with multiprocessing.Pool(cpu_num) as p:
        result += p.starmap(start_task, task_list)

    df = pd.concat(result, axis=1)
    return df.drop(df.tail(5).index)


def main(data_dir, metadata_path, output, cpu_num):
    df = mergefile(data_dir, metadata_path, int(cpu_num))
    df.to_csv(output, sep="\t")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', '-d', help='TCGA data directory', required=True)
    parser.add_argument('--metadata', '-m', help='Metadata json file', required=True)
    parser.add_argument('--output', '-o', help='Output filename', default="./tcga.merge.tsv")
    parser.add_argument('--parallel', '-p', help="Used core numbers", default=multiprocessing.cpu_count())
    args = parser.parse_args()

    try:
        main(args.dir, args.metadata, args.output, args.parallel)
    except Exception as e:
        print(e)
