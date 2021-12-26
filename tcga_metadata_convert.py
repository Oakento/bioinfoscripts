import pandas as pd
import argparse
import json


def convert(file) -> pd.DataFrame:
    with open(file, 'r') as metadataf:
        metadata = json.load(metadataf)

    df = pd.DataFrame(
        [(
            item['file_id'],
            item['file_name'],
            item['associated_entities'][0]['entity_submitter_id'],
            item['associated_entities'][0]['case_id']
        ) for item in metadata],
        columns=['id', 'filename', 'TCGA_ID', 'case_id']
    )

    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', help='Metadata json file', required=True)
    parser.add_argument('--output', '-o', help='Output file', default='./tcga.metadata.csv')
    args = parser.parse_args()

    df = convert(args.input)
    df.to_csv(args.output)
