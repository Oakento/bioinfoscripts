# bioinfoscripts

Scripts that generate common files in bioinformatics that no one shares.

| file | description | usage | requirements |
| --- | --- | --- | --- |
| `gene_info_from_gtf.ipynb` | extract gene symbol & exon length from gtf file |  | `tqdm>=4.61.2` `pandas` |
|`gene_info_from_gtf.py` | extract gene symbol & exon length from gtf file, script | `python3 gene_info_from_gtf.py -f <gtf gzip file> -o <otuput file>` | `tqdm>=4.61.2` `pandas` |
| `tcga_metadata_convert.py` | convert tcga metadata from json to table | `python3 tcga_metadata_convert.py -i <metadata json file> -o <output file>` | `pandas` |
| `tcga_data_merge.py` | merge tcga data by columns | `python3 tcga_data_merge.py -d <tcga data directory> -m <metadata json file> -o <output file> -p <used core numbers>` | `pandas` |

