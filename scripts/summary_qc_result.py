import argparse
import json
from pathlib import Path


def main():
    args = parse_args()
    json_files = retrieve_json(args.indir)
    _total_length = []
    _total_reads = []
    _lines = []
    for json_file in json_files:
        with open(json_file) as j:
            sample_name = str(json_file).split("/")[-1].rstrip(".json")
            json_elements = json.load(j)
            json_elements = select_json_element(json_elements)
            _total_length.append(json_elements[6]*json_elements[3]/2+json_elements[7]*json_elements[3]/2)
            _total_reads.append(json_elements[3])
            json_elements.insert(0, sample_name)
            _lines.append(list(map(str, json_elements)))
    mean_read_length_across_samples = sum(_total_length)/sum(_total_reads)
    print("mean_read_length_across_samples: ", mean_read_length_across_samples)
    for l in _lines:
        print("\t".join(l))


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir", help="Path to directory contains QC report files")
    args = parser.parse_args()
    return args


def retrieve_json(dir: Path) -> list:
    path = Path(dir)
    return list(path.rglob('*.json'))


def select_json_element(json: dict) -> list:
    b1 = json["summary"]["before_filtering"]["total_reads"]
    b2 = json["summary"]["before_filtering"]["q20_rate"]
    b3 = json["summary"]["before_filtering"]["q30_rate"]
    a1 = json["summary"]["after_filtering"]["total_reads"]
    a2 = json["summary"]["after_filtering"]["q20_rate"]
    a3 = json["summary"]["after_filtering"]["q30_rate"]
    a4 = json["summary"]["after_filtering"]["read1_mean_length"]
    a5 = json["summary"]["after_filtering"]["read2_mean_length"]
    f1 = json["filtering_result"]["passed_filter_reads"]
    f2 = json["filtering_result"]["low_quality_reads"]
    f3 = json["filtering_result"]["too_many_N_reads"]
    f4 = json["filtering_result"]["too_short_reads"]
    f5 = json["filtering_result"]["too_long_reads"]
    d1 = json["duplication"]["rate"]
    return [b1, b2 ,b3 ,a1 ,a2 ,a3 ,a4 ,a5 ,f1 ,f2 ,f3 ,f4 ,f5 ,d1]


if __name__ == "__main__":
    main()
