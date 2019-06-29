import argparse
import kite
import csv

parser=argparse.ArgumentParser()
parser.add_argument("ref", help="feature csv file")
args=parser.parse_args()

def get_tags(filename):
    with open(filename, mode='r') as csv_file:
        csv_reader = csv.reader(csv_file)
        tags = {}
        next(csv_reader)
        for row in csv_reader:
            tags[row[1].strip()] = row[4].strip()
    return tags

tags = get_tags(args.ref)

t2g_path = "./Features_t2g.txt"

fasta_path= "./FeaturesMismatch.fa"

kite.kite_mismatch_maps(tags, t2g_path, fasta_path)

