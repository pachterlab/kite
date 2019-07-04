import argparse
import csv
import os
from collections import OrderedDict

parser=argparse.ArgumentParser()
parser.add_argument("FeatureRefCSV", help="path to csv file with whitelist Feature Barcodes")
parser.add_argument("--t2g", help="path to output t2g file, default ./FeaturesMismatch.t2g", default="./FeaturesMismatch.t2g", type=str)
parser.add_argument("--fa", help="path to output fasta file, default ./FeaturesMismatch.fa", default="./FeaturesMismatch.fa", type=str)
parser.add_argument("--header", help="include this flag if FeatureRefCSV has a header, omit if not", action="store_true")
parser.add_argument("-q", "--quiet", help="don't print processed Feature Barcodes", action="store_true")
args=parser.parse_args()

def version():
    print("0.0.2")

"""
get_tags requires a csv-formatted file of feature barcode names and sequences. If your CSV file contains a header, use the --header flag. 
"""    
def get_tags(filename):
    with open(filename, mode='r') as csv_file:
        csv_reader = csv.reader(csv_file)
        tags = {}
        if args.header: 
            next(csv_reader) #skip header unless the no_header option is specified
            if not args.quiet: print('\nCSV includes header row\n')
        else: 
            if not args.quiet: 
                if not args.quiet: print('\nCSV does not include header row.\n')
        for row in csv_reader:
            tags[row[0].strip()] = row[1].strip()
    return tags

def make_mismatch_map(FeatureDict):
    odict = OrderedDict()
    counter=0
    for item in FeatureDict:
        name=(item)
        seq=FeatureDict[item]
        if counter == 0:
            feature_barcode_length = len(seq)
            if not args.quiet: print("Feature Barcode Length: "+str(feature_barcode_length)+'\n')
            if not args.quiet: print('Read ' + str(len(FeatureDict)) +' Feature Barcodes:\n')
            counter+=1
        if not args.quiet:
            print(name)
            print(seq)
        odict[name+'-*-*'] = str(seq)[:feature_barcode_length]
        for pos in range(feature_barcode_length):
            letter =str(seq)[pos]
            barcode=list(str(seq)[:feature_barcode_length])
            if letter=='A':
                barcode[pos]='T'
                odict[name+'-'+str(pos)+'-1'] = "".join(barcode)
                barcode[pos]='G'
                odict[name+'-'+str(pos)+'-2'] = "".join(barcode)
                barcode[pos]='C'
                odict[name+'-'+str(pos)+'-3'] = "".join(barcode)
            elif letter=='G':
                barcode[pos]='T'
                odict[name+'-'+str(pos)+'-1'] = "".join(barcode)
                barcode[pos]='A'
                odict[name+'-'+str(pos)+'-2'] = "".join(barcode)
                barcode[pos]='C'
                odict[name+'-'+str(pos)+'-3'] = "".join(barcode)
            elif letter=='C':
                barcode[pos]='T'
                odict[name+'-'+str(pos)+'-1'] = "".join(barcode)
                barcode[pos]='G'
                odict[name+'-'+str(pos)+'-2'] = "".join(barcode)
                barcode[pos]='A'
                odict[name+'-'+str(pos)+'-3'] = "".join(barcode)
            else:
                barcode[pos]='A'
                odict[name+'-'+str(pos)+'-1'] = "".join(barcode)
                barcode[pos]='G'
                odict[name+'-'+str(pos)+'-2'] = "".join(barcode)
                barcode[pos]='C'
                odict[name+'-'+str(pos)+'-3'] = "".join(barcode)
                
    return odict


def write_mismatch_map(tag_map, mismatch_t2g_path, mismatch_fasta_path):
    tagmap_file = open(mismatch_t2g_path, "w+")
    tagmap_fasta = open(mismatch_fasta_path, "w+")
    for i in list(tag_map.keys()):
        if i[-4:]=='-*-*':
            #print(i[:-4]+'\t'+i[:-4]+'\t'+i[:-4])
            tagmap_file.write(i[:-4]+'\t'+i[:-4]+'\t'+i[:-4]+'\n')
            tagmap_fasta.write(">" + i[:-4] + "\n" +tag_map[i] + "\n")
        else:
            #print(i+'\t'+'-'.join(i.split('-')[:-2])+'\t'+'-'.join(i.split('-')[:-2]))
            tagmap_file.write(i+'\t'+'-'.join(i.split('-')[:-2])+'\t'+'-'.join(i.split('-')[:-2])+'\n')
            tagmap_fasta.write(">" + i + "\n" +tag_map[i] + "\n")
    tagmap_file.close()
    tagmap_fasta.close()
    
#wrapper function for make_mismatch_map and write_mismatch_map
def kite_mismatch_maps(FeatureDict, mismatch_t2g_path, mismatch_fasta_path):
    write_mismatch_map(make_mismatch_map(FeatureDict), mismatch_t2g_path, mismatch_fasta_path)
    if not args.quiet: print("\nThe t2g and fasta files are now ready \n")

tags = get_tags(args.FeatureRefCSV)

t2g_path = args.t2g

fasta_path= args.fa

kite_mismatch_maps(tags, t2g_path, fasta_path)

