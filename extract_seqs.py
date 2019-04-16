#!/home/b.bss81c/miniconda2/bin/python
#James Doonan
#arguments extract_seqs.py infile.gbk outfile.fasta


from Bio import SeqIO
import sys

product = 'DNA gyrase subunit B'

with open(sys.argv[2], 'w') as nfh:
        for rec in SeqIO.parse(sys.argv[1], "genbank"):
                if rec.features:
                        for feature in rec.features:
                                if feature.type == "CDS":
                                    if product in feature.qualifiers['product'][0]:
                                        nfh.write(">%s|%s %s  \n%s\n" % (
                                                  feature.qualifiers['gene'][0],
                                                  #feature.qualifiers['locus_tag'][0],
						  feature.qualifiers['product'][0],						  
						  rec.description,
						  #rec.name,
                                                  feature.location.extract(rec).seq))
