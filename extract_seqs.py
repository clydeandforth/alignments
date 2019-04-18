#!/home/b.bss81c/miniconda2/bin/python
#James Doonan
#arguments extract_seqs.py infile.gbk outfile.fasta


from Bio import SeqIO, Seq, AlignIO
import sys
import os

product = 'DNA gyrase subunit B'


with open(sys.argv[2], 'w') as t:
	for rec in SeqIO.parse(sys.argv[1], "genbank"):
       		if rec.features:
        		for feature in rec.features:
                      		if feature.type == "CDS":
                    			if product in feature.qualifiers['product'][0]:
              					t.write(
						">%s|%s %s \n%s\n" % (
                                    	 	feature.qualifiers['locus_tag'][0],
					 	feature.qualifiers['product'][0],						  
				 		rec.description,
					 		#rec.name,
                      				feature.location.extract(rec).seq))

input_file = 'outfile.fasta'
records = (rec.upper() for rec in SeqIO.parse(input_file, 'fasta'))
records = list(records) # make a copy, otherwise our generator
                        # is exhausted after calculating maxlen

maxlen = max(len(record.seq) for record in records)

# pad sequences so that they all have the same length
for record in records:
    if len(record.seq) != maxlen:
        sequence = str(record.seq).ljust(maxlen, '.')
        record.seq = Seq.Seq(sequence)
assert all(len(record.seq) == maxlen for record in records)
output_file = '{}_padded.fasta'.format(os.path.splitext(input_file)[0])
with open(output_file, 'w') as f:
    SeqIO.write(records, f, 'fasta')
#count = SeqIO.write(records, "upper.fas", "fasta")
#print("Converted %i records to upper case" % count)




