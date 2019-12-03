#!/scratch/b.bss81c/anaconda3/bin/python
#James Doonan
#arguments complete_script.py infile.gbk outfile.fasta 'gene of interest'


from Bio import SeqIO, Seq, AlignIO
import sys
import os
import re
from Bio.Phylo.TreeConstruction import *
from Bio.Phylo.Consensus import bootstrap_trees, get_support
from Bio.Phylo import draw
from Bio import Phylo, AlignIO
from Bio.Phylo.Consensus import *
from Bio.Align.Applications import ClustalwCommandline
import subprocess
import matplotlib
from Bio.Alphabet import Gapped, generic_dna
from Bio.Alphabet import IUPAC

#import matplotlib.pyplot as plt
import networkx
#import pylab

regex = re.compile(r"(_\d{5})")

gbk_filename = "/scratch/b.bss81c/Prokka/concat_file.gbk"
faa_filename = "outfile.fasta"
input_handle = open(gbk_filename, "r")
string_out = "string.fasta"
string_handle = open(string_out, "w+")
output_handle = open(faa_filename, "w+")
output_string= open(string_out, "w+")

gbk_filename_2 = "/scratch/b.bss81c/Prokka/concat_file.gbk"
faa_filename_2 = "outfile_2.fasta"
input_handle_2 = open(gbk_filename_2, "r")
output_handle_2 = open(faa_filename_2, "w+")
string_out_2 = "string2.fasta"
string_handle_2 = open(string_out_2, "w+")

gbk_filename_3 = "/scratch/b.bss81c/Prokka/concat_file.gbk"
faa_filename_3 = "outfile_3.fasta"
input_handle_3 = open(gbk_filename_3, "r")
output_handle_3 = open(faa_filename_3, "w+")
string_out_3 = "string3.fasta"
string_handle_3 = open(string_out_3, "w+")

gbk_filename_4 = "/scratch/b.bss81c/Prokka/concat_file.gbk"
faa_filename_4 = "outfile_4.fasta"
input_handle_4 = open(gbk_filename_4, "r")
output_handle_4 = open(faa_filename_4, "w+")
string_out_4 = "string4.fasta"
string_handle_4 = open(string_out_4, "w+")

gbk_filename_5 = "/scratch/b.bss81c/Prokka/concat_file.gbk"
faa_filename_5 = "outfile_5.fasta"
input_handle_5 = open(gbk_filename_5, "r")
output_handle_5 = open(faa_filename_5, "w+")
string_out_5 = "string5.fasta"
string_handle_5 = open(string_out_5, "w+")


gbk_filename_6 = "/scratch/b.bss81c/Prokka/concat_file.gbk"
faa_filename_6 = "outfile_6.fasta"
input_handle_6 = open(gbk_filename_6, "r")
output_handle_6 = open(faa_filename_6, "w+")
string_out_6 = "string6.fasta"
string_handle_6 = open(string_out_6, "w+")

def gene_find(gene, sub):
        for rec in SeqIO.parse(input_handle, "genbank"):
                        if rec.features:
                                for feature in rec.features:
                                        if feature.type == "CDS":

                                                if gene in feature.qualifiers['product'][0]:
                                                        if sub in feature.qualifiers['gene'][0]:
                                                                output_handle.write(
                                                                ">%s%s %s \n%s\n" % (
                                                                feature.qualifiers['locus_tag'][0],
                                                                feature.qualifiers['gene'][0],
                                                                rec.description,
                                                                feature.location.extract(rec).seq))
        input_handle.close()
        output_handle.close()

        def change_string():
            for rec in SeqIO.parse('outfile.fasta', 'fasta'):
                rec.description = re.sub(regex, "_", rec.description)
                string_handle.write(">" + rec.description+"\n" + str(rec.seq)+"\n")
            string_handle.close()
        change_string()

        def make_equal():

                    output_string = (string_out)
                    records = (rec.upper() for rec in SeqIO.parse(string_out, 'fasta'))
                    records = list(records) # make a copy, otherwise our generator
                    maxlen = max(len(record.seq) for record in records)
                    for record in records:
                        if len(record.seq) != maxlen:
                            sequence = str(record.seq).ljust(maxlen, '.')
                            record.seq = Seq.Seq(sequence)
                    assert all(len(record.seq) == maxlen for record in records)
                    output_file = '{}_padded.fasta'.format(os.path.splitext(string_out)[0])
                    with open(output_file, 'w') as f:
                        SeqIO.write(records, f, 'fasta')
        make_equal()

        def align_seq():

                cline = ClustalwCommandline("clustalw2", infile="string.fasta")
                stdout = cline()
        align_seq()



def gene_find_2(gene_2, sub_2):
        for rec in SeqIO.parse(input_handle_2, "genbank"):
                        if rec.features:
                                for feature in rec.features:
                                        if feature.type == "CDS":

                                                if gene_2 in feature.qualifiers['product'][0]:
                                                        if sub_2 in feature.qualifiers['gene'][0]:
                                                                output_handle_2.write(
                                                                ">%s%s %s \n%s\n" % (
                                                                feature.qualifiers['locus_tag'][0],
                                                                feature.qualifiers['gene'][0],
                                                                rec.description,
                                                                feature.location.extract(rec).seq))
        input_handle_2.close()
        output_handle_2.close()

        def change_string_2():
            for rec in SeqIO.parse('outfile_2.fasta', 'fasta'):
                rec.description = re.sub(regex, "_", rec.description)
                string_handle_2.write(">" + rec.description+"\n" + str(rec.seq)+"\n")
            string_handle_2.close()
        change_string_2()

        def make_equal_2():

                    output_string_2 = (string_out_2)
                    records = (rec.upper() for rec in SeqIO.parse(string_out_2, 'fasta'))
                    records = list(records) # make a copy, otherwise our generator
                    maxlen = max(len(record.seq) for record in records)
                    for record in records:
                        if len(record.seq) != maxlen:
                            sequence = str(record.seq).ljust(maxlen, '.')
                            record.seq = Seq.Seq(sequence)
                    assert all(len(record.seq) == maxlen for record in records)
                    output_file = '{}_2_padded.fasta'.format(os.path.splitext(string_out_2)[0])
                    with open(output_file, 'w') as f:
                        SeqIO.write(records, f, 'fasta')
        make_equal_2()

        def align_seq_2():

                cline = ClustalwCommandline("clustalw2", infile="string2.fasta")
                stdout = cline()
        align_seq_2()



def gene_find_3(gene_3, sub_3):
        for rec in SeqIO.parse(input_handle_3, "genbank"):
                        if rec.features:
                                for feature in rec.features:
                                        if feature.type == "CDS":

                                                if gene_3 in feature.qualifiers['product'][0]:
                                                        if sub_3 in feature.qualifiers['gene'][0]:
                                                                output_handle_3.write(
                                                                ">%s%s %s \n%s\n" % (
                                                                feature.qualifiers['locus_tag'][0],
                                                                feature.qualifiers['gene'][0],
                                                                rec.description,
                                                                feature.location.extract(rec).seq))
        input_handle_3.close()
        output_handle_3.close()

        def change_string_3():
            for rec in SeqIO.parse('outfile_3.fasta', 'fasta'):
                rec.description = re.sub(regex, "_", rec.description)
                string_handle_3.write(">" + rec.description+"\n" + str(rec.seq)+"\n")
            string_handle_3.close()
        change_string_3()






        def make_equal_3():

                    output_string_3 = (string_out_3)
                    records = (rec.upper() for rec in SeqIO.parse(string_out_3, 'fasta'))
                    records = list(records) # make a copy, otherwise our generator
                    maxlen = max(len(record.seq) for record in records)
                    for record in records:
                        if len(record.seq) != maxlen:
                            sequence = str(record.seq).ljust(maxlen, '.')
                            record.seq = Seq.Seq(sequence)
                    assert all(len(record.seq) == maxlen for record in records)
                    output_file = '{}_3_padded.fasta'.format(os.path.splitext(string_out_3)[0])
                    with open(output_file, 'w') as f:
                        SeqIO.write(records, f, 'fasta')
        make_equal_3()

        def align_seq_3():

                cline = ClustalwCommandline("clustalw2", infile="string3.fasta")
                stdout = cline()
        align_seq_3()


def gene_find_4(gene_4, sub_4):
        for rec in SeqIO.parse(input_handle_4, "genbank"):
                        if rec.features:
                                for feature in rec.features:
                                        if feature.type == "CDS":

                                                if gene_4 in feature.qualifiers['product'][0]:
                                                        if sub_4 in feature.qualifiers['gene'][0]:
                                                                output_handle_4.write(
                                                                ">%s%s %s \n%s\n" % (
                                                                feature.qualifiers['locus_tag'][0],
                                                                feature.qualifiers['gene'][0],
                                                                rec.description,
                                                                feature.location.extract(rec).seq))
        input_handle_4.close()
        output_handle_4.close()

        def change_string_4():
            for rec in SeqIO.parse('outfile_4.fasta', 'fasta'):
                rec.description = re.sub(regex, "_", rec.description)
                string_handle_4.write(">" + rec.description+"\n" + str(rec.seq)+"\n")
            string_handle_4.close()
        change_string_4()



        def make_equal_4():

                    output_string_4 = (string_out_4)
                    records = (rec.upper() for rec in SeqIO.parse(string_out_4, 'fasta'))
                    records = list(records) # make a copy, otherwise our generator
                    maxlen = max(len(record.seq) for record in records)
                    for record in records:
                        if len(record.seq) != maxlen:
                            sequence = str(record.seq).ljust(maxlen, '.')
                            record.seq = Seq.Seq(sequence)
                    assert all(len(record.seq) == maxlen for record in records)
                    output_file = '{}_4_padded.fasta'.format(os.path.splitext(string_out_4)[0])
                    with open(output_file, 'w') as f:
                        SeqIO.write(records, f, 'fasta')
        make_equal_4()

        def align_seq_4():

                cline = ClustalwCommandline("clustalw2", infile="string4.fasta")
                stdout = cline()
        align_seq_4()

def gene_find_5(gene_5, sub_5):
        for rec in SeqIO.parse(input_handle_5, "genbank"):
                        if rec.features:
                                for feature in rec.features:
                                        if feature.type == "CDS":

                                                if gene_5 in feature.qualifiers['product'][0]:
                                                        if sub_5 in feature.qualifiers['gene'][0]:
                                                                output_handle_5.write(
                                                                ">%s%s %s \n%s\n" % (
                                                                feature.qualifiers['locus_tag'][0],
                                                                feature.qualifiers['gene'][0],
                                                                rec.description,
                                                                feature.location.extract(rec).seq))
        input_handle_5.close()
        output_handle_5.close()

        def change_string_5():
            for rec in SeqIO.parse('outfile_5.fasta', 'fasta'):
                rec.description = re.sub(regex, "_", rec.description)
                string_handle_5.write(">" + rec.description+"\n" + str(rec.seq)+"\n")
            string_handle_5.close()
        change_string_5()

        def make_equal_5():

                    output_string_5 = (string_out_5)
                    records = (rec.upper() for rec in SeqIO.parse(string_out_5, 'fasta'))
                    records = list(records) # make a copy, otherwise our generator
                    maxlen = max(len(record.seq) for record in records)
                    for record in records:
                        if len(record.seq) != maxlen:
                            sequence = str(record.seq).ljust(maxlen, '.')
                            record.seq = Seq.Seq(sequence)
                    assert all(len(record.seq) == maxlen for record in records)
                    output_file = '{}_5_padded.fasta'.format(os.path.splitext(string_out_5)[0])
                    with open(output_file, 'w') as f:
                        SeqIO.write(records, f, 'fasta')
        make_equal_5()

        def align_seq_5():

                cline = ClustalwCommandline("clustalw2", infile="string5.fasta")
                stdout = cline()
        align_seq_5()

def gene_find_6(gene_6, sub_6, rm_6):
        for rec in SeqIO.parse(input_handle_6, "genbank"):
                        if rec.features:
                                for feature in rec.features:
                                    if feature.type == "CDS":
                                            if gene_6 in feature.qualifiers['product'][0]:
                                                    if sub_6 in feature.qualifiers['gene'][0]:
                                                        if rm_6 not in feature.qualifiers['gene'][0]:
                                                                output_handle_6.write(
                                                                ">%s%s %s \n%s\n" % (
                                                                feature.qualifiers['locus_tag'][0],
                                                                feature.qualifiers['gene'][0],
                                                                rec.description,
                                                                feature.location.extract(rec).seq))
        input_handle_6.close()
        output_handle_6.close()

    

        def change_string_6():
            for rec in SeqIO.parse('outfile_6.fasta', 'fasta'):
                rec.description = re.sub(regex, "_", rec.description)
                string_handle_6.write(">" + rec.description+"\n" + str(rec.seq)+"\n")
            string_handle_6.close()
        change_string_6()

        def make_equal_6():

                    output_string_6 = (string_out_6)
                    records = (rec.upper() for rec in SeqIO.parse(string_out_6, 'fasta'))
                    records = list(records) # make a copy, otherwise our generator
                    maxlen = max(len(record.seq) for record in records)
                    for record in records:
                        if len(record.seq) != maxlen:
                            sequence = str(record.seq).ljust(maxlen, '.')
                            record.seq = Seq.Seq(sequence)
                    assert all(len(record.seq) == maxlen for record in records)
                    output_file = '{}_5_padded.fasta'.format(os.path.splitext(string_out_6)[0])
                    with open(output_file, 'w') as f:
                        SeqIO.write(records, f, 'fasta')
        make_equal_6()

        def align_seq_6():

                cline = ClustalwCommandline("clustalw2", infile="string6.fasta")
                stdout = cline()
        align_seq_6()


def convert():
                        alignment = AlignIO.read(open("string.aln"), "clustal", alphabet=IUPAC.unambiguous_dna)
                        g = open("outfile_padded.nexus", "w")
                        g.write (alignment.format("nexus"))

def convert_2():
                        alignment = AlignIO.read(open("string2.aln"), "clustal", alphabet=IUPAC.unambiguous_dna)
                        g = open("outfile_padded_2.nexus", "w")
                        g.write (alignment.format("nexus"))
def convert_3():
                        alignment = AlignIO.read(open("string3.aln"), "clustal", alphabet=IUPAC.unambiguous_dna)
                        g = open("outfile_padded_3.nexus", "w")
                        g.write (alignment.format("nexus"))

def convert_4():
                        alignment = AlignIO.read(open("string4.aln"), "clustal", alphabet=IUPAC.unambiguous_dna)
                        g = open("outfile_padded_4.nexus", "w")
                        g.write (alignment.format("nexus"))

def convert_5():
                        alignment = AlignIO.read(open("string5.aln"), "clustal", alphabet=IUPAC.unambiguous_dna)
                        g = open("outfile_padded_5.nexus", "w")
                        g.write (alignment.format("nexus"))

def convert_6():
                        alignment = AlignIO.read(open("string6.aln"), "clustal", alphabet=IUPAC.unambiguous_dna)
                        g = open("outfile_padded_6.nexus", "w")
                        g.write (alignment.format("nexus"))

def main():
        gene_find('DNA gyrase subunit B', 'gyrB')
        gene_find_2('DNA-directed RNA polymerase subunit beta', 'rpoB')        
        gene_find_3('Translation initiation factor IF-2', 'infB')
        gene_find_4('Transcription termination/antitermination protein NusA', 'nusA')
        gene_find_5('ATP synthase subunit beta', 'atpD')
        gene_find_6('Leucine--tRNA ligase', 'leuS', 'leuS_2')
        convert()
        convert_2()
        convert_3()
        convert_4()
        convert_5()
        convert_6()
main()
