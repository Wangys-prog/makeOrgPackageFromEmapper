#!/usr/bin/python
# File created on 2019/8/19
__author__ = "Wang,Yansu"
__copyright__ = "Copyright 2019, JSNU"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Wang,yansu"
__email__ = "Wangys_c@hotmail.com"

from optparse import OptionParser
from Bio import SeqIO

"""
-i Bg_contraction.txt -f Orthogroups.txt -s genome_all.fasta -o Bg_contraction.fasta
"""
def MakeOption():
    # make option
    parser = OptionParser(usage="%prog [-h] [-v] -i[--input=] -f [--ref=] -s [--sequence=]-o[--output=]",
                          version="%prog 1.0")

    parser.add_option("-i", "--input", action="store", dest="input",
                      help="the generalists or specialists Orthogroup id  ",
                      default=False)
    parser.add_option("-f", "--ref", action="store", dest="ref",
                      help="Orthogroups containing Orthogroup id and sequence id",
                      default=False)
    parser.add_option("-s", "--sequence", action="store", dest="seq",
                      help="all sequence with fasta format",
                      default=False)
    parser.add_option("-o", "--output", action="store", dest="output",
                      help="the protein sequences with fasta format after pick out",
                      default=False)
    (options, args) = parser.parse_args()
    # extract option from command line
    input = options.input
    ref = options.ref
    seq = options.seq
    output = options.output
    return (input,ref,seq,output)

def get_sequence_id(input, ref):
   """
   根据specialist and generalist orthogroups id 得到序列的序列号
   """
   dict = {}
   list1 = []
   for line1 in open(input,"r"):
       each1 = line1.split()[0].strip()
       dict[each1] = {}
   for line2 in open(ref,"r"):
       each2 = line2.split(":")[0].strip()
       if each2 in dict:
          dict[each2]=line2.split(":")[1].strip()
   for value in dict.values():
       list1.append(str(value).split())
       list2 = [i for k in list1 for i in k]
   return list2

def get_seq_according_to_id(list2,seq,output):
    output = open(output, "w")
    seq_dict = {}
    a_dict = {}
    extract_seq_dict = {}
    for eachline in SeqIO.parse(seq, "fasta"):
        ids = eachline.id
        ids_description = eachline.description
        seqs = eachline.seq
        seq_dict[ids] = str(seqs)
        a_dict[ids_description] = str(seqs)
    for i in range(len(list2)):
            if str(list2[i].strip())in seq_dict :
               extract_seq_dict[str(list2[i].strip())] = seq_dict[str(list2[i].strip())]
               new_dict = {v: k for k, v in a_dict.items()}
               query_seq = extract_seq_dict[str(list2[i].strip())]
               if query_seq in new_dict:
                   ids3 = new_dict[query_seq]
                   output.write(">%s\n%s\n" % (ids3,query_seq))
    output.close()

def main():
    input, ref,seq,output = MakeOption()
    list2 = get_sequence_id(input, ref)
    get_seq_according_to_id(list2, seq, output)
if __name__ == '__main__':
    main()