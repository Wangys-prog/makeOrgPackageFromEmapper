#python
import sys
from Bio import SeqIO

from optparse import OptionParser

def MakeOption():
    # make option
    parser = OptionParser(usage="%prog [-h] [-v] -i[--input=] -o[--output=]",
                          version="%prog 1.0")
    parser.add_option("-i", "--input", action="store", dest="input",
                      help="the generalists or specialists Orthogroup id  ",
                      default=False)
    parser.add_option("-o", "--output", action="store", dest="output",
                      help="the protein sequences with fasta format after pick out",
                      default=False)
    (options, args) = parser.parse_args()
    input = options.input
    output = options.output
    return (input,output)

def main():
    input,output = MakeOption()
    output=open(output,"w")
    for eachline in SeqIO.parse(input,"fasta"):
        ids=eachline.description
        ids2 = ids.strip().split(">")[0].split(" ")[0]
        output.write("%s\n"%(ids2))
    output.close()
if __name__ == '__main__':
    main()

