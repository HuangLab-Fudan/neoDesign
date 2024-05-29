#main
import util
import re
import threading
import subprocess
import argparse
import multiprocessing

parser = argparse.ArgumentParser()
parser.add_argument("-p","--pep",help="neoantigen peptides file")
parser.add_argument("-l","--parameter",help="the parameter balancing the library numbers and structure")
parser.add_argument("-s","--species",help="species")
parser.add_argument("-v","--version",help="the prediction method of netchop:[0] Cterm3.0 [1] 20S-3.0",default=0)
parser.add_argument("-e1","--cutoff1",help="the cutoff of netMHCpanEL,defaule=2")
parser.add_argument("-e2","--cutoff2",help="the cutoff of netMHCpanBA,default=2")
parser.add_argument("-e3","--cutoff3",help="the cutoff of MHCflurry,default=2")
parser.add_argument("-e4","--cutoff4",help="the cutoff of netchop,default=0.5")
parser.add_argument("-e","--cutoff",help="the cutoff of hmmerscan,recommended:1e-5")
parser.add_argument("-c","--cpu",help="cpu numbers")
parser.add_argument("-d","--directory",help="Pfam directory")
parser.add_argument("-d2","--directory2",help="NetMHCpan4.0 directory")
parser.add_argument("-d3","--directory3",help="Netchop3.1 directory")
args = parser.parse_args()

if __name__=='__main__': 
    import add_linker
    print(add_linker.x)
    print(add_linker.x_front)
    print(add_linker.front_sorted_length)
    print(add_linker.sorted_length)
    if add_linker.result_linker:
        print(add_linker.result_linker)
    else:
        print(add_linker.final)
