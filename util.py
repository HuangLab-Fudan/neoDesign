import re
import threading
import subprocess
import argparse
import multiprocessing
from gor4 import GOR4
#自定义列表
default_linker_library = ["GGGGSGGGGS","SGGGGSGGGG","GSGSGSGSGS","GGSGGSGGSGGS","GGGGSGGGGSAAA","GGGGSGAAAGGSGGGG","GSGSSGSGSS"]
parser = argparse.ArgumentParser()
parser.add_argument("-p","--pep",help="neoantigen peptides file")
parser.add_argument("-l","--parameter",help="the parameter balancing the library numbers and structure")
parser.add_argument("-s","--species",help="species")
parser.add_argument("-v","--version",help="the prediction method of netchop:[0] Cterm3.0 [1] 20S-3.0",default=0)
parser.add_argument("-e1","--cutoff1",help="the cutoff of netMHCpanEL,defaule=2")
parser.add_argument("-e2","--cutoff2",help="the cutoff of netMHCpanBA,default=2")
parser.add_argument("-e3","--cutoff3",help="the cutoff of MHCflurry,default=2")
parser.add_argument("-e4","--cutoff4",help="the cutoff of netchop,default=0.5")
parser.add_argument("-e5","--cutoff5",help="the cutoff of pepsickle,default=0.5")
parser.add_argument("-e","--cutoff",help="the cutoff of hmmerscan,recommended:1e-5")
parser.add_argument("-c","--cpu",help="cpu numbers")
parser.add_argument("-d","--directory",help="Pfam directory")
parser.add_argument("-d2","--directory2",help="NetMHCpan4.0 directory")
parser.add_argument("-d3","--directory3",help="Netchop3.1 directory")
parser.add_argument(
        '--linker_library',
        type=str,
        nargs='*',  # 允许多个值
        default=default_linker_library,
        help='List of linker sequences (default: use built-in ones).'
    )
args = parser.parse_args()

#basic functions
#split_string
def split_string(string, positions):
    result = []
    start = 0
    for pos in positions:
        if pos <= len(string):
            result.append(string[start:pos])
            start = pos
    result.append(string[start:])
    return result
#SB_antigen:find neoantigens
def SB_antigen(seq,num,e1=2,e2=2,e3=2,e4=0.5,e5=0.5,v=0):
    filename="./result/new_antigen"+str(num)+".fas"
    output="./result/antigen_out"+str(num)+".txt"
    output2="./result/predictions"+str(num)+".csv"
    output3="./result/chop"+str(num)+".out"
    output4="./result/pepsickle"+str(num)+".out"
    with open(filename,"w") as f:
        f.write(">seq")
        f.write("\n")
        f.write(seq)
    if args.species=="mouse":
        cmd = args.directory2 + " -f " + filename + " -BA -xls -a H-2-Db,H-2-Kb > " + output
    else:
        cmd = args.directory2 + " -f " + filename + " -BA -xls >" + output
    res = subprocess.Popen(cmd,shell=True)
    res.wait()
    antigen=[]
    with open(output,"r") as f:
        for line in f.readlines()[52:]:
            if line.startswith("  ") and line.strip():
                modified_line = re.sub(r'\s+',' ',line.strip())
                x=modified_line.split(" ")
            if float(x[12])<float(e1) or float(x[14])<float(e2):
                antigen.append(x[9])
    if args.species=="mouse":
        cmd2="mhcflurry-predict-scan  --sequences "+seq+ " --alleles H-2-Kb,H-2-Db --threshold-affinity-percentile "+ e3 + " --out " +output2
    else:
        cmd2="mhcflurry-predict-scan  --sequences "+seq+ " --alleles HLA-A*02:01,HLA-A*03:01,HLA-B*57:01,HLA-B*45:01,HLA-C*02:02,HLA-C*07:02 --threshold-affinity-percentile "+ e3 +" --out " +output2
    res = subprocess.Popen(cmd2,shell=True)
    res.wait()
    with open(output2,"r") as f:
        for line in f.readlines()[1:]:
            antigen.append(line.strip().split(",")[2])
    cmd3 = args.directory2 + " " + filename + " -t " + str(e4) + " -v " + str(v) + " -inptype 1 > " + output3
    res = subprocess.Popen(cmd3,shell=True)
    res.wait()
    positions=[]
    with open(output3,"r") as f:
        for line in f.readlines()[20:-5]:
            modified_line = re.sub(r'\s+',' ',line.strip())
            x=modified_line.split(" ")
            if x[2]=="S":
                positions.append(int(x[0]))
    for i in split_string(seq, positions):
        if len(i)>7 and len(i)<12:
            antigen.append(i)
    cmd4 =  "pepsickle -f " + filename + " -t " + str(e5) + " -o " + output4
    res = subprocess.Popen(cmd4,shell=True)
    res.wait()
    positions=[]
    with open(output4,"r") as f:
        next(f)
        for line in f:
            modified_line = re.sub(r'\s+', ' ', line.strip())
            x=modified_line.split(" ")
            if len(x) > 4 and x[3] == "True":
                positions.append(int(x[0]))
    for i in split_string(seq, positions):
        if len(i) > 7 and len(i) < 12:
            antigen.append(i)                   
    return antigen
#find unexpected neoantigens
def find_element(A,B):
	not_in_B=[]
	for element in A:
		if all(element not in b for b in B):
			not_in_B.append(element)
	return not_in_B

def gor4_cal(seq):
	secondary = {}
	gor4 = GOR4()
	result = gor4.predict(seq)
	structure = result["predictions"]
	total_count=len(seq)
	H_count=((structure.count("H")+structure.count("G")+structure.count("I"))/total_count)*10
	C_count=((structure.count("C")+structure.count("S"))/total_count)*10
	E_count=((structure.count("E")+structure.count("B"))/total_count)*10
	secondary["H"]=H_count
	secondary["C"]=C_count
	secondary["E"]=E_count
	return secondary


#optimal path filtering
def move_first_to_last(dictionary):
    first_key = next(iter(dictionary))
    first_value = dictionary[first_key]
    del dictionary[next(iter(dictionary))]
    dictionary[first_key]=first_value

#add linker
def in_B(A, B):
    in_B_substrings = []
    for element in A:
        longest_substring = ''
        substrings = [element[i:j] for i in range(len(element)) for j in range(i + 1, len(element) + 1)]
        for substr in substrings:
            if any(substr in b for b in B) and len(substr) > len(longest_substring):
                longest_substring = substr
        in_B_substrings.append(longest_substring)
    list2 = list(set(in_B_substrings))
    return list2

def check_position(A, B):
    x_values = {}
    for item_A in A:
        found = False
        for item_B in B:
            if item_B.startswith(item_A):
                x_values[item_A] = 0
                found = True
                break
            elif item_B.endswith(item_A):
                x_values[item_A] = 1
                found = True
                break
        if not found:
            x_values[item_A] = -1
    return x_values

def check_index(A, sequence):
    index_value = {}
    save = []
    for item_A in A:
        index = sequence.find(item_A)
        index2 = sequence.find(item_A) + len(item_A)
        if index not in save and index2 not in save:
            index_value[item_A] = 0
            save.append(index)
            save.append(index2)
        else:
            index_value[item_A] = 1
    return index_value

def find_and_replace(A, x_values, index_values, linker, sequence):
    index_list = []
    for element in A:
        index = sequence.find(element)
        index_list.append(index)
        if x_values[element] == 0 and index_values[element] == 0:
            sequence = sequence[:index] + linker + sequence[index:]
        elif x_values[element] == 1 and index_values[element] == 0:
            sequence = sequence[:index + len(element)] + linker + sequence[index + len(element):]
    return sequence
#get modules
with open(args.pep,"r") as f:
    modules=[]
    i=0
    for line in f:
        line = line.strip()
        modules.append(line)
        i+=1
        if i>30:
            break   


d=args.directory
l=int(args.parameter)
s=args.species
e=args.cutoff
c=args.cpu
e1=args.cutoff1
e2=args.cutoff2
e3=args.cutoff3
e4=args.cutoff4
e5=args.cutoff5
v=int(args.version)
p=args.pep
linker_library = args.linker_library if args.linker_library else default_linker_library