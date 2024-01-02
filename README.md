# neoDesign
NeoDesign: a program for optimizing tumor polyvalent neoantigen vaccine
--------------------------------------------
NeoDesign manuscript now published in XXX. Open Access: www.XXX

What is the NeoDesign
----------------
Desiging effective tumor polyvalent neoantigen vaccines poses challenges due to the random combination of peptide segments.
Here, we present neoDesign,a tool developed to address these challenges and select an ideal protein sequence for tumor polyvalent neoantigen vaccines. NeoDesign contains three modules, Library Construction, Optimal Path Filtering and Linker Addition.It provides valuable guidance for the selection of ideal sequences. NeoDesign holds promising applications in the design of tumor polyvalent neoantigen vaccines for different types of vaccines, contributing to improved immunotherapy approaches.

Implementation and Dependencies
-------------------------------

neoDesign was developed with python (recommend>3.9) and shell (bash) language. Before running the program, it is necessary to check or download python packages and local functions as follow:
* gor4
* mhcflurry
* NetMHCpan4.1
* hmmer(>3.4)

NeoDesign Installation
------------

## environment construction
## conda create
conda create -n neoDesign
conda activate neoDesign

## install GOR4
pip install gor4

## install mhcflurry(details in https://github.com/openvax/mhcflurry)
pip install mhcflurry
mhcflurry-downloads fetch

## download NetMHCpan4.1（details in the website）
https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/

## download Pfam-A.hmm(new release is 36.0)
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam36.0/Pfam-A.hmm.gz
gzip -d Pfam-A.hmm.gz

## install hmmer (new version is 3.4)
wget http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz
tar xzf hmmer-3.4.tar.gz
./configure --prefix=/yout/install/path
make
make install
## be sure to add hmmer to the environment variable
export PATH=/yout/install/path/bin:$PATH 


Files Needed:
------------
1.peptides.txt(neoantigen peptides file)

Commands and arguments
--------------------------
* The path in options must be absolute path.
--------------------------
Use neoDesign
usage: python main.py [OPTIONS] 

	-p <peptides file>  
               neoantigen peptides file
	-l <parameter>  
               the parameter for balancing selected numbers and structure(recommand 1, the bigger the numbers, the less linkers the sequence would generate;the smaller the numbers,the less complex structure of the sequence)
	-s <species> 
               species(only support for human and mouse)
	-e <cutoff> 
                the cutoff for domain prediction: recommand (1e-5)
	-e1 <cutoff1> 
                the cutoff of netMHCpanEL,defaule=2
	-e2 <cutoff2>   
                the cutoff of netMHCpanBA,default=2
	-e3 <cutoff3>
                the cutoff of MHCflurry,default=2
	-c <cpu>  
                the kernel numbers of cpu
	-d <path>  
                Pfam directory
	-d2 <path>   
                NetMHCpan4.0 directory


### example ##
	usage: 
 	* The results folder must be created in the current directory before running for the first time.
	mkdir result
	python3 main.py -p MC38.txt -l 1 -s mouse -e 1e-5 -c 2 -d /home/neoDesign/Pfam-A.hmm -d2 /home/neoDesign/netMHCpan-4.1/netMHCpan
	* attention:
	* After each program execution, delete the generated intermediate files in the current program directory and subdirectory result, such as sequence....txt，seq....txt
 	rm ./result/*
  	rm seq*.txt
   	rm sequence*.txt

Output:  
		 result.txt(optimized sequence)
		 ./result(A directory that stores neoantigen predictions)
--------------------------------------------------------------------------------------------------

