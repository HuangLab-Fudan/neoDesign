
NeoDesign: A Computational Tool for Optimal Selection of Polyvalent Neoantigen Combinations
--------------------------------------------
NeoDesign manuscript now published in Bioinformatics. Open Access: https://doi.org/10.1093/bioinformatics/btae585

What is the NeoDesign
----------------
Desiging effective tumor polyvalent neoantigen vaccines poses challenges due to the random combination of peptide segments.
Here we present neoDesign, a tool developed to address these challenges and select an ideal protein sequence for polyvalent tumor neoantigen vaccines. In addition, neoDesign provides a scheme for balancing mRNA stability and protein expression. NeoDesign consists of four modules: Library Construction, Optimal Path Filtering and Linker Addition, and λ-Evaluation. It provides valuable guidance for the selection of ideal protein sequences that can be used in the design of different types of vaccines such as DNA, RNA, peptides and DC. Based on the optimal protein sequences, it can also generate a recommended λ value that can be used in LinearDesign for further mRNA design to achieve a better balance between mRNA stability and protein expression. NeoDesign holds promising applications in the design of tumor polyvalent neoantigen vaccines for different types of vaccines, contributing to improved immunotherapy approaches.

Implementation and Dependencies
-------------------------------

neoDesign was developed with python (recommend>3.9) and shell (bash) language. Before running the program, it is necessary to check or download python packages and local functions as follow:
* gor4
* mhcflurry
* NetMHCpan4.1
* NetChop3.1
* pepsickle
* hmmer(>3.4)

NeoDesign Installation
------------
Follow the installation steps
-------------
# environment construction
# conda create
conda create -n neoDesign
conda activate neoDesign

# install GOR4
pip install gor4

# install mhcflurry(details in https://github.com/openvax/mhcflurry)
pip install mhcflurry
mhcflurry-downloads fetch

# install NetMHCpan4.1
# download NetMHCpan4.1
https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/

# install NetChop3.1
# download NetChop3.1
https://services.healthtech.dtu.dk/services/NetChop-3.1/


# install pepsickle
pip install pepsickle


# download Pfam-A.hmm(new release is 36.0)
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam36.0/Pfam-A.hmm.gz
gzip -d Pfam-A.hmm.gz

# install hmmer (new version is 3.4)
wget http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz
tar xzf hmmer-3.4.tar.gz
./configure --prefix=/yout/install/path
make
make install
# be sure to add hmmer to the environment variable
export PATH=/yout/install/path/bin:$PATH 


Files Needed:
------------
1.peptides.txt(neoantigen peptides file)
2.target_protein_sequence.txt(the optimal protein sequences to generate the recommended λ)

Commands and arguments
--------------------------
* The path in options must be absolute path.
--------------------------
Use neoDesign to generate the optimal protein sequences
usage: python main.py [OPTIONS] 

	-p <peptides file>  
               neoantigen peptides file
	-l <parameter>  
               the parameter for balancing selected numbers and structure(recommand 1, the bigger the numbers, the less linkers the sequence would generate;the smaller the numbers,the less complex structure of the sequence)
	-s <species> 
               species(only support for human and mouse)
	-v <version>
			   the prediction method of netchop:[0] Cterm3.0 [1] 20S-3.0",default=0
	-e <cutoff> 
                the cutoff for domain prediction: recommand (1e-5)
	-e1 <cutoff1> 
                the cutoff of netMHCpanEL,default=2
	-e2 <cutoff2>   
                the cutoff of netMHCpanBA,default=2
	-e3 <cutoff3>
                the cutoff of MHCflurry,default=2
	-e4 <cutoff4>
				the cutoff of netchop,default=0.5
	-e5 <cutoff5>
				the cutoff of pepsickle,default=0.5
	-c <cpu>  
                the kernel numbers of cpu
	-d <path>  
                Pfam directory
	-d2 <path>   
                NetMHCpan4.0 directory
	-d3 <path>
				NetChop3.1 directory



### example ##
	usage: 
	* generate a subdirectory first
	mkdir result
	python3 main.py -p example/TCGA-US-A77E-01A.txt -l 1 -s human -v 0 -e 1e-5 -e1 0.05 -e2 0.05 -e3 0.05 -e4 0.5 -e5 0.5 -c 2 -d /home/neoDesign/Pfam-A.hmm -d2 /home/neoDesign/netMHCpan-4.1/netMHCpan -d3 /home/neoDesign/netchop-3.1/netchop
	attention:
	After each program execution, delete the generated intermediate files in the current program directory and subdirectory result, such as sequence....txt

Output:  
		 result.txt(optimized sequence,the name of the file determined by your input peptides file)
		 ./result(A directory that stores neoantigen predictions)
--------------------------------------------------------------------------------------------------


Use neoDesign to generate the recommended λ
usage: 
	1.Save your optimal protein sequences in a txt file. Each line is a protein sequence and the file is named :target_protein_sequence.txt
	2. Execute the following code:
		python lambda_evaluation.py

Output:  
		 output_lambda.txt(The λ-value corresponds one-to-one with the input protein sequence)
----------------------------------------------------------------------------------------

Supplementary note:
We prepare some example data in the example folder, during execution we need to place the target_protein_sequence.txt file in the same directory as the code, the peptide file can be passed with an absolute path in the input parameter.
