import re
import threading
import subprocess
import argparse
import multiprocessing
import util

#x_front,front library construction
x_front={}
def library_construction_front(pep,num):
	save=[]
	for i in util.modules:
		if pep!=i:
			seq=i+pep
			antigen_list=util.SB_antigen(seq,num,util.e1,util.e2,util.e3,util.e4,util.e5)
			no_need=util.find_element(antigen_list,util.modules)
			if len(no_need)==0:
				save.append(i)
	x_front[pep]=save
#x,end library construction
x={}
def library_construction_end(pep,num):
	save=[]
	for i in util.modules:
		if pep!=i:
			seq=pep+i
			antigen_list=util.SB_antigen(seq,num,util.e1,util.e2,util.e3,util.e4,util.e5)
			no_need=util.find_element(antigen_list,util.modules)
			if len(no_need)==0:
				save.append(i)
	x[pep]=save

#distribute threads,suitable for less kernels of cpu
def process_data_threaded(m):
	d={}
	for i in range(len(m)):
		d["t"+str(i)]=threading.Thread(target=library_construction_front,args=(m[i],i))
	for i in range(len(m),2*len(m)):
		d["t"+str(i)]=threading.Thread(target=library_construction_end,args=(m[i-len(m)],i))
	for i in d.keys():
		d[i].start()
	for i in d.keys():
		d[i].join()
		
#process threads for front library construction and end library construction
process_data_threaded(util.modules)

#front_sorted_length，sort the front library
length = {}
for k in x_front.keys():
	length[k]=len(x_front[k])
front_sorted_length = dict(sorted(length.items(), key=lambda x: x[1]))

#sorted_length，sort for the end library
length={}
for k in x.keys():
	length[k]=len(x[k])
sorted_length = dict(sorted(length.items(), key=lambda x: x[1]))