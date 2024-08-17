#optimal path filtering
import re
import threading
import subprocess
import argparse
import multiprocessing
from gor4 import GOR4
import util
import LC_all

#main domain command
cmd = "hmmscan -o seq.txt --tblout seq.tbl --noali -E " + util.e + " --cpu " + util.c + " " + util.d + " sequence.txt"
#insert x_front,x,front_sorted_length,sorted_length
x_front=LC_all.x_front
x=LC_all.x
front_sorted_length=LC_all.front_sorted_length
sorted_length=LC_all.sorted_length
#choose the zero situation
first_front_key = list(front_sorted_length.keys())[0]
first_end_key = list(sorted_length.keys())[0]
if front_sorted_length[first_front_key]==0:
    util.move_first_to_last(front_sorted_length)
if sorted_length[first_end_key]==0:
    util.move_first_to_last(sorted_length)

if len(util.modules)%2==1:
	if front_sorted_length[first_front_key]==1 or sorted_length[first_end_key]==1:
		if front_sorted_length[first_front_key] < sorted_length[first_end_key]:
			best_module=None
			best_score=999
			current_sequence=x_front[first_front_key][0]+first_front_key
			module1=x_front[first_front_key][0]
			for key in x_front.keys():
				if first_front_key in x_front[key]:
					x_front[key].remove(first_front_key)
			for key in x.keys():
				if first_front_key in x[key]:
					x[key].remove(first_front_key)
			for key in x_front.keys():
				if module1 in x_front[key]:
					x_front[key].remove(module1)
			for key in x.keys():
				if module1 in x[key]:
					x[key].remove(module1)
			del front_sorted_length[first_front_key]
			del sorted_length[first_front_key]
			del front_sorted_length[module1]
			del sorted_length[module1]
			for key in x[first_front_key]:
				with open ("sequence.txt","w") as f:
					seq=current_sequence+key
					f.write(">seq1")
					f.write("\n")
					f.write(seq)
				res = subprocess.Popen(cmd, shell=True)
				res.wait()
				with open("seq.txt","r") as f:
					i=0
					for line in f:
						if line.startswith(">>"):
							i+=1
				score = util.l*sorted_length[key]+i+util.gor4_cal(seq)["H"]+util.gor4_cal(seq)["E"]-util.gor4_cal(seq)["C"]
				print(score)
				if score <= best_score:
					best_module=key
					score=best_score
			current_sequence=current_sequence+best_module
			module2=best_module
			for key in x_front.keys():
				if module2 in x_front[key]:
					x_front[key].remove(module2)
			for key in x.keys():
				if module2 in x[key]:
					x[key].remove(module2)
			del front_sorted_length[module2]
			del sorted_length[module2]
		
		else:
			best_module=None
			best_score=999
			current_sequence=first_end_key+x[first_end_key][0]
			module2=x[first_end_key][0]
			for key in x_front.keys():
				if first_end_key in x_front[key]:
					x_front[key].remove(first_end_key)
			for key in x.keys():
				if first_end_key in x[key]:
					x[key].remove(first_end_key)
			for key in x_front.keys():
				if module2 in x_front[key]:
					x_front[key].remove(module2)
			for key in x.keys():
				if module2 in x[key]:
					x[key].remove(module2)
			del front_sorted_length[module2]
			del sorted_length[module2]
			del front_sorted_length[first_end_key]
			del sorted_length[first_end_key]
			for key in x_front[first_end_key]:
				with open ("sequence.txt","w") as f:
					seq=key+current_sequence
					f.write(">seq1")
					f.write("\n")
					f.write(seq)
				res = subprocess.Popen(cmd, shell=True)
				res.wait()
				with open("seq.txt","r") as f:
					i=0
					for line in f:
						if line.startswith(">>"):
							i+=1
				score = util.l*front_sorted_length[key]+i+util.gor4_cal(seq)["H"]+util.gor4_cal(seq)["E"]-util.gor4_cal(seq)["C"]
				print(score)
				if score <= best_score:
					best_module=key
					score=best_score
			current_sequence=best_module+current_sequence
			module1=best_module
			for key in x_front.keys():
				if module1 in x_front[key]:
					x_front[key].remove(module1)
			for key in x.keys():
				if module1 in x[key]:
					x[key].remove(module1)
			del front_sorted_length[module1]
			del sorted_length[module1]
	else:
		if front_sorted_length[first_front_key]<sorted_length[first_end_key]:
			best_module=None
			best_score=999
			for key in x_front[first_front_key]:
				with open ("sequence.txt","w") as f:
					seq=key+first_front_key
					f.write(">seq")
					f.write("\n")
					f.write(seq)
				res = subprocess.Popen(cmd, shell=True)
				res.wait()
				with open("seq.txt","r") as f:
					i=0
					for line in f:
						if line.startswith(">>"):
							i+=1
				score = util.l*front_sorted_length[key]+i+util.gor4_cal(seq)["H"]+util.gor4_cal(seq)["E"]-util.gor4_cal(seq)["C"]
				print(score)
				if score <= best_score:
					best_module=key
					score=best_score
			current_sequence=best_module+first_front_key
			module1=best_module
			for key in x_front.keys():
				if first_front_key in x_front[key]:
					x_front[key].remove(first_front_key)
			for key in x.keys():
				if first_front_key in x[key]:
					x[key].remove(first_front_key)
			for key in x_front.keys():
				if module1 in x_front[key]:
					x_front[key].remove(module1)
			for key in x.keys():
				if module1 in x[key]:
					x[key].remove(module1)
			del front_sorted_length[first_front_key]
			del sorted_length[first_front_key]
			del front_sorted_length[module1]
			del sorted_length[module1]
			best_module=None	
			best_score=999
			for key in x[first_front_key]:
				with open ("sequence.txt","w") as f:
					seq=current_sequence+key
					f.write(">seq")
					f.write("\n")
					f.write(seq)
				res = subprocess.Popen(cmd, shell=True)
				res.wait()
				with open("seq.txt","r") as f:
					i=0
					for line in f:
						if line.startswith(">>"):
							i+=1
				score = util.l*sorted_length[key]+i+util.gor4_cal(seq)["H"]+util.gor4_cal(seq)["E"]-util.gor4_cal(seq)["C"]
				print(score)
				if score <= best_score:
					best_module=key
					score=best_score
			current_sequence=current_sequence+best_module
			module2=best_module
			for key in x_front.keys():
				if module2 in x_front[key]:
					x_front[key].remove(module2)
			for key in x.keys():
				if module2 in x[key]:
					x[key].remove(module2)
			del front_sorted_length[module2]
			del sorted_length[module2]
		else:
			best_module=None
			best_score=999
			for key in x[first_end_key]:
				with open ("sequence.txt","w") as f:
					seq=first_end_key+key
					f.write(">seq")
					f.write("\n")
					f.write(seq)
				res = subprocess.Popen(cmd, shell=True)
				res.wait()
				with open("seq.txt","r") as f:
					i=0
					for line in f:
						if line.startswith(">>"):
							i+=1
				score = util.l*sorted_length[key]+i+util.gor4_cal(seq)["H"]+util.gor4_cal(seq)["E"]-util.gor4_cal(seq)["C"]
				print(score)
				best_module=key
				score=best_score
			current_sequence=first_end_key+best_module
			module2=best_module
			for key in x_front.keys():
				if first_end_key in x_front[key]:
					x_front[key].remove(first_end_key)
			for key in x.keys():
				if first_end_key in x[key]:
					x[key].remove(first_end_key)
			for key in x_front.keys():
				if module2 in x_front[key]:
					x_front[key].remove(module2)
			for key in x.keys():
				if module2 in x[key]:
					x[key].remove(module2)
			del front_sorted_length[module2]
			del sorted_length[module2]
			del front_sorted_length[first_end_key]
			del sorted_length[first_end_key]
			best_module=None
			best_score=999
			for key in x_front[first_end_key]:
				with open ("sequence.txt","w") as f:
					seq=key+current_sequence
					f.write(">seq")
					f.write("\n")
					f.write(seq)
				res = subprocess.Popen(cmd, shell=True)
				res.wait()
				with open("seq.txt","r") as f:
					i=0
					for line in f:
						if line.startswith(">>"):
							i+=1
				score = util.l*front_sorted_length[key]+i+util.gor4_cal(seq)["H"]+util.gor4_cal(seq)["E"]-util.gor4_cal(seq)["C"]
				print(score)
				best_module=key
				score=best_score
			current_sequence=best_module+current_sequence
			module1=best_module
			for key in x_front.keys():
				if module1 in x_front[key]:
					x_front[key].remove(module1)
			for key in x.keys():
				if module1 in x[key]:
					x[key].remove(module1)
			del front_sorted_length[module1]
			del sorted_length[module1]


else:
	if front_sorted_length[first_front_key]==1 or sorted_length[first_end_key]==1:
		if front_sorted_length[first_front_key] < sorted_length[first_end_key]:
			best_module=None
			best_score=999
			current_sequence=x_front[first_front_key][0]+first_front_key
			module1=x_front[first_front_key][0]
			module2=first_front_key
			for key in x_front.keys():
				if first_front_key in x_front[key]:
					x_front[key].remove(first_front_key)
			for key in x.keys():
				if first_front_key in x[key]:
					x[key].remove(first_front_key)
			for key in x_front.keys():
				if module1 in x_front[key]:
					x_front[key].remove(module1)
			for key in x.keys():
				if module1 in x[key]:
					x[key].remove(module1)
			del front_sorted_length[first_front_key]
			del sorted_length[first_front_key]
			del front_sorted_length[module1]
			del sorted_length[module1]
		
		else:
			best_module=None
			best_score=999
			current_sequence=first_end_key+x[first_end_key][0]
			module2=x[first_end_key][0]
			module1=first_end_key
			for key in x_front.keys():
				if first_end_key in x_front[key]:
					x_front[key].remove(first_end_key)
			for key in x.keys():
				if first_end_key in x[key]:
					x[key].remove(first_end_key)
			for key in x_front.keys():
				if module2 in x_front[key]:
					x_front[key].remove(module2)
			for key in x.keys():
				if module2 in x[key]:
					x[key].remove(module2)			
			del front_sorted_length[first_end_key]
			del sorted_length[first_end_key]
			del front_sorted_length[module2]
			del sorted_length[module2]
	else:
		if front_sorted_length[first_front_key]<sorted_length[first_end_key]:
			best_module=None
			best_score=999
			for key in x_front[first_front_key]:
				with open ("sequence.txt","w") as f:
					seq=key+first_front_key
					f.write(">seq")
					f.write("\n")
					f.write(seq)
				res = subprocess.Popen(cmd, shell=True)
				res.wait()
				with open("seq.txt","r") as f:
					i=0
					for line in f:
						if line.startswith(">>"):
							i+=1
				score = util.l*front_sorted_length[key]+i+util.gor4_cal(seq)["H"]+util.gor4_cal(seq)["E"]-util.gor4_cal(seq)["C"]
				print(score)
				if score <= best_score:
					best_module=key
					score=best_score
			current_sequence=best_module+first_front_key
			module1=best_module
			module2=first_front_key
			for key in x_front.keys():
				if first_front_key in x_front[key]:
					x_front[key].remove(first_front_key)
			for key in x.keys():
				if first_front_key in x[key]:
					x[key].remove(first_front_key)
			for key in x_front.keys():
				if module1 in x_front[key]:
					x_front[key].remove(module1)
			for key in x.keys():
				if module1 in x[key]:
					x[key].remove(module1)
			del front_sorted_length[first_front_key]
			del sorted_length[first_front_key]
			del front_sorted_length[module1]
			del sorted_length[module1]
		else:
			best_module=None
			best_score=999
			for key in x[first_end_key]:
				with open ("sequence.txt","w") as f:
					seq=first_end_key+key
					f.write(">seq")
					f.write("\n")
					f.write(seq)
				res = subprocess.Popen(cmd, shell=True)
				res.wait()
				with open("seq.txt","r") as f:
					i=0
					for line in f:
						if line.startswith(">>"):
							i+=1
				score = util.l*sorted_length[key]+i+util.gor4_cal(seq)["H"]+util.gor4_cal(seq)["E"]-util.gor4_cal(seq)["C"]
				print(score)
				best_module=key
				score=best_score
			current_sequence=first_end_key+best_module
			module1=first_end_key
			module2=best_module
			for key in x_front.keys():
				if first_end_key in x_front[key]:
					x_front[key].remove(first_end_key)
			for key in x.keys():
				if first_end_key in x[key]:
					x[key].remove(first_end_key)
			for key in x_front.keys():
				if module2 in x_front[key]:
					x_front[key].remove(module2)
			for key in x.keys():
				if module2 in x[key]:
					x[key].remove(module2)
			del front_sorted_length[module2]
			del sorted_length[module2]
			del front_sorted_length[first_end_key]
			del sorted_length[first_end_key]	
#greedy algorithm
if len(util.modules)%2==0:
	a=len(util.modules)/2-1
else:
	a=len(util.modules)/2-2

while a>0:
	best_module=None
	best_score=999
	#less situation
	if len(x_front[module1])<=len(x[module2]):
		for key in front_sorted_length.keys():
			if key in x_front[module1]:
				with open ("sequence.txt","w") as f:
					seq =key+current_sequence
					f.write(">seq1")
					f.write("\n")
					f.write(seq)
				res = subprocess.Popen(cmd, shell=True)
				res.wait()
				with open("seq.txt","r") as f:
					i=0
					for line in f:
						if line.startswith(">>"):
							i+=1
				score = util.l*front_sorted_length[key]+i+util.gor4_cal(seq)["H"]+util.gor4_cal(seq)["E"]-util.gor4_cal(seq)["C"]
				print(score)
				if score <= best_score:
					best_module=key
					score = best_score
		if best_module!=None:
			module1=best_module
			for key in x_front.keys():
				if best_module in x_front[key]:
					x_front[key].remove(best_module)
			for key in x.keys():
				if best_module in x[key]:
					x[key].remove(best_module)
			del front_sorted_length[best_module]
			del sorted_length[best_module]
		else:
			best_score=999
			for key in front_sorted_length.keys():
				with open ("sequence.txt","w") as f:
					seq =key+current_sequence
					f.write(">seq1")
					f.write("\n")
					f.write(seq)
				res = subprocess.Popen(cmd, shell=True)
				res.wait()
				with open("seq.txt","r") as f:
					i=0
					for line in f:
						if line.startswith(">>"):
							i+=1
				score = util.l*front_sorted_length[key]+i+util.gor4_cal(seq)["H"]+util.gor4_cal(seq)["E"]-util.gor4_cal(seq)["C"]
				print(score)
				if score <= best_score:
					best_module=key
					score = best_score
			module1=best_module
			for key in x_front.keys():
				if best_module in x_front[key]:
					x_front[key].remove(best_module)
			for key in x.keys():
				if best_module in x[key]:
					x[key].remove(best_module)
			del front_sorted_length[best_module]
			del sorted_length[best_module]	
		current_sequence=best_module+current_sequence
		best_module=None
		best_score=999
		for key in sorted_length.keys():
			if key in x[module2]:
				with open("sequence.txt","w") as f:
					seq=current_sequence+key
					f.write(">seq1")
					f.write("\n")
					f.write(seq)
				res = subprocess.Popen(cmd, shell=True)
				res.wait()
				with open("seq.txt","r") as f:
					i=0
					for line in f:
						if line.startswith(">>"):
							i+=1
				score = util.l*sorted_length[key]+i+util.gor4_cal(seq)["H"]+util.gor4_cal(seq)["E"]-util.gor4_cal(seq)["C"]
				print(score)
				if score <= best_score:
					best_module=key
					score = best_score
		if best_module!=None:
			module2=best_module
			for key in x_front.keys():
				if best_module in x_front[key]:
					x_front[key].remove(best_module)
			for key in x.keys():
				if best_module in x[key]:
					x[key].remove(best_module)
			del front_sorted_length[best_module]
			del sorted_length[best_module]
		else:
			best_score=999
			for key in sorted_length.keys():
				with open ("sequence.txt","w") as f:
					seq =current_sequence+key
					f.write(">seq1")
					f.write("\n")
					f.write(seq)
				res = subprocess.Popen(cmd, shell=True)
				res.wait()
				with open("seq.txt","r") as f:
					i=0
					for line in f:
						if line.startswith(">>"):
							i+=1
				score = util.l*sorted_length[key]+i+util.gor4_cal(seq)["H"]+util.gor4_cal(seq)["E"]-util.gor4_cal(seq)["C"]
				print(score)
				if score <= best_score:
					best_module=key
					score = best_score
			module2=best_module
			for key in x_front.keys():
				if best_module in x_front[key]:
					x_front[key].remove(best_module)
			for key in x.keys():
				if best_module in x[key]:
					x[key].remove(best_module)
			del front_sorted_length[best_module]
			del sorted_length[best_module]
		current_sequence=current_sequence+best_module
	else:
		best_module=None
		best_score=999
		for key in sorted_length.keys():
			if key in x[module2]:
				with open ("sequence.txt","w") as f:
					seq = current_sequence+key
					f.write(">seq1")
					f.write("\n")
					f.write(seq)
				res = subprocess.Popen(cmd, shell=True)
				res.wait()
				with open("seq.txt","r") as f:
					i=0
					for line in f:
						if line.startswith(">>"):
							i+=1
				score = util.l*sorted_length[key]+i+util.gor4_cal(seq)["H"]+util.gor4_cal(seq)["E"]-util.gor4_cal(seq)["C"]
				print(score)
				if score <= best_score:
					best_module=key
					score = best_score
		if best_module!=None:
			module2=best_module
			for key in x_front.keys():
				if best_module in x_front[key]:
					x_front[key].remove(best_module)
			for key in x.keys():
				if best_module in x[key]:
					x[key].remove(best_module)
			del front_sorted_length[best_module]
			del sorted_length[best_module]
		else:
			best_score=999
			for key in sorted_length.keys():
				with open ("sequence.txt","w") as f:
					seq =current_sequence+key
					f.write(">seq1")
					f.write("\n")
					f.write(seq)
				res = subprocess.Popen(cmd, shell=True)
				res.wait()
				with open("seq.txt","r") as f:
					i=0
					for line in f:
						if line.startswith(">>"):
							i+=1
				score = util.l*sorted_length[key]+i+util.gor4_cal(seq)["H"]+util.gor4_cal(seq)["E"]-util.gor4_cal(seq)["C"]
				print(score)
				if score <= best_score:
					best_module=key
					score = best_score
			module2=best_module
			for key in x_front.keys():
				if best_module in x_front[key]:
					x_front[key].remove(best_module)
			for key in x.keys():
				if best_module in x[key]:
					x[key].remove(best_module)
			if best_module!=None:
				del front_sorted_length[best_module]
				del sorted_length[best_module]
		current_sequence=current_sequence+best_module
		best_module=None
		best_score=999
		for key in front_sorted_length.keys():
			if key in x_front[module1]:
				with open ("sequence.txt","w") as f:
					seq = key+current_sequence
					f.write(">seq1")
					f.write("\n")
					f.write(seq)
				res = subprocess.Popen(cmd, shell=True)
				res.wait()
				with open("seq.txt","r") as f:
					i=0
					for line in f:
						if line.startswith(">>"):
							i+=1
				score = util.l*front_sorted_length[key]+i+util.gor4_cal(seq)["H"]+util.gor4_cal(seq)["E"]-util.gor4_cal(seq)["C"]
				print(score)
				if score <= best_score:
					best_module=key
					score = best_score
		if best_module!=None:
			module1=best_module
			for key in x_front.keys():
				if best_module in x_front[key]:
					x_front[key].remove(best_module)
			for key in x.keys():
				if best_module in x[key]:
					x[key].remove(best_module)
			del front_sorted_length[best_module]
			del sorted_length[best_module]
		else:
			best_score=999
			for key in front_sorted_length.keys():
				with open ("sequence.txt","w") as f:
					seq =key+current_sequence
					f.write(">seq1")
					f.write("\n")
					f.write(seq)
				res = subprocess.Popen(cmd, shell=True)
				res.wait()
				with open("seq.txt","r") as f:
					i=0
					for line in f:
						if line.startswith(">>"):
							i+=1
				score = util.l*front_sorted_length[key]+i+util.gor4_cal(seq)["H"]+util.gor4_cal(seq)["E"]-util.gor4_cal(seq)["C"]
				print(score)
				if score <= best_score:
					best_module=key
					score = best_score
			module1=best_module
			for key in x_front.keys():
				if best_module in x_front[key]:
					x_front[key].remove(best_module)
			for key in x.keys():
				if best_module in x[key]:
					x[key].remove(best_module)
			if best_module!=None:
				del front_sorted_length[best_module]
				del sorted_length[best_module]
		current_sequence=best_module+current_sequence
	a=a-1
print(current_sequence)
