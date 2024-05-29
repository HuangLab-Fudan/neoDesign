#add_linker
import re
import threading
import subprocess
import argparse
import multiprocessing
from gor4 import GOR4
import util
import optimal_path_filtering
x_front=optimal_path_filtering.x_front
x=optimal_path_filtering.x
front_sorted_length=optimal_path_filtering.front_sorted_length
sorted_length=optimal_path_filtering.sorted_length
sequence=optimal_path_filtering.current_sequence
new_antigen_list = util.find_element(util.SB_antigen(optimal_path_filtering.current_sequence,0,util.e1,util.e2,util.e3), util.modules)
linker_library = ["GGGGSGGGGS","SGGGGSGGGG","GSGSGSGSGS","GGSGGSGGSGGS","GGGGSGGGGSAAA","GGGGSGAAAGGSGGGG","GSGSSGSGSS"]
result_linker = {}
if len(new_antigen_list)==0:
        final=sequence
        print("result")
        print(final)
else:
        for i in linker_library:
                A = util.in_B(new_antigen_list,util.modules)
                x = util.check_position(A,util.modules)
                y = util.check_index(A,sequence)
                new_sequence = util.find_and_replace(A,x,y,i,sequence)
                if len(util.find_element(util.SB_antigen(new_sequence,0,util.e1,util.e2,util.e3),util.modules))==0:
                        final = new_sequence
                        result_linker[i]=new_sequence
                        print("linker",i)
                        print(final)
filename=str(util.p).split("/")[-1]
filename_2=filename.split(".")[0]+"_result.txt"
with open(filename_2,"w") as f:
        if result_linker:
            f.write("linker")
            f.write("\t")
            f.write("sequence")
            f.write("\n")
            for i in result_linker.keys():
                    f.write(i)
                    f.write("\t")
                    f.write(result_linker[i])
                    f.write("\n")
        else:
            f.write(final)
