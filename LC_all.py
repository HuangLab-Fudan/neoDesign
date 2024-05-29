import util
import re
import threading
import subprocess
import argparse
import multiprocessing

if multiprocessing.cpu_count()<len(util.modules)*10:
    print("less sockets parallel")
    import Library_Construction_parallel1
    print(Library_Construction_parallel1.front_sorted_length)
    print(Library_Construction_parallel1.sorted_length)
    x_front=Library_Construction_parallel1.x_front
    x=Library_Construction_parallel1.x
    front_sorted_length=Library_Construction_parallel1.front_sorted_length
    sorted_length=Library_Construction_parallel1.sorted_length
else:
    print("more sockets parallel")
    import Library_Construction_parallel2
    print(Library_Construction_parallel2.front_sorted_length)
    print(Library_Construction_parallel2.sorted_length)
    x_front=Library_Construction_parallel2.x_front
    x=Library_Construction_parallel2.x
    front_sorted_length=Library_Construction_parallel2.front_sorted_length
    sorted_length=Library_Construction_parallel2.sorted_length
