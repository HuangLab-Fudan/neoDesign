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
# 初始化逻辑说明：
# x_front[pep] 表示可以放在 pep 前面的模块，即 candidate + pep
# x[pep] 表示可以放在 pep 后面的模块，即 pep + candidate
# 因此：
# sorted_length[pep] == 0 不代表 pep 不能用，而是 pep 不能再接后端，适合作为右端终点；
# front_sorted_length[pep] == 0 不代表 pep 不能用，而是 pep 不能再接前端，适合作为左端起点。

def _domain_count(seq):
    with open("sequence.txt", "w") as f:
        f.write(">seq\n")
        f.write(seq)
    res = subprocess.Popen(cmd, shell=True)
    res.wait()
    count = 0
    with open("seq.txt", "r") as f:
        for line in f:
            if line.startswith(">>"):
                count += 1
    return count

def _score(seq, length_penalty):
    sec = util.gor4_cal(seq)
    return util.l * length_penalty + _domain_count(seq) + sec["H"] + sec["E"] - sec["C"]

def _best_left(anchor):
    # 选择一个模块放在 anchor 前面：candidate + anchor
    candidates = [k for k in x_front.get(anchor, []) if k in front_sorted_length]
    if not candidates:
        return None
    best_module = None
    best_score = 999999
    for k in candidates:
        seq = k + anchor
        score = _score(seq, front_sorted_length.get(k, 0))
        print(score)
        if score <= best_score:
            best_score = score
            best_module = k
    return best_module

def _best_right(anchor):
    # 选择一个模块放在 anchor 后面：anchor + candidate
    candidates = [k for k in x.get(anchor, []) if k in sorted_length]
    if not candidates:
        return None
    best_module = None
    best_score = 999999
    for k in candidates:
        seq = anchor + k
        score = _score(seq, sorted_length.get(k, 0))
        print(score)
        if score <= best_score:
            best_score = score
            best_module = k
    return best_module

def _remove_module(m):
    for key in list(x_front.keys()):
        if m in x_front[key]:
            x_front[key].remove(m)
    for key in list(x.keys()):
        if m in x[key]:
            x[key].remove(m)
    front_sorted_length.pop(m, None)
    sorted_length.pop(m, None)

def _choose_right_terminal():
    # 后端不能再接，但前端可接：作为右端终点
    candidates = [
        k for k in sorted_length.keys()
        if sorted_length.get(k, 0) == 0 and front_sorted_length.get(k, 0) > 0
    ]
    if not candidates:
        return None
    # 优先选前端可选数量少的，越受限越早固定
    return min(candidates, key=lambda k: front_sorted_length[k])

def _choose_left_terminal():
    # 前端不能再接，但后端可接：作为左端起点
    candidates = [
        k for k in front_sorted_length.keys()
        if front_sorted_length.get(k, 0) == 0 and sorted_length.get(k, 0) > 0
    ]
    if not candidates:
        return None
    # 优先选后端可选数量少的，越受限越早固定
    return min(candidates, key=lambda k: sorted_length[k])

def _choose_general_front_anchor():
    candidates = [k for k in front_sorted_length.keys() if front_sorted_length[k] > 0]
    if not candidates:
        return None
    return min(candidates, key=lambda k: front_sorted_length[k])

def _choose_general_end_anchor():
    candidates = [k for k in sorted_length.keys() if sorted_length[k] > 0]
    if not candidates:
        return None
    return min(candidates, key=lambda k: sorted_length[k])

# step 1：优先处理只有一侧可连接的端点
right_terminal = _choose_right_terminal()
left_terminal = _choose_left_terminal()

if right_terminal is not None:
    # 例如 SLAEYTDMI：后面不能接，但前面能接
    module2 = right_terminal
    module1 = _best_left(module2)
    if module1 is None:
        raise RuntimeError(f"{module2} has sorted_length=0 but no valid x_front candidate.")
    current_sequence = module1 + module2
    _remove_module(module1)
    _remove_module(module2)

elif left_terminal is not None:
    # 前面不能接，但后面能接
    module1 = left_terminal
    module2 = _best_right(module1)
    if module2 is None:
        raise RuntimeError(f"{module1} has front_sorted_length=0 but no valid x candidate.")
    current_sequence = module1 + module2
    _remove_module(module1)
    _remove_module(module2)

else:
    # step 2：没有端点约束时，按两侧候选数量更少的一侧开始
    first_front_key = _choose_general_front_anchor()
    first_end_key = _choose_general_end_anchor()

    if first_front_key is None or first_end_key is None:
        raise RuntimeError(
            f"No valid initial pair found. front_sorted_length={front_sorted_length}, "
            f"sorted_length={sorted_length}"
        )

    if front_sorted_length[first_front_key] <= sorted_length[first_end_key]:
        module2 = first_front_key
        module1 = _best_left(module2)
        if module1 is None:
            raise RuntimeError(f"No valid front candidate for {module2}.")
    else:
        module1 = first_end_key
        module2 = _best_right(module1)
        if module2 is None:
            raise RuntimeError(f"No valid end candidate for {module1}.")

    current_sequence = module1 + module2
    _remove_module(module1)
    _remove_module(module2)

# step 3：如果模块总数是奇数，初始化阶段需要先形成 3 个模块
# 原代码也是奇数时先构建 3 个模块，后面 while 每轮再左右各加一个。
if len(util.modules) % 2 == 1:
    # 如果当前右端是 terminal，就只能从左边补；
    # 如果当前左端是 terminal，就只能从右边补；
    # 否则按左右候选数量少的一侧补。
    if sorted_length.get(module2, 0) == 0:
        extra = _best_left(module1)
        if extra is None:
            raise RuntimeError(f"Cannot add extra module to the left of {module1}.")
        current_sequence = extra + current_sequence
        module1 = extra
        _remove_module(extra)

    elif front_sorted_length.get(module1, 0) == 0:
        extra = _best_right(module2)
        if extra is None:
            raise RuntimeError(f"Cannot add extra module to the right of {module2}.")
        current_sequence = current_sequence + extra
        module2 = extra
        _remove_module(extra)

    else:
        left_options = len([k for k in x_front.get(module1, []) if k in front_sorted_length])
        right_options = len([k for k in x.get(module2, []) if k in sorted_length])

        if left_options <= right_options:
            extra = _best_left(module1)
            if extra is None:
                extra = _best_right(module2)
                if extra is None:
                    raise RuntimeError("Cannot add third module for odd module count.")
                current_sequence = current_sequence + extra
                module2 = extra
                _remove_module(extra)
            else:
                current_sequence = extra + current_sequence
                module1 = extra
                _remove_module(extra)
        else:
            extra = _best_right(module2)
            if extra is None:
                extra = _best_left(module1)
                if extra is None:
                    raise RuntimeError("Cannot add third module for odd module count.")
                current_sequence = extra + current_sequence
                module1 = extra
                _remove_module(extra)
            else:
                current_sequence = current_sequence + extra
                module2 = extra
                _remove_module(extra)

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
					best_score = score
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
					best_score = score
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
					best_score = score
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
					best_score = score
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
					best_score = score
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
					best_score = score
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
					best_score = score
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
					best_score = score
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
