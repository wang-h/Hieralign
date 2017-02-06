#  -*- encoding: utf-8 -*-
# 　Copyright 2013 by Hao Wang  
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
import numpy as np
from matplotlib import font_manager
import matplotlib.pyplot as plt
from collections import defaultdict
from textwrap import wrap
import math
import random
import numpy as np
import linecache
import logging
logging.basicConfig( level=logging.WARN)
font_ja = font_manager.FontProperties(fname="PathTo/fonts-japanese-gothic.ttf")
font_zh = font_manager.FontProperties(fname="PathTo/chinese.msyh.ttf")

def line2list(line):
	if line:
		return [x.lower() for x in line.rstrip("\n").split()]
	else:
		return []
def align2matrix(alignment_line, ls, lt):
	l = alignment_line.rstrip("\n").split()
	matrix = [[1 for _ in range(lt)] for _ in range(ls)]
	for align_pair in l:
		i, j = map(int, align_pair.split('-'))
		#print (ls,lt,i, j)
		matrix[i][j] = 0
	return np.asarray(matrix)
def read_tt(tt_file):
	tt=defaultdict(dict)
	for line in tt_file:
		parts = line.rstrip("\n").split() 
		if len(parts)==3: 
			sw,tw,score=line.rstrip("\n").split() 
			tt[sw][tw]=float(score)

	return tt 
def read_e2f_f2e(e2f_file, f2e_file):
	tt=defaultdict(dict)
	for line1,line2 in zip(e2f_file, f2e_file):
		sw,tw,pts=line1.rstrip("\n").split()
		pst=float(line2.rstrip("\n").split()[2]) 
		tt[sw][tw]=(float(pts)*pst)

	return tt 
def build_softmatrix(sws, tws, tt, theta):
	matrix = []
	for sw in sws:
		row = []
		for tw in tws:
			if tw in tt[sw]:
				row.append(1-math.exp(math.log(tt[sw][tw])/theta))
			else:
				row.append(1)
		matrix.append(row ) 
	image_matrix = np.asarray(matrix)
	return image_matrix
def build_softmatrix_with_phi(sws, tws, tt, theta, phi):
	matrix = []
	ls = len(sws)
	lt = len(tws)
	for i,sw in enumerate(sws):
		row = []
		for j,tw in enumerate(tws):
			if tw in tt[sw]:
				row.append(1-math.exp(math.log(tt[sw][tw])/theta)*math.exp(math.log(1-abs((i+1)/ls-(j+1)/lt))/phi))
			else:
				row.append(1)
		matrix.append(row ) 
	image_matrix = np.asarray(matrix)
	return image_matrix

def readLine(corpus, alignFiles, n=-1): 
	if alignFiles:
		for i, ((ss, ts ), aligns) in enumerate(zip(corpus, zip(*alignFiles))):
			if i == n:
				return (ss, ts , aligns)
				break
	else:
		for i, (ss, ts ) in enumerate(corpus):
			if i == n:
				return (ss, ts , None)
				break
	return (None, None , None)
def readParallelLine(corpus, alignFiles, n=-1): 
	print(corpus, alignFiles, n)
	if alignFiles:
		for i, (pair, aligns) in enumerate(zip(*corpus, zip(*alignFiles))):
			if i == n:
				ss, ts = pair.split("\t")
				return (ts, ss , aligns)
				break
	else:
		for i, pair in enumerate(corpus[0]):  
			if i == n:

				ss, ts = pair.split("\t")
				return (ss, ts , None)
				break
	return (None, None , None)
def plot(tt_file, corpus, index, alignFiles):
	subgraph_titles =  [x.name for x in  args.alignFiles]
	if len(corpus)==2:
		ss, ts , aligns = readLine(corpus, alignFiles, index)
	elif len(corpus)==1:
		ss, ts , aligns = readParallelLine(corpus, alignFiles, index)
	print(ss, ts , aligns)
	matrices_count  = len(alignFiles)
	sws, tws = map(line2list, [ss, ts ]) 
	ls = len(sws)
	lt = len(tws)
	if aligns:
		align_matrices = [y for y in map(lambda x:align2matrix(x, ls, lt), aligns )]
	else:
		align_matrices =[]
	#print (tt_file)
	if tt_file :
		if len(tt_file)==1:
			tt    = read_tt(tt_file[0])

		elif len(tt_file)==2:
			tt    = read_e2f_f2e(tt_file[0],tt_file[1])
		subgraph_titles = ["soft alignment matrix"] + subgraph_titles
		soft_matrix    = build_softmatrix(sws, tws, tt, 1)
		matrices = [soft_matrix] + align_matrices
		

		soft_matrix1    = build_softmatrix(sws, tws, tt, 3)
		soft_matrix2    = build_softmatrix(sws, tws, tt, 5) 
		
		matrices.extend([soft_matrix1,soft_matrix2])
		soft_matrix3    = build_softmatrix_with_phi(sws, tws, tt, 5, 0.1) 
		soft_matrix4    = build_softmatrix_with_phi(sws, tws, tt, 5, 0.5) 
		soft_matrix5    = build_softmatrix_with_phi(sws, tws, tt, 5, 1) 

		matrices.extend([soft_matrix3,soft_matrix4,soft_matrix5])
		for matrix in matrices:
			print(matrix)
	else:
		matrices = align_matrices
	 
	if len(matrices) == 1:
		fig, ax = plt.subplots(nrows=1, ncols=1)
		subgraph_list = [ax]
	elif len(matrices) == 2:
		fig, ((ax, ax1) ) = plt.subplots(nrows=1, ncols=2)
		subgraph_list = [ax, ax1]
	elif len(matrices) == 3:
		fig, ((ax, ax1, ax2) ) = plt.subplots(nrows=1, ncols=3)
		subgraph_list = [ax, ax1, ax2]
	elif len(matrices) == 4:
		fig, ((ax, ax1),(ax2, ax3) ) = plt.subplots(nrows=2, ncols=2)
		subgraph_list = [ax, ax1, ax2, ax3]
	elif len(matrices) == 6:
		fig, ((ax, ax1, ax2),(ax3, ax4, ax5) ) = plt.subplots(nrows=2, ncols=3)
		subgraph_list = [ax, ax1, ax2, ax3, ax4, ax5]
	else:
		logging.WARN("TOO MANY ALIGNMENT FILES!" )	
	#print 
	#subgraph_titles = ["soft alignment matrix", "GIZA++", "fast_align", "Hieralign"]
	subgraph_titles = ["[$\sigma_{\Theta}=1$]", "[$\sigma_{\Theta}=3$]", "[$\sigma_{\Theta}=5$]"]
	subgraph_titles += ["[$\sigma_{\Theta}=5, \sigma_{\delta}=0.1$]", "[$\sigma_{\Theta}=5, \sigma_{\delta}=0.5$]", "[$\sigma_{\Theta}=5, \sigma_{\delta}=1$]"]
	subplot(ls, lt, sws, tws, matrices, subgraph_list, subgraph_titles)	
	fig.tight_layout()
	#plt.grid()
	plt.show()
def subplot(ls, lt, sws, tws, align_matrices, subgraph_list, subgraph_titles):	
	for matrix, ax, title in zip(align_matrices, subgraph_list, subgraph_titles):
		#print (matrix, ax, title)
		ax.imshow(matrix, cmap= plt.cm.gray, interpolation='nearest')
		#ax.spines['left'].set_position(('outward', 10))
		#ax.spines['bottom'].set_position(('outward', 12))
		# Hide the right and top spines
		#ax.spines['right'].set_visible(False)
		#ax.spines['top'].set_visible(False)
		ax.set_title(title) 
		# 	# Only show ticks on the left and bottom spines
		ax.yaxis.set_ticks_position('left')
		ax.xaxis.set_ticks_position('bottom')
		ax.xaxis.set_ticks([i for i in range(lt)])
		ax.yaxis.set_ticks([i for i in range(ls)])
		ax.xaxis.set_ticks([i+0.52 for i in range(lt)], minor = True)
		ax.yaxis.set_ticks([i+0.52 for i in range(ls)], minor = True)
		ax.yaxis.set_ticklabels([ '\n'.join(wrap(l, 8)) for l in  sws], rotation=340)
		ax.xaxis.set_ticklabels(tws, fontproperties = font_ja, rotation=90)
		#ax.grid(which='both')                                                            

		# or if you want differnet settings for the grids:                                             
		ax.grid(which='minor', alpha=0.5)  	
		
if __name__ == '__main__': 
	import argparse
	parser = argparse.ArgumentParser(description="""
	This is an python script for extracting ngram entries for building a memory-based machine translation system (also, ebmt and abmt).
	""")
	parser.add_argument('--tt', action='store', dest='tt', type=argparse.FileType('r'),
					help='translation table')
	parser.add_argument('--e2f', action='store', dest='e2f', type=argparse.FileType('r'),
					help='e2f table')
	parser.add_argument('--f2e', action='store', dest='f2e', type=argparse.FileType('r'),
					help='f2e table')
	parser.add_argument('-s', action='store', dest='source', type=argparse.FileType('r'),
					help='source file')
	parser.add_argument('-t', action='store', dest='target', type=argparse.FileType('r'),
					help='target file')
	parser.add_argument('--parallel', action='store', dest='parallel', type=argparse.FileType('r'),
					help='parallel corpus file')
	parser.add_argument('-a', action='append', dest='alignFiles', type=argparse.FileType('r'),
                    default=[],
                    help='Add alignment file to a list')
	parser.add_argument('-i', action='store', dest='index', type=int,
					help='plot the i-th alignment matrix')
	parser.add_argument('--version', action='version', version='%(prog)s 0.1')
	args   = parser.parse_args() 
	corpus = None 
	
	if (args.source and args.target):
		corpus =[open(args.source),open(args.source)]
	if (args.parallel):
		corpus = [args.parallel]
	print(corpus)
	if args.e2f and  args.f2e:
		plot([args.e2f,args.f2e], corpus, args.index, args.alignFiles)
	if args.tt: 
		plot([args.tt], corpus, args.index, args.alignFiles)
