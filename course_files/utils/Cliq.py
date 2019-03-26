#!/usr/bin/python

############################################################
# developed in Python 2.7.3 by Chen Xu  cxu3@uncc.edu      #
# Department of Bioinformatics and Genomics, UNC-Charlotte #       
# Aug 2014                                                 #
############################################################

import re, os, sys, getopt
from collections import defaultdict
import itertools

### step4. format output 
### infile: one line for one clique
### outfile: cluster ID in a column; line Number is the cell index
def output(clique_list=None, outfile=None, number_cells=None):
	if (clique_list is None) or (outfile is None) or (number_cells is None): 
		sys.stderr.write("wrong argument number from output()")
		sys.exit(1)
	cell_clqID={}

	clqID=1
	for line in clique_list:
		cells=line.split(' ')
		cells=map(int, cells)  # index
		if len(cells)<3:
			for c in cells:
				cell_clqID[c]="-1"
		else:
			for c in cells:
				cell_clqID[c]=str(clqID)
			clqID+=1
	try:
		outfn=open(outfile, 'w')
		for c in range(1,number_cells+1):
			if cell_clqID.has_key(c):
				outfn.write(cell_clqID[c]+'\n')
			else:   # some cells have not been covered because no clique found for them.
				outfn.write('-1\n')
		outfn.close()
	except IOError:
		sys.stderr.write("cannot open outfile\n")
		sys.exit(1)
# end 4


### 3. delete overlap. After merging, if one cell appears in multiple cliques (and the cliques do not satisfy merging), select one cliq for the cell to be in. At the same time, delete the cell from the other cliques.
### criteria: select the cliq who has bigger average link weights to the node (by checking the SNN graph)
def uniq(clique_list=None, edgeFile=None, outfile=None):
	# clique_list: clique output from merge() 
	# edgeFile: pair of nodes and weigh of edge (3 column file)
	# outfile: clique output file, each line is a clique
 	if (clique_list is None) or (edgeFile is None): 
		sys.stderr.write("wrong argument number from uniq()")
		sys.exit(1)

	# read in the edge file (node1\tnode2\tweight)
	A=defaultdict(dict)
	try:
		infn=open(edgeFile,'r')
		for line in infn:
			line=line.rstrip('\r\n')
			lst=line.split()
			A[lst[0]][lst[1]]=lst[2]
			A[lst[1]][lst[0]]=lst[2]
		infn.close()
	except IOError:
		sys.stderr.write("cannot open "+edgeFile+"\n")
		sys.exit(1)

	# parse the clique list
	cell_cliq=defaultdict(list)
	cliq_cell=defaultdict(list)
	cliqNum=1
	for line in clique_list:
		cells=line.split(' ')
		for c in cells:
			cell_cliq[c].append(str(cliqNum))
			cliq_cell[str(cliqNum)].append(c)
		cliqNum+=1

#------------------- unique assign---------------------------------------- 
	for c in cell_cliq.keys(): # for each node
		if len(cell_cliq[c])>1:  # if it is in more than one clique
			# count links in each clique associated with the cell
			count_link=defaultdict(list) # key: cliq, value: weighted link in the cliq to the node
			for connect_node in A[c].keys(): # each neighbor of node c
				for cl in cell_cliq[connect_node]: # for each clique that node c in
					if cl in cell_cliq[c]:
						count_link[cl].append(A[c][connect_node])
			# calculate the average of link weights connect from the cliq to the node c		
			max_link={}
			for cl in count_link.keys():
				max_link[cl]=sum(map(float,count_link[cl]))/len(count_link[cl])
			# select the clique to assign the node c
			assign_cliq=max(max_link, key=max_link.get)		
			#print c+" in "+" ".join(cell_cliq[c])+" assign: "+assign_cliq	
			
			# assign to one and delete in other cliqs
			for cl in cell_cliq[c]:
				if not cl==assign_cliq:
					cliq_cell[cl].remove(c)

#--------------- output------------------------------------------
	uniqCliq_list=[]
	for cl in sorted(map(int,cliq_cell.keys())):
		if len(cliq_cell[str(cl)])>0:
			uniqCliq_list.append(" ".join(cliq_cell[str(cl)])+'\n')	

	if outfile:
		try:
			outfn=open(outfile,'w')
			outfn.write("\n".join(uniqCliq_list))
			outfn.close()
		except IOError:
			sys.stderr.write("cannot open "+outfile+"\n")	
			sys.exit(1)
	
	return uniqCliq_list

## end 3.


### step2. merge quasi-clique. Some cluster are long and extended. Need to merge quasi-cliques to find such cluster.
### be aware that after merging, the new clusters generated are bigger and have the potential to be merged again.
### cutoff: the shared must exceed the cutoff percent in one clique to be merged
def merge(clique_list=None, cutoff=None):
	# clique_list: clique output from findquasicliq(). one line is space deliminated indices of a cluster.
	# cutoff: threshold of overlapping rate for merging
	if (clique_list is None):
		sys.stderr.write("wrong argument number from merge()")
		sys.exit(1)
	if cutoff is None:	cutoff=0.5 # 
	cell_cliq=defaultdict(list) # key: cell ID, value: a list of cliq ID(str) that the cell belongs to
	cliq_cell=defaultdict(list) # key: cliq ID(str), value: a list of cells ID in the cliq
	cliqNum=1

	for line in clique_list:
		cells=line.split(' ')
		for c in cells:
			cell_cliq[c].append(str(cliqNum))
			cliq_cell[str(cliqNum)].append(c)
		cliqNum+=1

		
#------------- inner merge cycle---------------------------
	merged_bl=mergeOnline(cutoff, cell_cliq, cliq_cell)
	while merged_bl:
		#for cl in cliq_cell.keys():
		#	print str(cl)+" cliq: "+" ".join(cliq_cell[cl])
		#for c in cell_cliq.keys():
		#	print c+" cell: "+" ".join(map(str,cell_cliq[c]))
		#print "------ new merge begins ------"
		merged_bl=mergeOnline(cutoff, cell_cliq, cliq_cell)

#------------- output --------------------------------------
	mergedCliq_list=[]
	for cl in cliq_cell.keys():
		mergedCliq_list.append(" ".join(cliq_cell[cl]))
	return mergedCliq_list

# 2a. on-line merging algorithm. refresh the cliques and re-do the candidate searching after each merge happens. 
# sort the candidate of merging by the size of cliques, this will allow bigger clique join small cliques first, instead of small cliques join first.
def mergeOnline(cutoff=None, cell_cliq=None, cliq_cell=None):
	# the function will make one possible merge and return

	overlapCliq={} # pairs of cliques that have overlap in between
	# overlapCliq: key is pair of clique names, value is Number of cells overlapping between them
	for c in cell_cliq.keys():
		if len(cell_cliq[c])>1:
			if len(cell_cliq[c])>2: ## pair-wise if 1 cell appear in more than 2 cliques.
				s=list(itertools.combinations(cell_cliq[c],2))
				for pair in s:
					grp=" ".join(pair)
					if overlapCliq.has_key(grp):
						overlapCliq[grp]+=1
					else:
						overlapCliq[grp]=1
			else:
				grp=" ".join(cell_cliq[c])
				if overlapCliq.has_key(grp):
					overlapCliq[grp]+=1
				else:
					overlapCliq[grp]=1

	
	#select overlapped cliques to merge
	# Since one cliq can overlap with multiple cliqs, create a score to sort all possible mergings
	# score is the size of the cliques to merge
	cliq2Merge={}
	for grp in overlapCliq.keys():
		#print "overlap "+grp+": "+str(overlapCliq[grp])
		clqs=grp.split(' ')
		ol=overlapCliq[grp]
		bl=0 # boolean to see if it qualify a merge
		ss=0 # sum of score (sum of overlap proportion in each cliq of the group)
		sum_size=0
		for cl in clqs:
			sum_size+=len(cliq_cell[cl])
			# algo2: overlapping rate only need to exceed cutoff in one cliq
			score=float(ol)/len(cliq_cell[cl])
			ss+=score
			if score>=cutoff:
				bl=1
		if bl==1:
			cliq2Merge[grp]=sum_size
	
	# online merge: merge only one candidate group, then refresh.
	# merge the one with the biggest score
	currCliqNum=max(map(int, cliq_cell.keys()))+1
	if len(cliq2Merge.keys()):
		grp=max(cliq2Merge, key=cliq2Merge.get)	
		#print "merge "+grp+": "+str(cliq2Merge[grp])	
		new={}
		for cl in grp.split(' '):
			for c in cliq_cell[cl]:
				new[int(c)]=''
				# refresh the cell_cliq for on-line learning 
				index=cell_cliq[c].index(cl)
				del cell_cliq[c][index]   # delete old cliq index
				if str(currCliqNum) not in cell_cliq[c]: # add new cliq index
					cell_cliq[c].append(str(currCliqNum)) 
			del cliq_cell[cl]
		cliq_cell[str(currCliqNum)]=map(str,sorted(new.keys()))

	return len(cliq2Merge.keys())
		
# end of mergeInner
## end 2. merge()

### step 1. find quasi-clique
# find a r-quasi-clique associated with each node greedily
# inFile: pair-wise similarity #format: "node_index node_index edge_similarity\n"
# r: cutoff for quasi-clique
def findQuasiCliq(infile=None, r=None):
	if infile is None: 
		sys.stderr.write("wrong argument number for findQuasiCliq()")
		sys.exit(1)
	if r is None:	r=0.7 # 
	minCliqueSize=3 # minimum acceptable clique size

	nodehash=defaultdict(list)
	try:
		infn=open(infile,'r')
		for line in infn:
			line=line.rstrip('\r\n')
			lst=line.split()
			if len(lst)<3 or (not lst[0].isdigit()) or (not lst[1].isdigit()):
				sys.stderr.write("Error: input file not in the right format:\nIndex_of_node1 Index_of_node2 Weight_of_edge\n")
				sys.exit(1)
			try:
				float(lst[2])
			except ValueError:
				sys.stderr.write("Error: input file not in the right format:\nIndex_of_node1 Index_of_node2 Weight_of_edge\n")
				sys.exit(1)
			if float(lst[2])>0:
				nodehash[lst[0]].append(lst[1])
				nodehash[lst[1]].append(lst[0])
		infn.close()
	except IOError:
		sys.stderr.write("Error: cannot open input file "+infile+"\n")	
		sys.exit(1)

#-----------------find clique for each node -------------------
	cliques={}
	
	for node in sorted(nodehash): 
		neighList=list(nodehash[node]) ### get all neighbors
		while 1:
			if len(neighList)<(minCliqueSize-1):
				break
			degreeHash={} ### each neighbor connect to how many of the other neighbors
			for neighNode in neighList:
				neighNeighList=list(nodehash[neighNode])
				degree=len(set(neighNeighList) & set(neighList))+1
				degreeHash[neighNode]=degree

			### greedy: start from the neighbor with the least connection to other neighbors
			minDegNeigh=min(degreeHash, key=degreeHash.get)
            ## %NeighHash length is the number of all neighbors to the Node (including the minDegNeigh itself
			if degreeHash[minDegNeigh]>=int(len(neighList)*r+1):
				###  the local cluster satisfy to be a quasi-clique
				break
			else:
				### does not satisfy, delete the neighbor
				neighList.remove(minDegNeigh)
			
		# quasi-cliq found	
		if len(neighList)>=(minCliqueSize-1):	
			neighList.append(node)
			neighList.sort(key=int)
			cliques[" ".join(neighList)]=len(neighList)

	#for clq in cliques.keys():
	#	print clq

#------------------ delete redundant clique -----------------------
### delete cliques that are part of other cliques
	cliques_clean={}
	for currCliq in sorted(cliques, key=cliques.get, reverse=True):
		bl=1 # bl=1: unique ; bl=0: redundant
		for existCliq in sorted(cliques_clean, key=cliques_clean.get, reverse=True):
			if (cliques_clean[existCliq] <= cliques[currCliq]):
				break
			else:
				bl=0
				currCliqAry=currCliq.split(" ")
				existCliqAry=existCliq.split(" ")
				for win in currCliqAry:
					if win not in existCliqAry:
						bl=1  # unique
						break
				if bl==0:
					break
		if bl:
			cliques_clean[currCliq]=cliques[currCliq]

	max_ID=max(map(int, nodehash.keys()))
	return (cliques_clean.keys(),max_ID)
# end 1.
################       end of algorithm         #########

def usage():
	help_msg="\n  usage: SNNgraph_clustering.py -i <edge_file> -o <out_file> [options]\n\
  -i,--input\tinput file path and name\n\
  -o,--output\toutput file path and name\n\
  optional arguments:\n\
  -r,--r-quasi-cliq\tquasi-clique parameter. A number in the range of (0  1]. Default is 0.7.\n\
  -m,--merging\tmerging parameter. A number in the range of (0  1]. Default is 0.5.\n\
  -n,--number\tnumber of objects. Default is the maximum index in the input file.\n\
  -h,--help\tprint help message.\n"
	print help_msg
	sys.exit(0)

def main(argv):

	edgeFile=None
	outfile=None
	merge_cutoff=0.5
	r_cutoff=0.7
	number_cells=None

	try:
		opts, args=getopt.getopt(argv, "i:o:m:r:n:h", ["input=","output=","merging=","r-quasi-cliq=","number=","help"])
	except getopt.GetoptError as err:
		print str(err)
		usage()
		sys.exit(0)

	for opt, arg in opts:
		if opt in ('-h', '--help'):
			usage()
		elif opt in ('-i','--input'):
			edgeFile=arg
		elif opt in ('-o', '--output'):
			outfile=arg
		elif opt in ('-m','--merging'):
			merge_cutoff=float(arg)
		elif opt in ('-r','--r-quasi-cliq'):
			r_cutoff=float(arg)
		elif opt in ('-n','--number'):
			number_cells=int(arg)
		else:
			print "Error in arguments: unknown option\n"
			usage()
			sys.exit(0)
			
	if (edgeFile is None) or (outfile is None):
		sys.stderr.write("Error in arguments: must specify -i -o\n")
		usage()
	else:
		if not os.path.exists(edgeFile):
			sys.stderr.write("Error: "+edgeFile+" does not exists.\n")
			sys.exit(0)
		try:
			outfn=open(outfile, 'w')
			outfn.close()
		except IOError:
			sys.stderr.write("Error: file path not exists "+outfile+"\n")
			sys.exit(1)
	if (r_cutoff>1) or (r_cutoff<=0) or (merge_cutoff>1) or (merge_cutoff<=0):
		sys.stderr.write("Error: parameters are not in the right range.\n")
		usage()
		sys.exit(0)
	cliques=[]
	results=findQuasiCliq(edgeFile, r_cutoff)
	cliques=results[0]
	max_ID=results[1]
	if number_cells is None:
		number_cells=max_ID
	print "input file "+edgeFile
	print "find "+str(len(cliques))+" quasi-cliques"
	merge_cliques=[]
	merge_cliques=merge(cliques, merge_cutoff)
	print "merged into "+str(len(merge_cliques))+" clusters"
	uniq_cliques=[]
	uniq_cliques=uniq(merge_cliques, edgeFile)
	print "unique assign done \n"
	output(uniq_cliques, outfile, number_cells)

if __name__=="__main__":
	main(sys.argv[1:])

