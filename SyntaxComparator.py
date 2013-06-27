#IN:  geneprediction obj and cutoff score
#DOES:  goes through the alignment and checks if the reference/input disagrees with the related sequences
#OUT: updates geneprediction.syntaxscoring  this is a list of potentially wrong positions in the predicted ORF
import operator, re, math, os, csv
from operator import mul
from pprint import pprint
from StringIO import StringIO
def runSIFT(geneprediction):
	siftpath = "sift"
	##write file
	tempfile = "stripfile.aln"
	f = open(tempfile, "w")
	i = 0
	for k in geneprediction.cleanedaln:
		print >> f, geneprediction.alignment[i].header
		print >> f, k
		i += 1
	f.close()
	##run sift
	siftpath = "sift5.0.3/"
	cmd = "env BLIMPS_DIR="+siftpath+"blimps/blimps-3.8/ "+siftpath+"bin/info_on_seqs "+tempfile+" - siftout"
        os.system(cmd)
	##read sift output
	filename = "siftout"
	with open(filename, 'rb') as csvfile:
      		#get rid of muliple spaces
       		io = StringIO(csvfile.read().replace('\s+', ' '))
       		reader = csv.reader(io, delimiter=' ')
        	header={}
        	score_matrix = []
        	aln_pos = 0
        	for row in reader:
                	if row[0] == '':
                        	if row[2] == 'A': #header row
                                	continue
                        	else:
                                	score_list = []
                                	for cell in row:
                                       		if cell != '':
                                                	score_list.append(float(cell))
                                	score_matrix.append(score_list)
	##correlate with input sequence
	#order of scores  in scores_list is   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z   *   -
	score_key ={
	'A':0,
	'B':1,
	'C':2,
	'D':3,
	'E':4,
	'F':5,
	'G':6,
	'H':7,
	'I':8,
	'K':9,
	'L':10,
	'M':11,
	'N':12,
	'P':13,
	'Q':14,
	'R':15,
	'S':16,
	'T':17,
	'V':18,
	'W':19,
	'X':20,
	'Y':21,
	'Z':22,
	'*':23,
	'-':24,
	}
	ref_seq = geneprediction.sequence
	sift_scores = []
	for i,c in enumerate(ref_seq):
	        try:
	                sift_scores.append(score_matrix[i][score_key[c]])
	        except:
	                print "failed at", i, c
	                continue
	return sift_scores

def SyntaxCompare(geneprediction, cutoff):
	#sift scores has a list of residue by residue scores, layout is relative to the non-gapped input ref sequence
	sift_scores = runSIFT(geneprediction)
	print sift_scores
	tempfile = "tempfile.aln"
	f = open(tempfile, "w")
	for k in geneprediction.alignment:
		print >> f, k.header
		print >> f, k.sequence
	f.close()

	#for each position in the ALIGNMENT, we will calculate a score
	debug = False
	alignment = []
	for k in geneprediction.alignment:
		alignment.append(k.sequence)
	position = 0
	ref_gene_position = 0
	markedupstring = ""
	scoring_track = [1]*len(alignment[0]) #keeps the score for the alignment, will add to the geneprediction obj
        ##substitute non-gaps with sift scores
        idx = 0
        sift_idx = 0
        for c in alignment[0]:
                if c != '-':
                        scoring_track[idx] = sift_scores[ sift_idx ]
                        sift_idx += 1
                idx += 1

	error_flag = False #becomes True if there is an error anywhere
	errors = {} #dict with key = position, value = the correct amino acid
	variance_score = 0 #used to compute how divergent the sequences are
	aln_idx = 0
	for c in alignment[0]:
		#look at the other characters and compute something
		D = 1
		N_c = len(alignment)
		reference_char = c
		N_ca_counts = {
		'A': 1,
		'X': 1, 
		'R': 1, 
		'N': 1, 
		'D': 1, 
		'C': 1, 
		'Q': 1, 
		'E': 1, 
		'G': 1, 
		'H': 1, 
		'I': 1, 
		'L': 1, 
		'K': 1, 
		'M': 1, 
		'F': 1, 
		'P': 1, 
		'S': 1, 
		'T': 1, 
		'W': 1, 
		'Y': 1, 
		'V': 1, 
		'-': 1}
		
		for h in range(1, len(alignment)): #iterate through each aligned sequence
			#formula is Prob(amino acid a @ position c) = N_c/(N_c+B_c)*g_ca + B_c/(N_c+B_c)
			# N_c = # seq in aln
			# B_c = larger if more variance, also some empirical scaling factor
			# B_c = e^D
			# D = (product of N_ca, N_ca are padded by pseudocounts of 1)^(1/N_c)
			ca_char = alignment[h][position]
			if ca_char in N_ca_counts:
				N_ca_counts[ca_char] += 1
			else:
				#ignore
				1+1
				print "ignored ", ca_char
			if debug:
				print "in aln, char is ", ca_char
		counts = N_ca_counts.values()
		if debug:
			print "----------===========position", position
			print "counts", counts
		#count of most common character
		max_count = max(counts)
		#most common character
		mostcommonchar = ""
		for char, count in N_ca_counts.iteritems():
		    if count == max_count:
        		mostcommonchar = char
		
		product = reduce(mul, counts) #calculates the product of the list of counts
		D = product**(1/(N_c)) 
		B_c = 2.718**(D)
		#simple ratio g_ca
		g_ca = float( N_ca_counts[reference_char] - 1 ) / N_c
		if g_ca < 0.05:
			print "g_ca is ", g_ca, " reference_char is ", reference_char, " position is ", position, " num1 ", N_ca_counts[reference_char] - 1, " N_c ", N_c
		pre = ""
		post = ""
		if position > 0:
			pre = alignment[0][position-1]
		if position < len(alignment[0])-1:
			post = alignment[0][position+1]
		if debug:
			print "reference_char", pre.lower() + reference_char + post.lower()
			print alignment[0][position]
			print "g_ca", g_ca
			print "B_c", B_c
			print "D", D
			print "product", product
			print "that was N_ca_counts"
	
		if(reference_char is not "-"):
			p_ca = scoring_track[aln_idx] #retain sift score
		else:
			#for now, use simple ratio for gaps
			p_ca = g_ca
		p_ca = round(p_ca, 3)
		scoring_track[aln_idx] =p_ca
		if debug:
			print "p_ca", p_ca
			print "====================================="
		if(reference_char != '-'):
			ref_gene_position = ref_gene_position +1 
		position += 1
			#if non-gap character, use sift score instead
		#the first alignment is always the reference input
		##seeing if the current position is incorrect
		if p_ca < cutoff:
			error_flag = True
			majority_char = max(N_ca_counts.iteritems(), key=operator.itemgetter(1))[0]
			if(ref_gene_position in errors):
				errors[ref_gene_position] = errors[ref_gene_position] + majority_char+str(p_ca)
			else:
				errors[ref_gene_position] = majority_char+str(p_ca)
			#also, we have to build a string that is used for representing the errors
			#in this string, the error characters are lower case and wrong gaps are = instead of -
			if re.match("^[A-Z]$", reference_char):
				markedupstring += reference_char
			else:
				markedupstring += "-"
		else:
			#no error, copy over the char as normal
			markedupstring += reference_char
		aln_idx += 1
	ref_len = len(alignment[0])
	variance_score = variance_score / ref_len
	variance_score = math.pow(10, variance_score)
	geneprediction.markedupstring = markedupstring					
	##substitute non-gaps with sift scores
	#idx = 0
	#sift_idx = 0
	#for c in alignment[0]:
	#	if c != '-':
	#		scoring_track[idx] = sift_scores[ sift_idx ]
	#		sift_idx += 1
	#	idx += 1	 
	geneprediction.scoring_track = scoring_track 
	geneprediction.errors = errors
	if debug:
		print "scoring_track", scoring_track
		print "error_flag", error_flag
		print "errors"
		pprint(errors)	
