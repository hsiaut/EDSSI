import os
from GenePredictionReader import NCBIProteinEntry

#IN:  takes in list of geneprediction obj`s, each obj has a list of similarORFs
# aligns using muscle or your favorite alinger
#OUT:  updates the geneprediction obj     geneprediction.alignment  contains the alignment

def AlignHomologs(geneprediction):
	stablepy_path = "stable.py"
	buffer = geneprediction.name + "\n" + geneprediction.sequence
	for h in geneprediction.homologs:
		buffer += "\n" + h.header + "\n" + h.sequence
	#write out to temp file and call muscle
	f = open("muscle.tempin", "w")
	f.write(buffer)
	f.close()
	os.system("muscle -quiet -in muscle.tempin -out muscle.tempout")	
	#call stable.py 
	os.system("python "+stablepy_path+" muscle.tempin muscle.tempout > muscle.stable")
	#read in the alignment
	f = open("muscle.stable")
	alignment = []
	for l in f:
		if l.startswith(">"):
			header = l.rstrip()
		else:
                        fields = header.split("|")
			ncbi_id = ""
			description = ""
			if len(fields)>3:
				ncbi_id = fields[1]
				description = fields[4]	
			fastaline = l.rstrip()
			q = NCBIProteinEntry(header, fastaline)
			q.ncbi_id = ncbi_id
			q.description = description
			alignment.append(q)
	os.remove("muscle.stable")
	os.remove("muscle.tempin")
	os.remove("muscle.tempout")
	return alignment

#removes positions in alignment where only one non-ref sequence is found
def CleanAlignment(geneprediction):
	#the first sequence in alignment array is the input seq
	truncationposition = 0
	ctr = 0
	for char in geneprediction.alignment[0].sequence:
		if char != "-":
			break
		sumcount = 0
		for h in geneprediction.alignment[1:]:
			if h.sequence[ctr] != "-":
				sumcount += 1
		ctr += 1
		if sumcount > 1:
			truncationposition = ctr-1
			break
	#copy over everything in alignment after truncation position
	newalign = []
	for h in geneprediction.alignment:
		h.sequence = h.sequence[truncationposition:]
	return geneprediction.alignment

#strips every position in alignment where there is a gap in the ref seq
def StripAlignment(geneprediction):
	tempaln = [""]*len(geneprediction.alignment)
	ctr = 0
	for char in geneprediction.alignment[0].sequence:
		if char != "-":
			#non gap character in ref seq, copy over to temp alignment
			p = 0
			for m in geneprediction.alignment:
				tempaln[p] += m.sequence[ctr]
				p += 1
		ctr += 1
	return tempaln
