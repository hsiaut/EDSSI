import sys, CommonDNA
#IN:   reads in genemark input
#########gmm called with ./gmhmmp -m hmm/Escherichia_coli_K12.mod -o test_out -k -r ../blastulator/12568252 -a
#OUT:  returns List of GenePrediction objects

class GenePrediction:
	#attributes are seq, start, stop, len, similarORFs (list), alignment (obj or matrix?), syntaxscoring (list)
	def __init__( self, name, sequence, start, stop, strand):
		self.name = name
		self.sequence = sequence
		self.start = int(start)
		self.stop = int(stop)
		#self.len = stop-start
		self.strand = strand
		self.ncbiId = ""
	def __str__(self):
		return "\n".join([self.name, self.sequence, str(self.start), str(self.stop), self.strand])

class NCBIProteinEntry:
	def __init__(self, header, sequence):
		self.sequence = sequence
		self.header = header
	def __str__(self):
		return "\n".join([self.sequence, self.header])

def FastaToGenePrediction(genemarkFasta):
	contents = genemarkFasta.split("\n")
	seq = ""
	start = 0
	stop = 0
	strand = ""
	for line in contents:
		if line.startswith(">"):
			fields = line.split("|")
			name = fields[0]
			strand = fields[3]
			start = fields[4]
			tmp = fields[5].split(" ")[0]
			tmp2 = tmp.split("\t")[0]
			stop = tmp2
		else:
			seq = seq+line
	seq = seq.strip()
	g = GenePrediction(name, seq, start, stop, strand)
	g.description = "predicted"
	return g

def ExtendNonMetORFs(predictions, DNAseq):
	for p in predictions:
		if (p.sequence[0] != 'M' and p.sequence[0] != 'V'): #starting amino acid is not M
			#extend upstream in DNAseq until Met found
			nextAA = ""
			i=0
			while(nextAA != 'M'):
				nextcodon = ''
				if p.strand == '+':
					if p.start-i-3 < 0:
						break
					nextcodon = DNAseq[p.start-i-3 : p.start - i]
				else:
					if p.stop + i + 3 > len(DNAseq):
						break
					nextcodon = DNAseq[p.stop + i : p.stop + i + 3]
					nextcodon= CommonDNA.reverseComplement(nextcodon)
                                nextAA = CommonDNA.translate(nextcodon)
				if(nextAA == "_"):
					break
				i+=3
				p.sequence = nextAA + p.sequence
			if p.strand == "+":
				p.start = p.start-i-3
			else:
				p.stop = p.stop + i + 3
	return predictions

def ExtendTruncatedORFs(predictions, DNAseq):
	#takes in a gene prediction, along with its homologs, and sees if there is an upstream sequence in the DNA that matches
	for p in predictions:
		if len(p.homologs) < 1:
			#skip if no homologs
			continue
		else:
			firstTenAA = ""
			try:
				firstTenAA = p.sequence[0:10]
			except IndexError:
				continue
			else:
				#find first ten AA in other homologs
				for h in p.homologs:
					offset = h.sequence.find(firstTenAA)
					if h.sequence.find(firstTenAA) > 1:
						#find extra sequence in our DNA input
						inputMatch = ""
						if p.strand == "+":
							inputMatch = DNAseq[p.start: p.start - offset*3]
						else:
							inputMatch = DNAseq[p.stop: p.stop + offset*3]
							inputMatch = CommonDNA.reverseComplement(inputMatch)
						
	return predictions	

def MergeOverlapping(geneA, geneB):
	#assumed that geneA starts earlier than geneB
	# add on beginning of A to B, and return that object
	loc = geneA.sequence.find(geneB.sequence[:20])
	geneB.sequence = geneA.sequence[:loc]+geneB.sequence
	geneB.start = min(geneA.start, geneB.start)
	geneB.stop = max(geneA.stop, geneB.stop)
	return geneB

#want to take the union of list A and list B.
#rules for merging two together are 
def MergeGenePredObjects(listA, listB):
	results = []
	allmembers = listA+listB
	allmembers.sort(key= lambda x: len(x.sequence), reverse= True)
	ctr = 0
	ignore = []
	for a in allmembers:
		flag = False
		y = ctr + 1
		z = y
		if ctr in ignore:
			ctr = ctr+1
			continue
		for b in allmembers[y:]:
			merge_flag = False
			attempt_merge = False
			if b.sequence in a.sequence:
				merge_flag = True
			#if they overlap, have the same frame, and on the same strand
			elif ((a.start < b.start and  b.start < a.stop)): #or (a.start > b.start and a.start < b.stop)):
				attempt_merge = True
			elif(a.start > b.start and a.start < b.stop):
				attempt_merge = True
			if attempt_merge:
				if(b.sequence[:20] in a.sequence):
					a = MergeOverlapping(a, b)
					merge_flag = True
				elif(a.sequence[:20] in b.sequence): #oppostie strand
					a = MergeOverlapping(b, a)
					merge_flag = True
			if merge_flag:
				ignore.append(z)
				#flag = True  #TODO: copy over the annotations for b onto a
				if(a.description == "predicted"):
					a.description = b.description
				else:
					if(b.description not in a.description):
						a.description += " "+b.description
			z = z +1
		if not flag:
			results.append(a)
		ctr = ctr + 1
	return results
	
def GenemarkReader(genemarkFile):
	listOfGenePredictions = []
	contents=open(genemarkFile)
	buffer = ""
	flag = False
	for line in contents:
		if flag == True and not(line.startswith(">")):
			buffer = buffer+line.strip()
		if line.startswith(">"):
			if flag == True:
				#this is the 2nd or 3rd or.. gene
				#process the buffer and clear it
				listOfGenePredictions.append(FastaToGenePrediction(buffer))
				buffer = ""
			buffer = buffer+line
			flag = True
	listOfGenePredictions.append(FastaToGenePrediction(buffer))
	#for gp in listOfGenePredictions:
		#print gp
		#print "====================\n"
	return listOfGenePredictions

#example usage, import and then call with
#GenemarkReader(sys.argv[1])
