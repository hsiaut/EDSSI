from GenePredictionReader import GenemarkReader, MergeGenePredObjects, ExtendNonMetORFs, ExtendTruncatedORFs
from BlastHandler import FindHomologs, FindConservedDomains, efetchLinkedPMIDs, efetchLinkedPubmedInfos, convertPMIDtoJSON
from AlignerHandler import AlignHomologs, CleanAlignment, StripAlignment
from SyntaxComparator import SyntaxCompare
import sys, os, json, threading, Queue
from ReadFasta import fasta_read
from time import time
from collections import defaultdict
#master process
start = time()
debug = False
genemark_path = "genemark_suite_linux_64/"
inputfile = sys.argv[1]
outputdir = sys.argv[2]

inputData = fasta_read(inputfile)
inputSeq = inputData[0].seq
if not os.path.exists(outputdir):
        os.makedirs(outputdir)

def printRunTime(start):
	return round(time()-start, 4), "s"
if debug:
	print "running genemark...", printRunTime(start)
tempgeneout = outputdir +"/temp.geneout"
genemarkcmd = genemark_path + "gmhmmp -m " + genemark_path + "hmm/Escherichia_coli_K12.mod -o "+tempgeneout+" -k -a -r "+inputfile
os.system(genemarkcmd)
#run genemark
if debug:
	print "parsing genemark file...", printRunTime(start)
predictions = GenemarkReader(tempgeneout)

#run cdd (rpstblastn)

conservedDomains = FindConservedDomains(inputSeq, False, outputdir)
#list of ncbiproteinEntry objs, #seq has protein sequence
#for c in conservedDomains:
#	print "found conserved domain, ", c.ncbiId, c.description
#run BlastHandler
if debug:
	print "found ", len(predictions), " predicted genes"
	print "blasting for homologs...", printRunTime(start)
outStr = ""
for g in conservedDomains:
        outStr += str(g) + "\n dnaseq" + g.dnaSeq +"\n"+g.dnaPreSeq+"\n"+g.description+"\n"
f = open(outputdir+"/out.consDomain", "w")
f.write(outStr)
f.close()



predictions = MergeGenePredObjects(predictions, conservedDomains)
#if any ORFs do not start with a M, extend 5'
predictions = ExtendNonMetORFs(predictions, inputSeq)

outStr = ""
for g in predictions:
	outStr += str(g)
f = open(outputdir+"/out.genepreds", "w")
f.write(outStr)
f.close()

class ThreadRemoteBlast(threading.Thread):
	def __init__(self, queue):
		threading.Thread.__init__(self)
		self.queue = queue
	def run(self):
		while True:
			#grab obj from queue
			g = self.queue.get()
			#runs FindHomologs
			homologs = FindHomologs(g.sequence, False, outputdir, "nr")
			g.homologs = homologs
			#signal task is done
			self.queue.task_done()
queue = Queue.Queue()
for i in range(5):
	t = ThreadRemoteBlast(queue)
	t.setDaemon(True)
	t.start()
for g in predictions:
	queue.put(g)
queue.join()

#extend truncated ORFs
predictions = ExtendTruncatedORFs(predictions, inputSeq)

if debug:
	print "aligning homologs...", printRunTime(start)
for g in predictions:
	g.alignment = AlignHomologs(g)
	g.alignment = CleanAlignment(g)
	g.cleanedaln = StripAlignment(g)
for g in predictions:
	if debug:
		print "syntaxcompare call\n"
	SyntaxCompare(g, 0.05)
	
	#at this point, if there are errors, we go back and look at the input DNA to see if the gene prediction made an error

#detect if there is an exact match
exactMatchFlag = 0
for g in predictions:
	g.exactMatchFlag = 0 #if 0, no exact match, else there is an exact match
	ctr = 1
	for h in g.homologs:
		if h.sequence == g.sequence:
			g.exactMatchFlag = ctr
		ctr+=1

#now each predicton has a scoring_track associated with it, as well as a has_error flag
#build the JSON out
#TODO: parcel this out to a different function or file
jsons = []
js_codes = []
for g in predictions:
	alignment_fastas = []
	alignment_fasta_names = []
	fasta_names=[]
	ncbi_ids=[]
	ctr = 0
	for k in g.alignment:
		if ctr <1:
			alignment_fastas.append(g.markedupstring)
			ctr+= 1
		else:	
			alignment_fastas.append(k.sequence)
			ctr += 1
		alignment_fasta_names.append(k.description)
		fasta_names.append(k.header)
		ncbi_ids.append(k.ncbi_id)
	protein_filtered_ids = filter(None, ncbi_ids) #get rid of empty string/False values
	if(protein_filtered_ids):
		pmids = efetchLinkedPMIDs(protein_filtered_ids)
		#if(pmids):
		#	pubmedInfos = convertPMIDtoJSON(pmids)
		#else:
		#	pubmedInfos =[]
		pubmedInfos = []
	else:
		pubmedInfos = []
	#find most common description
	nameFreq = defaultdict(int)
	for n in alignment_fasta_names:
                description = n.split("[")[0].rstrip()
                nameFreq[description] += 1
	result = max(nameFreq.iteritems(), key=lambda x: x[1])
	g.consensus_name = result[0].strip()
	#1 if exact match exists, 2 if db hits but no exact match, 3 if no db match
	ORFtype = 0
	if g.exactMatchFlag:
		ORFtype = 1
	elif len(g.alignment) > 1:
		ORFtype = 2
	else:
		ORFtype = 3
	data = [ { 'ORFtype': ORFtype, 'exact_match': g.exactMatchFlag,'description': result[0], 'start' : g.start, 'stop' : g.stop,'orientation': g.strand, 'sequence': g.sequence, 'errors': g.errors, 'alignment' : {'scoring_track':g.scoring_track, 'pubmedInfos':pubmedInfos, 'num_seqs': len(g.alignment), 'sequences':alignment_fastas, 'names':alignment_fasta_names, 'protein_ids':ncbi_ids, 'fasta_names':fasta_names}} ]
	jsons += (data)

if not os.path.exists(outputdir):
	os.makedirs(outputdir)



####output functios
def write_gbfile(predictions, outputdir, inputSeq):
	featuresText = ""
	featureColors={'chloramphenicol acetyltransferase':'blue', 'beta-lactamase':'red', "aminoglycoside 3'-phosphotransferase":'green', 'streptomycin3''-adenylyltransferase': 'yellow'}
	for prediction in predictions:
		color="cyan"
		try:
			color = featureColors[prediction.consensus_name]
		except:
			color = "cyan"
		location_info = "     CDS             "+str(prediction.start)+".."+str(prediction.stop)
		label = "                     /label="+prediction.consensus_name
		colors = "                     /ApEinfo_fwdcolor="+color
		colors += "\n                     /ApEinfo_revcolor="+color
		featuresText += location_info+"\n"+label+"\n"+colors+"\n"
	templateText="LOCUS       pdc261                  1008 bp ds-DNA     linear       01-MAY-2013\nDEFINITION  .\nACCESSION   \nVERSION     \nSOURCE      .\n  ORGANISM  .\nCOMMENT     .\nFEATURES             Location/Qualifiers\n"
	outputtext = templateText + featuresText + "ORIGIN\n" + inputSeq +"\n//"
	f = open(outputdir+"/out.gbk", "w")
	f.write(outputtext)
	f.close()
#make fasta file
buffer = ""
for geneprediction in predictions:
	buffer += "\n======================"
	#buffer += "\n>Input\n"+geneprediction.sequence
	for h in geneprediction.alignment:
		buffer += "\n" + h.header + "\n" + h.sequence
f = open(outputdir+"/out.fa", "w")
f.write(buffer)
f.close()

#write annotated "genbank" file
write_gbfile(predictions, outputdir, inputSeq)

#write out.js
f = open(outputdir+"/out.js", "w")
f.write('syntaxcheckeroutput=')
#write main part of json
inputLength = len(inputData[0].seq)
inputName = inputData[0].header
hasErrorFlag = True
f.write('{"inputName":"'+inputName+'", "length":' + str(inputLength) + ', "has_error":"' + str(hasErrorFlag) + '", "ORFs":')
#dump list portion of json
f.write(json.dumps(jsons))
f.write("}")
f.write(';')
f.close()
#output main html page
#main_text = "<html> hello</html>"
#q = open(outputdir+"/index.html", "w")
#q.write(main_text)
#q.close()

sys.stdout.write('\a')
sys.stdout.flush()
