import os, hashlib, sys, re, xmltodict
from urllib import urlopen
from GenePredictionReader import NCBIProteinEntry, GenePrediction
from CommonDNA import gencode, reverseComplement, translate, complement_alphabet
#IN:  list of geneprediction Objs

#OUT: updated geneprediction Objs.  geneprediction.similarORFs  is a list containing blast hits to that prediction


def md5sum(contents):
	d = hashlib.md5()
	d.update(contents)
	return d.hexdigest()

def FastaToObj(fasta):
        contents = fasta.split("\n")
	name = ""
        seq = ""
        for line in contents:
                if line.startswith(">"):
			name = line
                else:
                        seq = seq+line
        seq = seq.strip()
	z = NCBIProteinEntry(name, seq)
	return z

def efetchProteinFasta(ncbiId):
	url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id="
        url += ncbiId + "&rettype=fasta&retmode=text"
	return urlopen(url).read()

def efetchLinkedPMIDs(proteinIds):
	url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=protein&db=pubmed&id="
	for id in proteinIds:
		url += str(id) + ","
	url = url[:-1]
	try:
		xml = urlopen(url).read()
		xml = xml.split("LinkSetDb",1)[1]
		hound = re.compile('([0-9]+)')
		m = hound.findall(xml)
		m = list(set(m))
		return m
	except:
		print "error in efetchLinkedPMIDs"
		return ""

def efetchLinkedPubmed(pmids):
	searchList = ""
	for pmid in pmids:
		searchList += str(pmid) +","
	searchList = searchList[:-1]
	#url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=text&email=hsiaut@gmail.com&tool=paris&id="+searchList
	url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=text&rettype=xml&id="+searchList
	xml = urlopen(url).read()
	try:
		doc = xmltodict.parse(xml)
	except:
		doc = False
	return doc

def efetchLinkedPubmedInfos(pmids):
	print "these are my pmids", pmids
	x = efetchLinkedPubmed(pmids)
	#print x
	results = []
	if(not x):
		return results
	for d in  x['PubmedArticleSet']['PubmedArticle']:
        	try:
			pmid = d['MedlineCitation']['PMID']['#text']
	        	try:
				title = d['MedlineCitation']['Article']['ArticleTitle']
			except KeyError:
				title = "Unavailable"
			try:
				abstract = d['MedlineCitation']['Article']['Abstract']['AbstractText']
			except KeyError:
				abstract = "Unavailable"
		except:
			pmid="Error"
			title="Unavailable"
			abstract="Unavailable"
		r = {'pmid':pmid, 'title':title, 'abstract':abstract}
		results.append(r)
	return results
def convertPMIDtoJSON(pmids):
	results = []
	for pmid in pmids:
		r = {'pmid':pmid}
		results.append(r)
	return results
#other method is to use blastdbcmd, if locally installed, this is usually much faster,
# strategy is to call this first, then resort to eutils if blast is not locally installed
# or protein is not in downloaded copy of db
def blastdbcmdProteinFasta(ncbiId):
	filename = "blastdbcache/"+ncbiId
	cmd = "./blastdbcmd -db nr -target_only -entry " + ncbiId + " > "+filename
	os.system(cmd)
	f = open(filename, "r")
	contents = f.read()
	return contents
	#protein not found error
	#Error: 51: OID not found
	#return an Exception
	#try/catch no command found exception
	
#IN sequence and blast commands
#OUT list of lists, each line is a row in the output (sans comments)
#cachetype is a textual descriptor of which blast search was done, it will showup in the cache file name
def blastCmdHelper(sequence, command, test_flag, outputdir, filename, cachetype):
	#todo: change the tempin tempout names to a unique hash
	results = []
	uniqueName = md5sum(sequence)
	blastcache = "blastcache"
        tempin = blastcache+"/tempin"+uniqueName+filename
        tempout = blastcache+"/"+cachetype+uniqueName+filename
	#check if tempout exists, else do the blast
	if not os.path.exists(tempout):
       		f = open(tempin, "w")
        	f.write(">temp\n")
        	f.write(sequence)
        	f.close()
		command_string = command+" -outfmt 7 -query "+tempin+" > "+tempout
		if not(test_flag):
			os.system(command_string)
	contents=open(tempout)
        for line in contents:
                if line.startswith("#"):
			#todo, add code that regex's for the RID
                        continue
                k = line.split("\t")
		results.append(k)
        #os.remove(tempin)
        #os.remove(tempout)
	return results

#IN:  sequence
#OUT: list of homologs
def FindHomologs(sequence, test_flag, outputdir, database):
	homologs = []
	num_hits = 1000
	blastresults = blastCmdHelper(sequence, "./blastp -db "+database+" -remote -evalue 1e-25 -max_target_seqs "+str(num_hits), test_flag, outputdir, "findhomologs", database+"blast")
	for k in blastresults:
		alnlen = int(k[3]) #length of blast alignment
		alnid = float(k[2]) #percent identity of blast alignment
		overallmatch = alnid*alnlen/len(sequence)
		if k[1].startswith("gi") and (overallmatch> 35):
			ncbiId = k[1].split("|")[1]
			description = k[1].split("|")[4]
			homologs.append(ncbiId)

        subsampleflag = False
	subsamples =[]
	#take even subsamples
	if subsampleflag:
		subsamples =[]
       		num_comparators = 30
       		stepsize = int( len(homologs) / 30 )
		if stepsize <1:
			stepsize = 1
	        idx = 0
       		while (idx < len(homologs) ):
               		subsamples.append( homologs[idx] )
                	idx += stepsize	
	else:
		slice_size = 55
		if len(homologs) > slice_size:
			subsamples=homologs[:slice_size]
		else:
			subsamples=homologs
	homologs = []
	for z in subsamples:
		ncbiId = z
		fasta_text = blastdbcmdProteinFasta(ncbiId)
		if "Error" in fasta_text or len(fasta_text) < 10:
			#not found locally, check ncbi server
			print >> sys.stderr, 'checking ncbi server for fasta ', ncbiId
			fasta_text = efetchProteinFasta(ncbiId)
		if len(fasta_text) < 10:
			print >> sys.stderr, 'bad data for ', ncbiId
		z = FastaToObj(fasta_text)
		z.ncbiId = ncbiId
		#z.description = description.rstrip()
		homologs.append( z)
	return homologs


#looks up description of conserved domains from local file
def CDLookUp(ncbiId):
	cmd = "cat cddid_all.tbl | awk '{if(/^"+str(ncbiId)+"\t/) print $3;}'"
	p = os.popen(cmd)
	s = p.readline()
	p.close()
	return s.rstrip()

#IN: sequence
#OUT: list of gene prediction objects, generated from rps blast (identifies CDD hits)
def FindConservedDomains(sequence, test_flag, outputdir):
	conservedDomains=[]
	cmd = "./rpstblastn -remote -db cdd -evalue 1e-25"
	blastresults = blastCmdHelper(sequence, cmd, test_flag, outputdir, "conserveddomains", "cddblast")
	for k in blastresults:
		if k[1].startswith("gnl"):
			cddId = k[1].split("|")[2]
			start = int(k[6])
			stop = int(k[7])  #outfmt 7 reports start > stop if on reverse strand
                        #for matching to a protein entry, you must take the blast coordinates minus 1, example is 207:1058  becomes fastastr[206:1057]
			if (start > stop):
				temp = int(start)
				start = int(stop-1)
				stop = int(temp)
				strand = "-"
			else:
				start = start-1
				strand = "+"
			#we should get rid of a couple of the terminal amino acids, most of the time they don't match
			start = start + 9
			stop = stop - 9
			gene = sequence[start:stop]
			dnaPreSeq = str(gene)
			if(strand == "-"):
				gene = reverseComplement(gene)
			protein = translate(gene)
			z = GenePrediction(">"+k[1], protein, start, stop, strand)
			z.ncbiId = cddId
			z.dnaPreSeq = dnaPreSeq
			z.dnaSeq = gene
			z.description = CDLookUp(cddId)
			#ignore or split results that contain stop codon
			if(not "_" in protein):
				conservedDomains.append(z)
	return conservedDomains
			
#example usage
#homologs = FindHomologs("MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRIDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFFHNMGDHVTRLDRWEPELNEAIPNDERDTTMPVAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW")
#for h in homologs:
#	print h
