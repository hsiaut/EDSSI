import BlastHandler, json

x = BlastHandler.efetchLinkedPubmed([12429098,10995770])

#x = json.dumps(x)

#print x 

for d in  x['PubmedArticleSet']['PubmedArticle']:
	#print d
	pmid = d['MedlineCitation']['PMID']['#text']
	print pmid
	title = d['MedlineCitation']['Article']['ArticleTitle']
	print title
	abstract = d['MedlineCitation']['Article']['Abstract']['AbstractText']
	print abstract
#print x	
