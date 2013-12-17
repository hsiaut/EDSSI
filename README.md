EDSSI
=====
Paper abstract:
DNAs encoding for polypeptides often contain design errors that cause experiments to prematurely fail.  One class of design errors is incorrect or missing elements in the DNA, referred to here as syntax errors.  We have identified three major causes of syntax errors: point mutations from sequencing or manual data entry, gene structure misannotation, and unintended ORFs.  EDSSI is an online bioinformatics pipeline that checks for syntax errors through three steps.   First, ORF prediction in input DNA sequences is done by GeneMark; next, homologous sequences are retrieved by BLAST; and finally, syntax errors in the protein sequence are predicted by using the SIFT algorithm.  We show that EDSSI is able to identify previously published examples of syntatical errors and also show that our indel addition to the SIFT program is 97% accurate on a test set of E. coli proteins.  EDSSI is available at http://andersonlab.qb3.berkeley.edu/Software/EDSSI/


An API service is available at  http://andersonlab.qb3.berkeley.edu/Software/EDSSI/syntaxchecker.php

example curl command

curl -X POST -F "api=1" -F "inputSeq=&lt;myfilename" http://andersonlab.qb3.berkeley.edu/Software/EDSSI/syntaxchecker.php

Where myfilename is an input fasta file. The < in front of the file name will tell curl to read the file contents and append it to the POST string.  The API will time out after 10 minutes.




Software Notes:

Required binaries (not supplied):

blastp  (BLAST 2.2.25+)
rpstblastn (BLAST 2.2.25+)
http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

gmhmmp  (GeneMark version 2.5m) and E. coli model (hmm/Escherichia_coli_K12.mod)
http://www.genepro.com/Products/PrGeneMarkhmm.aspx

SIFT (5.0.3) 
http://sift.bii.a-star.edu.sg/
(requirement: BLIMPS) 
http://blocks.fhcrc.org/blocks/uploads/blimps/
see http://blocks.fhcrc.org/blocks/uploads/blimps/blimps.patch.txt  for notes on compiling


MIT License:
Copyright (c) 2013 Tim Hsiau, J. Christopher Anderson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
