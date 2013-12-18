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


BSD License:
Copyright (c) 2013, Tim Hsiau, J. Christopher Anderson, Regents of UC.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
