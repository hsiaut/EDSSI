
gencode = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

# Description: converts DNA string to amino acid string
def translate( sequence ):
        sequence = sequence.upper()
        """Return the translated protein from 'sequence' assuming +1 reading frame"""
        return ''.join([gencode.get(sequence[3*i:3*i+3],'X') for i in range(len(sequence)//3)])

complement_alphabet = {'A':'T', 'T':'A', 'C':'G', 'G':'C','R':'Y', 'Y':'R',
                                              'W':'W', 'S':'S', 'M':'K', 'K':'M', 'H':'D', 'D':'H',
                                              'B':'V', 'V':'B', 'N':'N','a':'t', 'c':'g', 'g':'c',
                                              't':'a', 'r':'y', 'y':'r', 'w':'w', 's':'s','m':'k',
                                              'k':'m', 'h':'d', 'd':'h', 'b':'v', 'v':'b', 'n':'n'}
# Description: case preserving reverse complementation of nucleotide sequences
def reverseComplement(sequence):
        return "".join([complement_alphabet.get(nucleotide, '') for nucleotide in sequence[::-1]])

