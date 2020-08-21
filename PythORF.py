##########################
# MAIN ORF FINDER FUNCTION
##########################
def ORFinder(DNA):
    try:
        Startpos=int(input('Desired start position (default: 0): '))
        if Startpos not in range(len(DNA)-1):
            print('Start position out of range. Defaulting to 0.')
            Startpos = 0
    except:
        Startpos=0
    try:
        print('Default end position:', len(DNA)-1)
        Endpos=int(input('Desired end position: '))
        if Endpos not in range(len(DNA)-1):
            print('End position out of range. Defaulting to', len(DNA)-1,'.')
            Endpos=len(DNA)-1
    except:
        Endpos=len(DNA)-1
    print('Selected DNA between positions:', Startpos,'-', Endpos,'. Selected DNA sequence is', len(DNA[Startpos:Endpos]), 'bp long.')
    print('75 pb means peptides at least 25 aminoacids long. Default: 75 pb.')
    try:
        ORFlength=int(input('Minimal length of ORF? (pb) '))
    except:
        ORFlength = 75
    Compl = {
            "A":"T",
            "T":"A",
            "C":"G",
            "G":"C",
            }
    Transcription = {
            "A":"U",
            "C":"G",
            "G":"C",
            "T":"A",
            }
    Translation = {
            "UUU":"F",
            "UUC":"F",
            "UUA":"L",
            "UUG":"L",
            "CUU":"L",
            "CUC":"L",
            "CUA":"L",
            "CUG":"L",
            "AUU":"I",
            "AUC":"I",
            "AUA":"I",
            "AUG":"M",
            "GUU":"V",
            "GUC":"V",
            "GUA":"V",
            "GUG":"V",
            "UCU":"S",
            "UCC":"S",
            "UCA":"S",
            "UCG":"S",
            "CCU":"P",
            "CCC":"P",
            "CCA":"P",
            "CCG":"P",
            "ACU":"T",
            "ACC":"T",
            "ACA":"T",
            "ACG":"T",
            "GCU":"A",
            "GCC":"A",
            "GCA":"A",
            "GCG":"A",
            "UAU":"Y",
            "UAC":"Y",
            "UAA":"*",
            "UAG":"*",
            "UGA":"*",
            "CAU":"H",
            "CAC":"H",
            "CAA":"G",
            "CAG":"G",
            "AAU":"B",
            "AAC":"B",
            "AAA":"K",
            "AAG":"K",
            "GAU":"D",
            "GAC":"D",
            "GAA":"Q",
            "GAG":"Q",
            "UGU":"C",
            "UGC":"C",
            "UGG":"W",
            "CGU":"R",
            "CGC":"R",
            "CGA":"R",
            "CGG":"R",
            "AGA":"R",
            "AGG":"R",
            "AGU":"S",
            "AGC":"S",
            "GGU":"G",
            "GGC":"G",
            "GGA":"G",
            "GGG":"G",
            }
    ORFs=0
    for i in range(3):
        print('Frame +',i+1,':')
        Codon=[]
        Protein=[]
        for ribo in [Transcription[Nuc] for Nuc in [Compl[D] for D in DNA[Startpos:Endpos]][i:]]:
            Codon.append(ribo)
            if len(Codon)==3:
                Protein.append(Translation[''.join(Codon)])
                del Codon[:]
        # Some ORFs take up the entire sequence, with no stop codons and just one start codon. These next 4 lines allow for detecting them easily.
        if 'M' in Protein and '*' not in Protein:
            ORFs+=1
            print('ORF', ORFs, ':', 3*(len(Protein[Protein.index('M'):]))+3, 'nts |', len(Protein[Protein.index('M'):]), 'aas:', ''.join(Protein[Protein.index('M'):]))
            del Protein[Protein.index('M'):]
        while 'M' in Protein and '*' in Protein:
            if Protein.index('M') < Protein.index('*'):
                if len(Protein[Protein.index('M'):Protein.index('*')])>int(ORFlength)/3:
                    ORFs+=1
                    print('ORF', ORFs, ':', 3*(len(Protein[Protein.index('M'):Protein.index('*')]))+3, 'nts |', len(Protein[Protein.index('M'):Protein.index('*')]), 'aas:', ''.join(Protein[Protein.index('M'):Protein.index('*')]))
                del Protein[Protein.index('M'):Protein.index('*')]
            else:
                Protein.remove('*')
    for i in range(3):
        print('Frame -',i+1,':')
        Codon=[]
        Protein=[]
        for ribo in [Transcription[Nuc] for Nuc in DNA[Startpos:Endpos][-i-1::-1]]:
            Codon.append(ribo)
            if len(Codon)==3:
                Protein.append(Translation[''.join(Codon)])
                del Codon[:]
        if 'M' in Protein and '*' not in Protein:
            ORFs+=1
            print('ORF', ORFs, ':', 3*(len(Protein[Protein.index('M'):]))+3, 'nts |', len(Protein[Protein.index('M'):]), 'aas:', ''.join(Protein[Protein.index('M'):]))
            del Protein[Protein.index('M'):]
        while 'M' in Protein and '*' in Protein:
            if Protein.index('M') < Protein.index('*'):
                if len(Protein[Protein.index('M'):Protein.index('*')])>int(ORFlength)/3:
                    ORFs+=1
                    print('ORF', ORFs, ':', 3*(len(Protein[Protein.index('M'):Protein.index('*')]))+3, 'nts |', len(Protein[Protein.index('M'):Protein.index('*')]), 'aas:', ''.join(Protein[Protein.index('M'):Protein.index('*')]))
                del Protein[Protein.index('M'):Protein.index('*')]
            else:
                Protein.remove('*')
    print(ORFs, 'ORFs were found.')
    return 0
###########################
# TARGET-FINDING FUNCTION
###########################
def Targets(DNA):
    Count=0
    Targets=['ATG','CGGCGG','GCCGCC']
    for Target in Targets:
        print('Targets found for', Target,':')
        Ind=0
        Hits=0
        for Nucl in DNA:
            if Nucl == list(Target)[Count]:
                Count+=1
                Ind+=1
            else:
                Count-=Count
                Ind+=1
            if Count== len(list(Target)):
                Count-=Count
                Hits+=1
                print(Ind-len(list(Target))+1, '-', ''.join(Target), '-', Ind)
        print(Hits,''.join(Target), 'targets found.')
        print('--------------')
    return 0
######################################################
# INITIAL LOADING AND PROCESSING OF THE FILE SEQUENCE
######################################################
try:
    Sequence=list(open(input('Insert exact filename of DNA sequence to be read: '),'r').read())
except:
    print('Invalid file/name. Exiting.')
    exit(1)
if '>' == Sequence[0]:
    print("FASTA format detected. Stripping first line as name: ", ''.join(Sequence[Sequence.index('>'):Sequence.index('\n')]))
    del Sequence[Sequence.index('>'):Sequence.index('\n')]
for a in [' ', '\t', '\n']:
    while a in Sequence:
        Sequence.remove(a)
print('The DNA sequence is', len(Sequence), 'bp long.')
for a in Sequence:
    if a not in ['A','T','C','G']:
        print('The sequence contains invalid characters. Exiting.')
        exit(1)
ORFinder(Sequence)
if input('Find genomic targets (yes/no)? ') in ['No','no','NO','N','n']:
    exit(0)
Targets(Sequence)
