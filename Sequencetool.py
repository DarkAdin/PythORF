DNA=list(open('IFNG.fas','r').read())
for a in [' ', '\t', '\n']:
    while a in DNA:
        DNA.remove(a) #Cleaning up the sequence
print('75 pb means peptides at least 25 aminoacids long')
ORFlength=input('Minimal length of ORF? (pb) ')
ORFlength=75
#Dictionaries for: creating a complementary DNA; transcribing the DNA; translating the RNA
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
for i in range(3):
    print('Frame +',i+1,':')
    Count=0
    Codon=[]
    Protein=[]
    for ribo in [Transcription[Nuc] for Nuc in [Compl[D] for D in DNA][i:]]:
        Count+=1
        Codon.append(ribo)
        if Count==3:
            Count-=3
            Protein.append(Translation[''.join(Codon)])
            del Codon[:]
    while 'M' in Protein and '*' in Protein:
        if Protein.index('M') < Protein.index('*'):
            if len(Protein[Protein.index('M'):Protein.index('*')])>int(ORFlength)/3:
                print(3*(len(Protein[Protein.index('M'):Protein.index('*')]))+3, 'nts |', len(Protein[Protein.index('M'):Protein.index('*')]), 'aas:', ''.join(Protein[Protein.index('M'):Protein.index('*')]))
            Protein.remove('M')
            Protein.remove('*')
        else:
            Protein.remove('*')
for i in range(3):
    print('Frame -',i+1,':')
    Count=0
    Codon=[]
    Protein=[]
    for ribo in [Transcription[Nuc] for Nuc in DNA[-i-1::-1]]:
        Count+=1
        Codon.append(ribo)
        if Count==3:
            Count-=3
            Protein.append(Translation[''.join(Codon)])
            del Codon[:]
    while 'M' in Protein and '*' in Protein:
        if Protein.index('M') < Protein.index('*'):
            if len(Protein[Protein.index('M'):Protein.index('*')])>int(ORFlength)/3:
                print(3*(len(Protein[Protein.index('M'):Protein.index('*')]))+3, 'nts |', len(Protein[Protein.index('M'):Protein.index('*')]), 'aas:', ''.join(Protein[Protein.index('M'):Protein.index('*')]))
            Protein.remove('M')
            Protein.remove('*')
        else:
            Protein.remove('*')
input('Press any key to find genomic targets...')
#########
# TARGETS
#########
Contador=0
Dianas_DNA=['ATG','CGGCGG','GCCGCC']
for Diana in Dianas_DNA:
    print('Targets found for ', ''.join(Diana),':')
    Ind=0
    Blancos=0
    for Nucl in DNA:
        if Nucl == list(Diana)[Contador]:
            Contador+=1
            Ind+=1
        else:
            Contador=0
            Ind+=1
        if Contador == len(list(Diana)):
            Contador=0
            Blancos+=1
            print(Ind-len(list(Diana)), '-', ''.join(Diana), '-', Ind)
    print(Blancos, 'targets', ''.join(Diana), 'found.')
    print('--------------')
