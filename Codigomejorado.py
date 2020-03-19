DNA=list(open('sequence2.txt', 'r', encoding='utf-8').read())
while '\n' in DNA:
    DNA.remove('\n')
while '\t' in DNA:
    DNA.remove('\t')
while ' ' in DNA:
    DNA.remove(' ')
print("(75 bp means peptides at least 25 aminoacids long)")
ORFlong=input('Minimal length of ORF? (bp) ')
Transcripcion = {
        "A":"U",
        "C":"G",
        "G":"C",
        "T":"A",
        }
ARN_aa = {
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
PlusIndex=0
MinusIndex=0
Dianas_DNA=['TATA', 'ATAT', 'TCCTCC', 'ATG', 'CCCTC', 'CCCCCC']
for Nuc in DNA:
    while PlusIndex < 3: #Indexing through all 3 'forward' frames
        print('Frame +', PlusIndex+1, ':')
        Contador=0
        Codon=[]
        Proteina=[]
        for ribo in [Transcripcion[Nuc] for Nuc in DNA[PlusIndex:]]:
            Contador +=1
            Codon.append(ribo)
            if Contador == 3:
                Contador = 0
                Proteina.append(ARN_aa[''.join(Codon)])
                Codon=[]
        while 'M' in Proteina and '*' in Proteina:
            try:
                if Proteina.index('M') < Proteina.index('*'):
                    if len(Proteina[Proteina.index('M'):Proteina.index('*')]) > int(ORFlong)/3:
                        print(3*(len(Proteina[Proteina.index('M'):Proteina.index('*')])), 'nts', '|', len(Proteina[Proteina.index('M'):Proteina.index('*')]), 'aas:',''.join(Proteina[Proteina.index('M'):Proteina.index('*')]))
                    del Proteina[Proteina.index('M'):Proteina.index('*')]
                if Proteina.index('M') > Proteina.index('*'):
                    Proteina.remove('*')
            except:
                pass
        PlusIndex+=1
    else:
        while MinusIndex < 3: #Indexin through 'reverse' frames
            print('Frame -', MinusIndex+1, ':')
            Contador=0
            Codon=[]
            Proteina=[]
            for ribo in [Transcripcion[Nuc] for Nuc in DNA[-MinusIndex-1::-1]]:
                Contador += 1
                Codon.append(ribo)
                if Contador==3:
                    Contador=0
                    Proteina.append(ARN_aa[''.join(Codon)])
                    Codon=[]
            while 'M' in Proteina and '*' in Proteina:
                try:
                    if Proteina.index('M') < Proteina.index('*'):
                        if len(Proteina[Proteina.index('M'):Proteina.index('*')]) > int(ORFlong)/3:
                            print(3*(len(Proteina[Proteina.index('M'):Proteina.index('*')])), 'nts', '|', len(Proteina[Proteina.index('M'):Proteina.index('*')]), 'aas:',''.join(Proteina[Proteina.index('M'):Proteina.index('*')]))
                        del Proteina[Proteina.index('M'):Proteina.index('*')]
                    if Proteina.index('M') > Proteina.index('*'):
                        Proteina.remove('*')
                except:
                    pass
            MinusIndex+=1
        else:
            break
input('Press any key to find genomic targets...')
#########
# DIANAS
#########
Contador=0
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
#input('Pulse Enter para hacer un BLAST...')
#############
# BLASTER: EN DESARROLLO
#############
#Contador=0
#for Nuc in DNA:
#    if Nuc == ADN2[Contador]:
#
#    Contador+=1
