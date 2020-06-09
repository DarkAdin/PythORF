#############
# ORF FINDER
#############
def ORFinder(DNA):
    try:
        Startpos=int(input('Desired start position (default: 0): ')) #User is prompted to input start and end positions. Incorrect values are defaulted to the actual start and end of the complete sequence.
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
    print('Selected DNA between positions:', Startpos,'-', Endpos,'. Selected DNA sequence is', len(DNA[Startpos:Endpos]), 'bp long.') #Whether the user inputs correct values or not, the sequence is set to be read in the decided range.
    print('75 pb means peptides at least 25 aminoacids long (default: 75 bp).')
    try:
        ORFlength=int(input('Minimal length of ORF? (pb) '))
    except:
        ORFlength = 75
    Compl = {
            "A":"T",
            "T":"A",
            "C":"G",
            "G":"C",
            } #Dictionary to create a complementary DNA strand
    Transcription = {
            "A":"U",
            "C":"G",
            "G":"C",
            "T":"A",
            } #Dictionary to transcribe the DNA
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
            } #Dictionary of triplets that contains the 'universal' genetic code to help translate every RNA into its corresponding protein
    ORFs=0 #Total number of ORFs counted
    for i in range(3):
        print('Frame +',i+1,':')
        Count=0
        Codon=[]
        Protein=[]
        for ribo in [Transcription[Nuc] for Nuc in [Compl[D] for D in DNA[Startpos:Endpos]][i:]]: #First, the DNA is turned into its complementary version. It is transcribed to RNA right after by passing each ribonucleotide through the dictionary for transcription.
            Count+=1
            Codon.append(ribo)
            if Count==3: #Making sure the codon is exactly 3 ribonucleotides long before translating it.
                Count-=3 #Resetting the counter. I feel unconfortable about these counters, but for now I can't seem to find an alternative...
                Protein.append(Translation[''.join(Codon)]) #Every codon is translated right before appending its corresponding aminoacid to the nascent protein.
                del Codon[:] #Resetting the codon. The variable is not re-initialized, only emptied.
        while 'M' in Protein and '*' in Protein: #Making sure the protein contains methionins and stops. As long as this simple condition is met, ORFs will be scanned.
            if Protein.index('M') < Protein.index('*'):
                if len(Protein[Protein.index('M'):Protein.index('*')])>int(ORFlength)/3: #The chosen minimal ORF Length finally comes into play. The length of all ORFs must obey its exact value, or higher. Otherwise they are not considered.
                    print(3*(len(Protein[Protein.index('M'):Protein.index('*')]))+3, 'nts |', len(Protein[Protein.index('M'):Protein.index('*')]), 'aas:', ''.join(Protein[Protein.index('M'):Protein.index('*')])) #The ORF itself! Not more than the two same indexes manipulated in different ways: the index for 'M' and the index for '*', meaning Methionins and stops.
                    ORFs+=1
                #Protein.remove('M') #Once the ORF has been shown, remove its beginning and end but not the section in between, which may contain more ORFs ('nested').
                #Protein.remove('*')
                del Protein[Protein.index('M'):Protein.index('*')] #This setting deactivates the search for nested ORFs. Comment it and uncomment the two lines above to obtain nested ORFs. The search for nested ORFs is turned off by default in NCBI's ORF Finder as well. Activating it means obtaining a greater number of ORFs than otherwise. Remember to change this setting on the next 'for' loop as well so the setting is fully active.
            else:
                Protein.remove('*') #Removal of stop codons that might come before any Methionin, interrupting the reading.
    for i in range(3):
        print('Frame -',i+1,':')
        Count=0
        Codon=[]
        Protein=[]
        for ribo in [Transcription[Nuc] for Nuc in DNA[Startpos:Endpos][-i-1::-1]]: #DNA is read backwards in all 3 remaining frames, and the rest of the process is exactly the same as the one above.
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
                    ORFs+=1
                #Protein.remove('M')
                #Protein.remove('*')
                del Protein[Protein.index('M'):Protein.index('*')]
            else:
                Protein.remove('*')
    print(ORFs, 'ORFs were found.')
    return 0
#########
# TARGETS
#########
def Targets(DNA):
    Count=0
    Targets=['ATG','CGGCGG','GCCGCC'] #List of targets to search for. It would be nice to make the user input all desired targets and do the search once they are done.
    for Target in Targets:
        print('Targets found for ', ''.join(Target),':')
        Ind=0 #Represents the index of the last position
        Hits=0
        for Nucl in DNA:
            if Nucl == list(Target)[Count]: #If the nucleotide in the sequence matches a nucleotide in the target. Each step in the loop checks the next nucleotide both in the sequence and the target.
                Count+=1
                Ind+=1
            else: #No matter if the previous nucleotides match, if one of them doesn't, then the homology is not complete and the target should not be considered.
                Count-=Count
                Ind+=1
            if Count== len(list(Target)):
                Count-=Count
                Hits+=1
                print(Ind-len(list(Target)), '-', ''.join(Target), '-', Ind)
        print(Hits,''.join(Target), 'targets found.')
        print('--------------')

########################
# PREPARING THE SEQUENCE
########################
try:
    Sequence=list(open(input('Insert exact filename of DNA sequence to be read: '),'r').read()) #DNA is read from a specific file. User is prompted to type the exact filename of the sequence, which must be in the same directory. FASTA format is not supported yet, since there needs to be a way of dealing with the first line that contains the metadata of the sequence. It will be included in time!
except: #If for some reason the file is invalid or its name is not correct, display an error message and exit.
    print('Invalid file. Exiting.')
    exit()
for a in [' ', '\t', '\n', '>']: #List of typical unwanted characters
    while a in Sequence:
        Sequence.remove(a) #Cleaning up the sequence
print('The DNA sequence is', len(Sequence), 'bp long.')
for a in Sequence:
    if a not in ['A','T','C','G']: #After removing newlines, tabs, spaces and such, the sequence is scanned in search of invalid characters. If so, the program exits.
        print('The sequence contains invalid characters. Exiting.')
        exit()
ORFinder(Sequence) #Calling the main ORF Finder function on the specified sequence.
if input('Find genomic targets (yes/no)? ') in ['No','no','NO','N','n']: #Allowing exiting without finding targets
    exit()
Targets(Sequence)
