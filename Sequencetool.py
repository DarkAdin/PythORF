try:
    DNA=list(open(input('Insert exact filename of DNA sequence to be read: '),'r').read()) #DNA is read from a specific file. User is prompted to type the exact filename of the sequence. FASTA format is not supported yet, since there needs to be a way of dealing with the first line that contains the metadata of the sequence. It will be implemented in time.
except: #If for some reason the file is invalid, doesn't exist or the path to it is not correct, display an error and exit.
    print('Invalid file/path. Exiting.')
    exit()
for a in [' ', '\t', '\n', '>']: #List of unwanted characters
    while a in DNA:
        DNA.remove(a) #Cleaning up the sequence
print('The DNA sequence is', len(DNA), 'bp long.')
Startpos=input('Start position (or leave blank): ') #The user is prompted to decide, if they want to, the start and end positions within the DNA sequence. The default option, of course, is scanning the entire DNA sequence.
Endpos=input('End position (or leave blank): ')
try:
    int(Startpos)
except:
    Startpos=0
try:
    int(Endpos)
except:
    Endpos=len(DNA)-1
if int(Startpos) not in range(len(DNA)): #If the user leaves them blank or inputs an invalid value outside the range of the DNA itself, they are defaulted to the start and end of the molecule, respectively.
    Startpos=0
if Endpos < Startpos or Endpos not in range(len(DNA)-1):
    Endpos=len(DNA)-1
DNA=DNA[int(Startpos):int(Endpos)] #Whether the user inputs correct values or not, the sequence is set to be read in the decided range.
print('Selected DNA between positions:', Startpos,'-', Endpos, '. Selected DNA sequence now is', len(DNA), 'bp long.')
print('75 pb means peptides at least 25 aminoacids long (default: 75 bp).')
ORFlength=input('Minimal length of ORF? (pb) ')
if ORFlength is not int: #If the user inserts anything other than an int, ORF length is defaulted to 75 bp
        ORFlength=75
#Dictionaries for: creating a complementary DNA; transcribing the DNA; translating the RNA
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
        } #Dictionary that contains the 'universal' genetic code to help translate every RNA into its corresponding protein
ORFs=0
for i in range(3):
    print('Frame +',i+1,':')
    Count=0
    Codon=[]
    Protein=[]
    for ribo in [Transcription[Nuc] for Nuc in [Compl[D] for D in DNA][i:]]: #First, the DNA is turned into its complementary version. It is transcribed to RNA right after
        Count+=1
        Codon.append(ribo)
        if Count==3: #Ensuring that the codon is exactly 3 ribonucleotids long
            Count-=3 #Resetting the counter. I feel unconfortable about these counters, but for now I can't seem to find an alternative...
            Protein.append(Translation[''.join(Codon)]) #Every codon is translated right before appending its corresponding aminoacid to the nascent protein
            del Codon[:] #Resetting the codon
    while 'M' in Protein and '*' in Protein: #Making sure the protein contains methionins and stops.
        if Protein.index('M') < Protein.index('*'):
            if len(Protein[Protein.index('M'):Protein.index('*')])>int(ORFlength)/3: #The chosen minimal ORF Length comes into play. The length of all ORFs must obey its exact value, or higher.
                print(3*(len(Protein[Protein.index('M'):Protein.index('*')]))+3, 'nts |', len(Protein[Protein.index('M'):Protein.index('*')]), 'aas:', ''.join(Protein[Protein.index('M'):Protein.index('*')])) #The ORF itself! Not more than the two same indexes manipulated in different ways: the index for 'M' and the index for '*', and nothing else.
                ORFs+=1
            #Protein.remove('M') #Once the ORF has been shown, remove its beginning and end but not the section in between, which may contain more ORFs ('nested').
            #Protein.remove('*')
            del Protein[Protein.index('M'):Protein.index('*')] #This setting deactivates the search for nested ORFs. Comment it and uncomment the two lines above to obtain nested ORFs. The search for nested ORFs is turned off by default in NCBI's ORF Finder as well. Activating it means obtaining a greater number of ORFs.
        else:
            Protein.remove('*') #Removal of stop codons that might come before any Methionin, interrupting the reading of any ORF.
for i in range(3):
    print('Frame -',i+1,':')
    Count=0
    Codon=[]
    Protein=[]
    for ribo in [Transcription[Nuc] for Nuc in DNA[-i-1::-1]]: #DNA is read backwards in all 3 remaining frames, and the rest of the process is exactly the same as above.
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
print(ORFs, 'ORFs found.')
input('Press any key to find genomic targets... ')
#########
# TARGETS
#########
def Targets(DNA): #Putting all this in a function to allow the user calling it in future versions, since it is an optional tool and should only prompt when the user chooses.
    Count=0
    Targets=['ATG','CGGCGG','GCCGCC']
    for Target in Targets:
        print('Targets found for ', ''.join(Target),':')
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
                print(Ind-len(list(Target)), '-', ''.join(Target), '-', Ind)
        print(Hits,''.join(Target), 'targets found.')
        print('--------------')
Targets(DNA)
