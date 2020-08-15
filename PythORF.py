############
# ORF FINDER
############
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
        Codon=[]
        Protein=[]
        for ribo in [Transcription[Nuc] for Nuc in [Compl[D] for D in DNA[Startpos:Endpos]][i:]]: #First, the DNA is turned into its complementary version. It is transcribed to RNA right after by passing each nucleotide through the dictionary for transcription. The DNA is indexed beginning from the 'i' position, being i 0, 1 or 2. After that indexing, which variates the starting position by one nucleotide each time, the DNA is passed through the Compl dictionary to obtain its complementary strand. After that, said complementary strand is passed through the Transcription dictionary to obtain the RNA. Perhaps in the future I might open this a little more for the sake of readability.
            Codon.append(ribo)
            if len(Codon)==3: #Making sure the codon is exactly 3 ribonucleotides long before translating it.
                Protein.append(Translation[''.join(Codon)]) #Every codon is translated right before appending its corresponding aminoacid to the nascent protein.
                del Codon[:] #Resetting the codon. The variable is not re-initialized, only emptied.
        while 'M' in Protein and '*' in Protein: #Making sure the protein contains methionins and stops. As long as this simple condition is met, ORFs will be scanned.
            if Protein.index('M') < Protein.index('*'):
                if len(Protein[Protein.index('M'):Protein.index('*')])>int(ORFlength)/3: #The chosen minimal ORF Length finally comes into play. The length of all ORFs must obey its exact value, or higher. Otherwise they are not considered.
                    ORFs+=1
                    print('ORF', ORFs, ':', 3*(len(Protein[Protein.index('M'):Protein.index('*')]))+3, 'nts |', len(Protein[Protein.index('M'):Protein.index('*')]), 'aas:', ''.join(Protein[Protein.index('M'):Protein.index('*')])) #The ORF itself! Not more than the two same indexes manipulated in different ways: the index for 'M' and the index for '*', meaning Methionins and stops.
                #Protein.remove('M') #Once the ORF has been shown, remove its beginning and end but not the section in between, which may contain more ORFs
                #Protein.remove('*')
                del Protein[Protein.index('M'):Protein.index('*')]
            else:
                Protein.remove('*') #Removal of stop codons that might come before any Methionin, interrupting the reading.
    for i in range(3):
        print('Frame -',i+1,':')
        Codon=[]
        Protein=[]
        for ribo in [Transcription[Nuc] for Nuc in DNA[Startpos:Endpos][-i-1::-1]]: #DNA is read backwards in all 3 remaining frames, and the rest of the process is exactly the same as the one above.
            Codon.append(ribo)
            if len(Codon)==3:
                Protein.append(Translation[''.join(Codon)])
                del Codon[:]
        while 'M' in Protein and '*' in Protein:
            if Protein.index('M') < Protein.index('*'):
                if len(Protein[Protein.index('M'):Protein.index('*')])>int(ORFlength)/3:
                    ORFs+=1
                    print('ORF', ORFs, ':', 3*(len(Protein[Protein.index('M'):Protein.index('*')]))+3, 'nts |', len(Protein[Protein.index('M'):Protein.index('*')]), 'aas:', ''.join(Protein[Protein.index('M'):Protein.index('*')]))
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
        print('Targets found for', Target,':')
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
    return 0

#####################################################
# PREPARING THE SEQUENCE AND CALLING OF THE FUNCTIONS
#####################################################
try:
    Sequence=list(open(input('Insert exact filename of DNA sequence to be read: '),'r').read()) #DNA is read from a specific text file. User is prompted to type the exact filename of the sequence, which must be in the same directory.
except: #If for some reason the file is invalid or its name is not correct, display an error message and exit.
    print('Invalid file/name. Exiting.')
    exit(1)
if '>' == Sequence[0]:
    print("FASTA format detected. Stripping first line as name:")
    print(''.join(Sequence[Sequence.index('>'):Sequence.index('\n')])) #The program assumes that the sequence's name is placed between the first '>' character and the next newline, typical of FASTA files. Prints it once on screen, deletes it from the sequence, and leaves the rest intact to be further processed.
    del Sequence[Sequence.index('>'):Sequence.index('\n')]
for a in [' ', '\t', '\n']: #List of typical unwanted characters
    while a in Sequence:
        Sequence.remove(a) #Cleaning up the sequence
print('The DNA sequence is', len(Sequence), 'bp long.')
for a in Sequence:
    if a not in ['A','T','C','G']: #After removing newlines, tabs, spaces and such, the sequence is scanned in search of invalid characters. If one is encountered, the program exits.
        print('The sequence contains invalid characters. Exiting.')
        exit(1)
ORFinder(Sequence) #Calling the main ORF Finder function on the sequence specified by the user.
if input('Find genomic targets (yes/no)? ') in ['No','no','NO','N','n']: #Allowing exiting without finding targets
    exit(0)
Targets(Sequence)
#if input('Find genomic alignments (yes/no)? ') in ['No','no','NO','N','n']: #Allowing exiting without finding targets
#    exit(0)
#Alignment(Sequence)
