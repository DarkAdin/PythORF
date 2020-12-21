from ORFinder import ORFinder
from Targets import Targets

try :
    Sequence = list ( open ( input ( 'Insert exact filename of DNA sequence to be read: ' ), 'r' ).read() )
except :
    print( 'Invalid file/name. Exiting.' )
    exit ( 1 )
# If the sequence contains a '>' character at the very beginning, it means it is a FASTA file and the first line should be treated as the name of the sequence, ignoring it in the analysis.
if '>' == Sequence [ 0 ] :
    print (
            "FASTA format detected. Stripping first line as name: ",
            ''.join ( Sequence [ Sequence.index ( '>' ) : Sequence.index ( '\n' ) ] )
            )
    del Sequence [ Sequence.index ( '>' ) : Sequence.index ( '\n' ) ]
for a in [ ' ', '\t', '\n' ] :
    while a in Sequence :
        Sequence.remove ( a )
    else :
        pass
print ( 'The DNA sequence is', len ( Sequence ), 'bp long.' )
# If the sequence contains characters other than A, T, C or G, exit.
for a in Sequence :
    if a not in [ 'A','T','C','G' ] :
        print ( '\nThe sequence contains invalid characters. Exiting.' )
        exit ( 1 )

# Calling the ORF Finder function on the sequence.
ORFinder ( Sequence )
if input ( 'Find genomic targets (yes/no)? ' ) in [ 'No','no','NO','N','n' ] :
    exit ( 0 )
Targets ( Sequence )
