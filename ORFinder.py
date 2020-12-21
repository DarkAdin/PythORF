# Dictionary for creating a complementary DNA strand
Compl = {
        "A":"T",
        "T":"A",
        "C":"G",
        "G":"C"
        }

# Dictionary for transcribing DNA to RNA
Transcription = {
        "A":"U",
        "C":"G",
        "G":"C",
        "T":"A"
        }

# The 'universal' genetic code, a dictionary for translating RNA triplets to aminoacids
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
        "CAA":"Q",
        "CAG":"Q",
        "AAU":"N",
        "AAC":"N",
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
        "GGG":"G"
        }

def ORFinder ( DNA ) :
    filename = input ( "Save ORF Finder results to file: ") + ".log"
    try :
        # Establish a starting position, otherwise it is defaulted to 0 (beginning position)
        Startposition = int ( input ( 'Desired start position (default: 0): ') )
        # If the user chooses a starting position outside from the actual range of the sequence or an invalid one, it is defaulted to 0, meaning the first position.
        if Startposition not in range ( len ( DNA ) - 1 ) :
            print('Start position out of range. Defaulting to 0.')
            Startposition = 0
    except :
        Startposition = 0
    # Opening the .log file to write results to it. No need to close the file since the 'with' argument closes it automatically. The file is opened in 'append' mode because it already exists and its contents should not be overwritten.
    with open ( filename, "w" ) as fileseq :
        try :
            # Try to establish an ending position, otherwise it is defaulted to the last possible position, the final position in the sequence
            print ( 'Default end position:' , len ( DNA ) )
            Endposition = int ( input ( 'Desired end position: ' ) )
            if Endposition not in range( len ( DNA ) - 1 ):
                print ( 'End position out of range. Defaulting to', len ( DNA ) - 1, '.' )
                Endposition = len ( DNA ) - 1
        except:
            Endposition = len ( DNA ) - 1
        # Both positions are written to the .log file.
        fileseq.write ( "\nStart position: " + str ( Startposition ) )
        fileseq.write ( "\nEnd position: " + str ( Endposition ) )
        # Summary of the selected positions and the length between them is printed on screen.
        print (
                'Selected DNA between positions:', Startposition,'-', Endposition,
                '. Selected DNA sequence is', len ( DNA [ Startposition:Endposition ] ), 'bp long.'
                )
        print ('75 pb means peptides at least 25 aminoacids long. Left blank: defaulted to 75 pb.' )
        # Allowing the user to specify the minimal length of an ORF. Open reading frames smaller than the specified value (in base pairs) are ignored in the analysis.
        try:
            ORFlength = int ( input ( 'Minimal length of ORF? (bp) ' ) )
        except:
            # If the user does not specify a value, the minimal length is defaulted to 75 bp.
            ORFlength = 75
        # The value is written to the .log file.
        fileseq.write ( '\nSelected minimal ORF length: ' + str ( ORFlength ) + ' (bp)' )
        # Set a counter of the total number of ORFs found.
        ORFs = 0
        # Counting to three to present the results of the three 'plus' frames
        for i in range( 3 ):
            # Since the range of three is actually counted from 0 to 2, adding 1 to 'i' in every loop results in the actual number of the current frame, which is printed on screen.
            print( 'Frame +', i + 1, ':' )
            # Writing the current frame to the .log file.
            fileseq.write ( "\nFrame +" + str ( i + 1 ) + ":\n" )
            # A list to contain the current codon. This variable will be emptied every three positions, see below.
            Codon = []
            # A list to contain the complete protein, translated directly from the RNA. This protein list is chopped into ORFs later, beginning from a Methionin ('M') and ending with a 'stop' (represented as '*').
            Protein = []
            # For every ribonucleotide in the RNA transcribed from the complementary DNA contained in the sequence file (and using the chosen starting and ending positions), which is read starting from the 'i' position (decided in the for loop)...
            for ribo in [ Transcription [ Nuc ]
                    for Nuc in [ Compl [ D ]
                                 for D in DNA [ Startposition : Endposition ] ] [ i : ] ] :
                # ... append said ribonucleotide to the Codon list.
                Codon.append ( ribo )
                # When the codon list reaches a length of exactly three ribonucleotides, translate it to an aminoacid through the dictionary, append it to the protein list declared above, and empty its content so it can be used again.
                if len ( Codon ) == 3:
                    # The codon is appended to the protein list in its translated joined form (string), since the dictionary translates strings with a length of three.
                    Protein.append ( Translation [ ''.join ( Codon ) ] )
                    # Codon list is just emptied, not re-declared. It is then available for the next loop.
                    del Codon [ : ]
            # Some ORFs take up the entire sequence, with no stop codons and just one start codon at the very beginning. These next 4 lines allow for detecting them easily.
            if 'M' in Protein and '*' not in Protein :
                ORFs += 1
                # Printing:
                    # - The current ORF number.
                    # - The number of nucleotides it contains (multiplying by 3 the number of aminoacids).
                    # - Its length in aminoacids.
                    # - The aminoacidic sequence itself in string form.
                print(
                        'ORF', ORFs, ':',
                        3 * ( len ( Protein [ Protein.index ( 'M' ) : ] ) ) + 3, 'nts |',
                        len ( Protein [ Protein.index ( 'M' ) : ] ), 'aas:',
                        ''.join ( Protein [ Protein.index ( 'M' ) : ] )
                        )
                # All of it is written to the .log file.
                fileseq.write (
                        'ORF' +  str ( ORFs ) +  ' : ' +
                        str ( 3 * ( len ( Protein [ Protein.index ( 'M' ) : ] ) ) + 3 ) +  'nts | ' +
                        str ( len ( Protein [ Protein.index ( 'M' ) : ] ) ) +  'aas: ' +
                        ''.join ( Protein [ Protein.index ( 'M' ) : ] )
                        )
                # And then, it is deleted from inside the protein list.
                #del Protein [ Protein.index ( 'M' ) : ]
                Protein [ Protein.index( 'M' ) ] = 'X'

            # This is the canonical loop, but it functions just like the one above. As long as the Protein list contains methionins and stops, it will be parsed for ORFs.
            while 'M' in Protein and '*' in Protein :
                # If there are no stops before any methionin, the path is clear and we can begin the search for ORFs.
                if Protein.index( 'M' ) < Protein.index ( '*' ) :
                    # Print only ORFs longer than the user-specified minimal length (or the default value of 75 base pairs).
                    if len ( Protein [ Protein.index ( 'M' ) : Protein.index ( '*' ) ] ) > int ( ORFlength ) / 3 :
                        ORFs += 1
                        print (
                                'ORF', ORFs, ':', 3 * ( len ( Protein [ Protein.index ( 'M' ) : Protein.index ( '*' ) ] ) ) + 3, 'nts |',
                                len ( Protein [Protein.index ( 'M' ) : Protein.index ( '*' ) ] ), 'aas:',
                                ''.join ( Protein [ Protein.index ( 'M' ) : Protein.index ( '*' ) ] )
                                )
                        # Write them to the .log file.
                        fileseq.write (
                                "\n" + 'ORF' +  str ( ORFs ) +  ': ' +  str ( 3 * ( len ( Protein [Protein.index ( 'M' ) : Protein.index ( '*' ) ] ) ) + 3 ) +  ' nts | ' +
                                str ( len ( Protein [ Protein.index ( 'M' ) : Protein.index ( '*' ) ] ) ) +  ' aas: ' +
                                ''.join ( Protein [ Protein.index ( 'M' ) : Protein.index ( '*' ) ] )
                                )
                        # And finally delete them from inside the protein list so they are not detected again in the next loop.
                    del Protein [ Protein.index ( 'M' ) : Protein.index ( '*' ) ]
                # If there are any stops before the first methionin, they are removed from the protein list first.
                else :
                    Protein.remove ( '*' )
        # Counting to three, the remaining three 'minus' frames. The rest of the process is identical to the 'plus' frames.
        for i in range ( 3 ) :
            print ( 'Frame -', i + 1, ':' )
            fileseq.write ( "\nFrame -" + str ( i + 1 ) + ":\n" )
            Codon = []
            Protein = []
            # Notice how this time we don't make a complementary DNA strand, since we are on 'minus' frames. Also, notice how the sequence is read backwards, so we don't need to declare a variable with the reversed sequence, we just parse it in the contrary direction.
            for ribo in [ Transcription [ Nuc ]
                    for Nuc in DNA [ Startposition : Endposition ] [ -i - 1: : -1 ] ] :
                Codon.append ( ribo )
                if len ( Codon ) == 3:
                    Protein.append ( Translation [ ''.join ( Codon ) ] )
                    del Codon [ : ]
            if 'M' in Protein and '*' not in Protein :
                ORFs += 1
                print (
                        'ORF', ORFs, ':', 3 * ( len ( Protein [ Protein.index ( 'M' ) : ] ) ) + 3,
                        'nts |', len ( Protein [ Protein.index( 'M' ) : ] ), 'aas:', ''.join ( Protein [ Protein.index ( 'M' ) : ] )
                        )
                del Protein [ Protein.index ( 'M' ) : ]
            while 'M' in Protein and '*' in Protein :
                if Protein.index ( 'M' ) < Protein.index ( '*' ) :
                    if len ( Protein [ Protein.index ( 'M' ) : Protein.index ( '*' ) ] ) > int ( ORFlength ) / 3:
                        ORFs += 1
                        print (
                                'ORF', ORFs, ':', 3 * ( len ( Protein [ Protein.index ( 'M' ) : Protein.index ( '*' ) ] ) ) + 3, 'nts |',
                                len ( Protein [ Protein.index ( 'M' ) : Protein.index ( '*' ) ] ), 'aas:',
                                ''.join ( Protein [ Protein.index ( 'M' ) : Protein.index ( '*' ) ] )
                                )
                        fileseq.write (
                                "\n" + 'ORF' +  str ( ORFs ) +  ': ' +  str ( 3 * ( len ( Protein [ Protein.index ( 'M' ) : Protein.index ( '*' ) ] ) ) + 3 ) +  ' nts | ' +
                                str ( len ( Protein [ Protein.index ( 'M' ) : Protein.index ( '*' ) ] ) ) +  ' aas: ' +
                                ''.join ( Protein [ Protein.index ( 'M' ) : Protein.index ( '*' ) ] )
                                )
                    del Protein [ Protein.index ( 'M' ) : Protein.index ( '*' ) ]
                else:
                    Protein.remove ( '*' )
        # Finally, print the total number of ORFs and write it to the .log file. The analysis is over.
        fileseq.write("\n" + str ( ORFs ) + " ORFs were found.\n")
        print ( ORFs, 'ORFs were found.' )
