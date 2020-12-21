def Targets ( DNA ) :
    # Opening the .log file in append mode, since it already contains previous results.
    filename = input ( "Save target-finding results to file: ") + ".log"
    with open ( filename, "a" ) as fileseq :
        Count = 0
        Targets = []
        print ( "Press Enter to skip introducing more targets.\n" )
        while True:
            trgt = input ( "Search for a specific target: ")
            if len ( trgt ) < 1:
                break
            else:
                Targets.append ( trgt )
                continue
        for Target in Targets:
            fileseq.write ( "Targets found for " + Target + ":\n")
            print ( 'Targets found for', Target, ':' )
            Ind = 0
            Hits = 0
            for Nucl in DNA :
                if Nucl == list ( Target ) [ Count ] :
                    Count += 1
                    Ind += 1
                else:
                    Count -= Count
                    Ind += 1
                if Count == len ( list ( Target ) ) :
                    Count -= Count
                    Hits += 1
                    print (
                            Ind - len ( list ( Target ) ) + 1, '-',
                            ''.join ( Target ), '-', Ind
                            )
                    fileseq.write (
                            str ( Ind - len ( list ( Target ) ) + 1 ) + ' - ' +
                            ''.join ( Target ) + ' - ' + str ( Ind ) + "\n"
                            )
            fileseq.write (
                    str ( Hits ) + " " +
                    ''.join ( Target ) + " targets found.\n"
                    )
            fileseq.write( "-------------\n" )
            print( Hits,''.join ( Target ), ' targets found.' )
            print( '--------------' )
