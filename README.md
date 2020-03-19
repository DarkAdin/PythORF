# Sequencetool
A little, simple, humble Python algorithm that finds all possible Open Reading Frames from within a given sequence and finds some arbitrary targets. Work still in progress!

This algorithm tries to emulate the results and reliability of the ORF Finder from NCBI, while avoiding the use of as many Python modules/libraries as possible to achieve it. This, although making the algorithm simple, has some obvious effect on its further capabilities.

## Functionalities

### Transcription and translation

Using Python's built-in method of creating dictionaries, the algorithm relies on two of them to fully transcribe and translate a sequence.

6 different mRNAs are generated, in a big loop, one by one, and on the fly, from just one DNA sequence.

While the transcription is pretty straightforward, the translation to protein is not. Every RNA sequence is entirely scanned in all 6 different reading frames: 3 in forward (plus, '+') and 3 in reverse (minus, '-').
For every one of these 6 iterations, the algorithm tries to find all possible start and stop codons from within the mRNA sequence, prints the result with the corresponding positions at its sides, and jumps to the next.
