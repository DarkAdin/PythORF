# Sequencetool
A little, simple, humble Python algorithm that finds all possible Open Reading Frames within a given DNA sequence, and finds some arbitrary targets.

This algorithm tries to emulate the results and reliability of the ORF Finder from NCBI, while avoiding the use of any Python module to keep it clean and self-reliant. It should be noticeably fast on any regular sequence and the CPU shouldn't even sweat, although it'd be advisable to avoid very large sequences (more than 100.000 bp).

## Functionalities

### Transcription and translation
Using Python's built-in method of creating dictionaries, the algorithm relies on two of them to fully transcribe and translate the sequence: one for transcription, and one for translation (which includes the full 'universal' genetic code).

Using these two dictionaries and a very compacted nested 'for' loop, two mRNAs are generated on the fly, one for 'plus', one for 'minus' frames, and translated into as many peptides as the current frame contains. This tries to avoid explicit list declaration, using as little memory as possible. In Python, this should go always first.

After finding a peptide that's long enough for the open reading frame's length (which can be set by the user), it is printed on the screen, and the loop jumps to the next.

The loop goes on until all 6 frames are finished. Empty frames are also printed on screen.

### Genomic targets
Since it's a common tool in genomics, there's also an option to find every possible genomic target on the DNA sequence. The entire DNA sequence is scanned and every target printed on screen in another 'for' loop.

## To do
This little algorithm is far from finished, far from polished. Things to be added:
  -Make the user decide as many parameters as possible, for example whether to just find open reading frames or just find genomic targets, without having to do one, then the other.
  -Include a system to rate every possible peptide and present the most probable one, given the context of translation start and ending sites.
  -Include a tool to design primers for the amplification of the sequence if possible. Again, the user here should have the last word on everything.
  -Make cleaner, even more direct code. Focusing on memory usage, the algorithm should be as efficient as possible, without a single variable, a single list, declared in vain. In Python, this should be one of our main concerns.
  -And more...
