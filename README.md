# Sequencetool
A little Python algorithm that finds all possible Open Reading Frames within a given DNA sequence, and finds some arbitrary targets.

This algorithm tries to emulate the results and reliability of the ORF Finder from NCBI, while avoiding the use of any additional Python module to keep it clean and self-reliant. It should be noticeably fast on any regular coding sequence and the CPU shouldn't even sweat, although it'd be advisable to avoid very large sequences (more than 100 Kb).

## Functionalities

### Transcription and translation
Using Python's built-in method of creating dictionaries, the algorithm relies on three of them to fully transcribe and translate the sequence: one for transcription, one for translation (which includes the full 'universal' genetic code). The third one serves the only purpose of generating a complementary DNA for the *plus* frames.

Using these three dictionaries and two very compacted nested 'for' loops, 6 mRNAs are generated on the fly, one for *plus*, one for *minus* frames, and translated into as many peptides as the current frame contains. This tries to avoid explicit list declaration, using as little memory as possible.

After finding a peptide that's long enough for the open reading frame's length (which can be set by the user), it is printed on the screen, and the loop jumps to the next.

The loop goes on until all 6 frames are finished. Empty frames are also printed on screen.

### Genomic targets
Since it's a common tool in genomics, there's also an option to find every possible genomic target on the DNA sequence. The entire DNA sequence is scanned and every target printed on screen in another 'for' loop along with its exact position.

## Purpose
The main purpose of this algorithm is to bring the most common tools both to unexperienced and experienced students and academics in the field of Molecular Biology and Genomics, without having to depend on some website, third party organization, or licensed application. The massive use of licensed applications plagues the entire academic field and personally is something I'd have loved to avoid.
Also, Python is a great language.

## To do
This little tool is far from finished, far from polished. Things to be added:
* Make the user decide as many parameters as possible, for example whether to just find open reading frames or just find genomic targets, without having to do one, then the other.
* Allow starting from RNA. Also, there are plenty of cases of ncRNAs that may contain little ORFs, some actually functional.
* Include a method of reverse transcripting said RNA.
* Include a system to rate every possible peptide and present the most probable one, given the context of translation start and ending sites at the mRNA. Perhaps include Shine-Dalgarno sites, and options for alternative poli-adenilation and termination.
* Include a tool to design primers for the amplification of the sequence if possible. Again, the user here should have the last word on deciding the primers best suited for their needs.
* Make cleaner, even more direct code. Focusing on memory usage, the algorithm should be as efficient as possible, without a single variable, a single list, declared in vain. When writing in Python, this always has to be one of the main concerns.
* And more...
