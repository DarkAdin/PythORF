# PythORF
A little Python algorithm that finds all possible Open Reading Frames within a given DNA sequence, and finds some arbitrary targets.

This algorithm tries to emulate the results and reliability of the ORF Finder from NCBI, while avoiding the use of any additional Python module to keep it clean and self-reliant. It should be noticeably fast on any regular coding sequence and the CPU shouldn't even sweat, although it'd be advisable to avoid very large sequences (more than 100 Kb).

## Functionalities

### Reading of the file
The user should enter the filename, which must be present in the same directory. FASTA is not supported yet.

### Transcription and translation
Using Python's built-in method of creating dictionaries, the algorithm relies on three of them to fully transcribe and translate the sequence: one for transcription, one for translation (which includes the full "universal" genetic code), and the third one serves the only purpose of generating a complementary DNA.

Using these three dictionaries and two very compacted nested *for* loops, 6 mRNAs are generated on the fly, three for *plus*, three for *minus* frames, and translated into as many peptides as the current frame contains. This tries to avoid explicit list declaration, using as little memory as possible.

After finding a peptide that's long enough for the open reading frame's length (which can be set by the user in the beginning), it is printed on the screen, and the loop jumps to the next.

The loop goes on until all 6 frames are finished. Empty frames are not ignored and are also presented on screen.

### Genomic targets
Since it's a common tool in genomics, there's also an option to find every possible genomic target on the DNA sequence from an easily modifiable arbitrary list. The entire DNA sequence is scanned and every target printed on screen in another 'for' loop along with its exact position within it.

### Sample sequences
The repository includes two sample sequences in *txt* format:
* Human Interferon - gamma (IFNG).
* Natriuretic Peptide A (NPPA).

## Purpose
The main goal of this algorithm is to bring some of the most common tools both to unexperienced and experienced students and academics in the field of Molecular Biology and Genomics, without having to depend on some website, third party organization, licensed application or bloated libraries. The massive use of licensed applications plagues the entire academic field and personally is something I'd have loved to avoid. Everyone knows that most software these days suffers from a serious lack of art in their design.
Also, Python is a great language.

The user should actively modify the code to suit their needs. In fact, there's another functionality which may be activated by this method: the search for *nested* ORFs. This can be done by uncommenting one pair of lines in the code and commenting another in both *for* loops within the ORF Finder function. There are enough clarifications in the form of comments in the code itself, so they should be easy to find.

## Usage
Just clone this repository or just download the Python code itself and execute it on any sequence you want. It will guide you through all the steps and prompt to enter useful information. If left blank, the values are defaulted. The only thing needed is a text file which contains a DNA sequence with only the four possible characters A,G,C,T in it. FASTA format is still not supported. If the file contains tabs, spaces or newlines, they are all deleted, so the file does not need to have just one line and they are not a problem.

## To do
This little tool is far from finished, far from polished. Things to be added:
* Make the user decide as many parameters as possible.
* Allow starting from RNA. There are plenty of cases of ncRNAs that may contain little ORFs, some actually functional.
* Include a method of reverse transcripting said RNA.
* Include a system to rate every possible peptide and present the most probable ones, given the context of translation start and ending sites at the mRNA. Perhaps include Shine-Dalgarno sites, and options for alternative poli-adenilation and termination.
* Include a tool to design primers for the amplification of the sequence if possible. Again, the user here should have the last word on deciding the primers best suited for their needs.
* Make cleaner, even more direct code. Focusing on minimal memory usage, the algorithm should be as efficient as possible, without a single variable, a single list, declared in vain. Every new tool should be prompted to be used in the form of a function and not entered by default; the user should query the tasks first, and their order should be easily set. This makes the program modular and introduces intentional pauses so the CPU usage is not as intensive.
* Allow writing the output into a log file to be further processed.
* Perhaps make it Ncurses-compatible, allowing for a pseudo-graphical interface within the terminal.
* And more...
In order to add useful functionalities, this program could benefit from using certain, built-in Python libraries like **os** or **system**. This possibility will be explored in the future, always keeping in mind that simpler is better.
