# Progress Log

**April 26th, 4-5pm**: I intended to setup the project itself (python environment) and download the 20 *E. coli* genomes. I accomplished this.

**April 26th, 5:30-6:30pm**: I intended to write the code to preprocess all the datasets into one JSON file. I accomplished this.

**April 26th, 8:00-9:00pm**: I intended to start work on the kmer filtering code to limit the number of pairs that I run Smith-Waterman on. I got an implementation but I realized the big-O complexity of it was untenable for this size of dataset.

**April 26th, 9:00-10:00pm**: I intended to figure out a way to improve the kmer filter implementation. I found a way to do many less comparisons by using a dictionary and flipping from iterating over proteins to iterating over kmers. This makes the algorithm run in ~1 minute.

**May 8th, 4:00-5:00pm**: My goal was to implement Smith-Waterman on the pairwise comparisons. I accomplished this, although it was naively way too slow, so I converted to Numba and now can run in ~1 hour.
