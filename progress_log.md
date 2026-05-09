# Progress Log

**April 26th, 4:00-5:00pm**: I intended to setup the project itself (python environment) and download the 20 *E. coli* genomes. I accomplished this.

**April 26th, 5:30-6:30pm**: I intended to write the code to preprocess all the datasets into one JSON file. I accomplished this.

**April 26th, 8:00-9:00pm**: I intended to start work on the kmer filtering code to limit the number of pairs that I run Smith-Waterman on. I got an implementation but I realized the big-O complexity of it was untenable for this size of dataset.

**April 26th, 9:00-10:00pm**: I intended to figure out a way to improve the kmer filter implementation. I found a way to do many less comparisons by using a dictionary and flipping from iterating over proteins to iterating over kmers. This makes the algorithm run in ~1 minute.

**May 8th, 4:00-5:00pm**: My goal was to implement Smith-Waterman on the pairwise comparisons. I accomplished this, although it was naively way too slow, so I converted to Numba and now can run in ~1 hour.

**May 8th, 5:00-6:00pm**: My goal was to get as far as possible in an implementation of MCL on the scored protein pairs, while at the same time letting the protein pair scoring run. I ran the scoring but realized I want to do normalization in stage 1 so I'm re-running it. I got most of MCL implemented.

**May 8th, 6:00-7:00pm**: My goal was to finish MCL implementation, run it, and validate the results. I finished implementing it but memory usage exploded beyond what my MacBook could handle.

**May 8th, 7:00-8:00pm**: My goal was to find a solution to the memory usage problem, run MCL, and look at the results. I solved it by pruning edges below a minimum threshold weight, so we drop from 11 million to ~1-2 million which is more managable. Running clustering was quick and the results seem biologically accurate.

**May 9th, 11:00am-12:00pm**: My goal was to get a cool visualization of the graph. I did this by exporting the nodes and edges as CSVs, loading them into Gephi, and playing around with layout and coloring techniques. I settled on sharing 10 colors between the MCL clusters and using OpenOrd to layout the graph. Layout takes ~5 minutes and gets good separation of clusters.

**May 9th, 2:00pm-8:00pm (6 hours)**: My goal was to make the poster. I did this!
