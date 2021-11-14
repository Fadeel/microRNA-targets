# microRNA-targets
This Python script is used to find potential miRNA seed binding sites. Basically, this script is used to find a small RNA sequence pattern (6-8 characters long) in a database of longer RNA seqences. For fast searching, longer RNA sequences are store in a hash table. 

## input arguments:
It takes 3 parameter: target sequences (Fasta file), miRNA sequences (Fasta file) and seed length(6,7 and 8).


### target sequences/miRNA files should be in FASTA format as follows
```
>sequence_id 
ATCTGTTTGTCAA
```

"id" will be used in the output files.


### seed length

6 => sequence starting at position 2 (from 5' side of miRNA seq.)  to position 7 

7 => sequence starting at position 2 (from 5' side of miRNA seq.)  to position 8  *** Recommended***

8 => sequence starting at position 1 (from 5' side of miRNA seq.)  to position 8


## Output files

1. `7mers.txt` (If you choose seed length = 7, each row represents a miRNA-3'UTR pair. It gives info. about how many sites are in the target sequence and p-Value indicating how likely to see this number of hits in this region.

2. `7mers_hits.txt`: gives info. about individual hit (start position of the hit in the target region).

## Example run on Linux command line:
```
python3 miRNA_targets_v0.3.py   3utr.fasta miRNA_sequences.fasta      7   
```
