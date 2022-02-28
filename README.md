# MutaGenie
When creating small mutations to a protein, there are many options on the DNA sequence level. However, it is ideal to create or disrupt a restriction site, such that a successful mutant can be identified on an agarose gel after an analytic digest. MutaGenie helps identify all DNA sequence variants that are easiest to mutate and identify.
MutaGenie uses an excel sheet as an easy user interface, where you can input your data, adjust parameters and choose which restriction enzymes to consider, based on your own stock of available enzymes.

## Before starting:

Install the Python module "xlrd".

## While using MutaGenie:

Paste the full plasmid sequence in cell G2. Paste the coding sequence of the gene in cell G5. Indicate the desired mutation in cell G8 (e.g. "S101A"). Make sure the spreadsheet and the script are in the same folder.

Optional: adjust parameters.

## Output:

### Table:
Rank: rank based on log10 separation.

Enzyme: the restriction enzyme used to distinguish mutated from non-mutated plasmids

Log10 separation: the largest log10 separation between a band that is present in mutated/non-mutated plasmid but not in the other, compared to the nearest band on the other plasmid.

Band changes: only shows bands that have changed.

Mutations: number of point mutations required.

Mutation range: the total length of mutated DNA sequence.

### More details:
Mutations: a list of mutations required using IUPAC codes and position in relation to the beginning of the plasmid sequence.

Old -> new bands: a complete list of all bands in WT and mutant digest.
