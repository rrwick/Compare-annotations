# Compare annotations

This repo contains a script (`compare_annotations.py`) for quantifying the improvement in an annotation when a genome is reassembled and/or reannotated.

For example, imagine you had an annotated bacterial genome that's a couple of years old. You've now come back to this genome with new versions of the assembler and annotator and made an updated version. This script can tell you how things changed at the gene level. Hopefully they got better!

This script works by doing a global alignment of the genes of one annotation to the other. This means the two genomes must be roughly aligned at the gene level – i.e. they should start and end at the same places. If your new genome contains structural rearrangements, that will break this script!

A few other things to note:
* The input annotated genomes must be in GenBank format.
* This script only looks at 'CDS' features in the genomes, nothing else.



## Requirements

This script uses Python 3 and [Biopython](https://biopython.org/). If you can run `python3 -c "import Bio"` without getting an error, you should be good to go!

No installation is required: just clone the repo and run the script:
```bash
git clone https://github.com/rrwick/Compare-annotations.git
Compare-annotations/compare_annotations.py --help
```



## Example usage

Here are two versions of a genome you can try this script on: [CP001172.1](https://www.ncbi.nlm.nih.gov/nuccore/CP001172.1) and [CP001172.2](https://www.ncbi.nlm.nih.gov/nuccore/CP001172.2).

Download them in GenBank format and then run the script like this:
```
compare_annotations.py CP001172.1.gb CP001172.2.gb > results
```



## Output

This script outputs a gene-by-gene analysis.

When two CDSs are identical, you'll see something like this:
```
Exact match
  old: Anhydro-N-acetylmuramic acid kinase(AnhMurNAc kinase) (14970-16098 -, 1128 bp)
  new: Anhydro-N-acetylmuramic acid kinase (14970-16098 -, 1128 bp)
```

Or if the CDSs are similar (i.e. the same gene) but not identical, you'll see something like this:
```
Inexact match (98.54% ID, old seq longer)
  old: tyrosyl-tRNA synthetase (16159-17392 +, 1233 bp)
  new: Tyrosine--tRNA ligase (16177-17392 +, 1215 bp)
```

And if a CDS is only in one of the two assemblies, you might see stuff like this:
```
In old but not in new:
  Glutathione S-transferase, C-terminal domain protein (24047-24413 +, 366 bp)

In new but not in old:
  Sodium:neurotransmitter symporter family protein (25765-27244 +, 1479 bp)
```



## Summarising output

Here are a few lines of Bash to generate a summary file:
```bash
printf "Features in old assembly: %5s\n" $(grep "Features in old assembly" results | grep -oP "\d+") >> summary
printf "Features in new assembly: %5s\n" $(grep "Features in new assembly" results | grep -oP "\d+") >> summary
printf "\n" >> summary
printf "Exact match:              %5s\n" $(grep -c "Exact match" results) >> summary
printf "\n" >> summary
printf "Inexact match:            %5s\n" $(grep -c "Inexact match" results) >> summary
printf "  same length:            %5s\n" $(grep -c "same length" results) >> summary
printf "  new seq longer:         %5s\n" $(grep -c "new seq longer" results) >> summary
printf "  old seq longer:         %5s\n" $(grep -c "old seq longer" results) >> summary
printf "\n" >> summary
printf "In new but not in old:    %5s\n" $(grep -c "In new but not in old" results) >> summary
printf "In old but not in new:    %5s\n" $(grep -c "In old but not in new" results) >> summary
printf "\n" >> summary
printf "No longer hypothetical:   %5s\n" $(grep -c "no longer hypothetical" results) >> summary
printf "Still hypothetical:       %5s\n" $(grep -c "still hypothetical" results) >> summary
printf "Became hypothetical:      %5s\n" $(grep -c "became hypothetical" results) >> summary
```

The result of which should look something like this:
```
Features in old assembly:  3458
Features in new assembly:  3427

Exact match:               2937

Inexact match:              354
  same length:                3
  new seq longer:           286
  old seq longer:            65

In new but not in old:      136
In old but not in new:      167

No longer hypothetical:     239
Still hypothetical:         638
Became hypothetical:         84
```



## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
