# RECode.py
Restriction-Enzyme-based genotyping enabled by re-enCODEd sequences

Ever wanted to use Cas9 to generate a point mutation, but could not figure out a way to generate a new restriction enzyme site to facilitate genotyping later? Perhaps you also wanted to mutate the original guide RNA recognition sequence (~20 nt) to ensure edited DNA did not get cut again?

Then RECode.py is your friend!

Download both the Python script and the restriction_enzymes.txt file to your working folder. Then run
```
RECode.py -n original_nucleotide_sequence -p mutated_peptide_sequence -o output_file -c cutter_recognition_length
```

For example, RECoding the DNA sequence 5'-CGTCAATCGTATGGA-3' (translation RQSYG) to give **C**QSYG instead, and see which 6-cutters can cut the recoded sequences.
```
RECode.py -n CGTCAATCGTATGGA -p CQSYG -c 6 -o recoded.6-cutters.txt
```

Requires only standard Python3 (no BioPython required)!
