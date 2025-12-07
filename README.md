# RECode.py
__R__estriction-__E__nzyme-based genotyping enabled by re-en__CODE__d sequences

Ever wanted to use Cas9 to generate a point mutation, but could not figure out a way to generate a new restriction enzyme site to facilitate genotyping later? Perhaps you also wanted to mutate the original guide RNA recognition sequence (~20 nt) to ensure edited DNA did not get cut again?

Then RECode.py is your friend!

Download both the Python script and the ```NEB_restriction_enzyme_list_and_recognition_sequences.txt```
 file to your working folder. Then run
```python
python3 RECode.py -n <ORIGINAL_NUCLEOTIDE_SEQ> -p <MUTATED_PEPTIDE_SEQ> -o <OUTPUT> -c <CUTTER_RECOGNITION_LENGTH>
```

| required | description |
|---|---|
|-n|Original nucleotide sequence. For ease of implementation, the first nucleotide defines the +1 reading frame and no other frames will be considered.|

|optional | description |
|---|---|
|-p |The translation you want to mutate to. If the argument is not given, RECode will only re-encode the WT sequence to look for RE sites|
|-o |Output file name. By default it's a long file name with your WT sequence in it.|
|-c |Specifies 'c'-cutters (eg. 6-cutters). If the argument is not given, RECode will consider all REs in the restriction enzymes text file.|


For example, RECoding the DNA sequence 5'-CGTCAATCGTATGGA-3' (translation RQSYG) to give **C**QSYG instead, and see which 6-cutters can cut the WT and recoded sequences, and save the results to output.txt.
```python
python3 RECode.py -n CGTCAATCGTATGGA -p CQSYG -c 6 -o output.txt
```

Requires only standard Python3 (no BioPython required)!
