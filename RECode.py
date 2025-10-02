# RECode.py
#
# Accepts a sequence in-frame, determines all possible synonymous encodings and
# figures out if any restriction enzymes can cut either the original or the re-encodings.

import getopt, sys, re
from Bio.Seq import Seq
from Bio.Restriction.PrintFormat import PrintFormat
from Bio.Restriction.Restriction import RestrictionBatch

reverse_codon_table = {"A": ["GCT", "GCC", "GCA", "GCG"],
                       "C": ["TGT", "TGC"],
                       "D": ["GAT", "GAC"],
                       "E": ["GAA", "GAG"],
                       "F": ["TTT", "TTC"],
                       "G": ["GGT", "GGC", "GGA", "GGG"],
                       "H": ["CAT", "CAC"],
                       "I": ["ATT", "ATC", "ATA"],
                       "K": ["AAA", "AAG"],
                       "L": ["TTA", "TTG", "CTT", "CTT", "CTC", "CTA", "CTG"],
                       "M": ["ATG"],
                       "N": ["AAT", "AAC"],
                       "P": ["CCT", "CCC", "CCA", "CCG"],
                       "Q": ["CAA", "CAG"],
                       "R": ["AGA", "AGG", "CGT", "CGC", "CGA", "CGG"],
                       "S": ["AGT", "AGC", "TCT", "TCC", "TCA", "TCG"],
                       "T": ["ACT", "ACC", "ACA", "ACG"],
                       "V": ["GTT", "GTC", "GTA", "GTG"],
                       "W": ["TGG"],
                       "Y": ["TAT", "TAC"],
                       "*": ["TAA", "TAG", "TGA"]}

codon_table = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
               "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
               "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
               "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
               "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
               "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
               "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
               "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
               "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
               "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
               "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
               "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
               "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
               "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
               "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
               "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}
               

help_text = """{0}
Accepts a sequence in reading frame, determines all the possible synonymous encodings, and
figures out if any restriction sites can be found.

{0} -s DNA_sequence -p AA_sequence 

version {1}
"""

version_number = "0.0.1"


nucleotide_regex_dict = {"A": "A",
                         "C": "C",
                         "G": "G",
                         "T": "T",
                         "N": ".",
                         "B": "[CGT]",
                         "D": "[AGT]",
                         "H": "[ACT]",
                         "V": "[ACG]",
                         "W": "[AT]",
                         "S": "[CG]",
                         "R": "[AG]",
                         "Y": "[CT]",
                         "M": "[AC]",
                         "K": "[GT]"}

def RE_list_to_regex_dict(RE_file):
    RE_regex_dict = {}
    max_RE_name_length = 0
    
    with open(RE_file, "r") as RE_fh:
        line = RE_fh.readline()
        
        while line:
            RE_name = line.split("\t")[0]
            RE_pattern = line.rstrip().split("\t")[1]
            
            regex_pattern = ""
            for nucleotide in RE_pattern:
                regex_pattern = regex_pattern + nucleotide_regex_dict[nucleotide]
            RE_regex_dict[RE_name] = regex_pattern
            
            if len(RE_name) > max_RE_name_length:
                max_RE_name_length = len(RE_name)
            else:
                pass
            
            line = RE_fh.readline()
    return RE_regex_dict, max_RE_name_length

def RE_list_to_IUPAC_dict(RE_file):
    RE_IUPAC_dict = {}
    max_RE_name_length = 0
    
    with open(RE_file, "r") as RE_fh:
        line = RE_fh.readline()
        
        while line:
            RE_name = line.split("\t")[0]
            RE_pattern = line.rstrip().split("\t")[1]
            
            RE_IUPAC_dict[RE_name] = RE_pattern
            
            if len(RE_name) > max_RE_name_length:
                max_RE_name_length = len(RE_name)
            else:
                pass
            
            line = RE_fh.readline()
    return RE_IUPAC_dict, max_RE_name_length

base_pairing_dict = {"A": "T",
                     "a": "t",
                     "C": "G",
                     "c": "g",
                     "G": "C",
                     "g": "c",
                     "T": "A",
                     "t": "a",
                     "N": "N",
                     "n": "n"}
                     
def complement(DNA_sequence):
    new_sequence = ""
    for nucleotide in DNA_sequence:
        new_sequence = new_sequence + base_pairing_dict[nucleotide]
    return new_sequence

def reverse_complement(DNA_sequence):
    return complement(DNA_sequence)[::-1]
    
def is_palindromic(pattern):
    if pattern.upper() == reverse_complement(pattern.upper()):
        return True
    else:
        return False

def enzyme_recognition_range(DNA_sequence, pattern):
    recognition_range_array = [" "]*len(DNA_sequence)
    
    matches = re.finditer(pattern, DNA_sequence)
    
    for match in matches:
        for i in range(match.start(), match.end()):
            recognition_range_array[i] = "-"

    revcomp_DNA_sequence = reverse_complement(DNA_sequence)
    matches = re.finditer(pattern, revcomp_DNA_sequence)
        
    for match in matches:
        for i in range(len(DNA_sequence)-match.end(), len(DNA_sequence)-match.start()):
            recognition_range_array[i] = "-"
    
    return "".join(recognition_range_array)
    
def reverse_translate(DNA_sequence, mutated_aa_sequence):
    aa_sequence = ""
    
    if mutated_aa_sequence == "":
        aa_sequence = translate(DNA_sequence.upper())
    else:
        aa_sequence = mutated_aa_sequence
    list_of_encodings = [""]
    for aa in aa_sequence:
        new_list_of_encodings = []
        for codon in reverse_codon_table[aa]:
            for encoding in list_of_encodings:
                new_list_of_encodings.append(encoding+codon)
        list_of_encodings = new_list_of_encodings
    
    # now lightly annotate the sequences:
    # lower case where mutations from original DNA_sequence are found
    new_list_of_encodings = []
    for encoding in list_of_encodings:
        annotated_encoding = ""
        if encoding == DNA_sequence:
            pass
        else:
            for index, nucleic_acid in enumerate(encoding):
                if nucleic_acid != DNA_sequence[index]:
                    annotated_encoding = annotated_encoding + nucleic_acid.lower()
                else:
                    annotated_encoding = annotated_encoding + nucleic_acid
            new_list_of_encodings.append(annotated_encoding)
    return new_list_of_encodings

def translate(DNA_sequence):
    sub_sequence = DNA_sequence
    aa_sequence = ""
    
    while len(sub_sequence) >= 3:
        codon = sub_sequence[0:3]
        aa_sequence = aa_sequence + codon_table[codon]
        
        sub_sequence = sub_sequence[3:len(sub_sequence)]
    return aa_sequence

def count_mismatches(DNA_sequence_1, DNA_sequence_2):
    num_of_mismatches = 0
    for index, nucleotide in enumerate(DNA_sequence_1):
        if nucleotide != DNA_sequence_2[index]:
            num_of_mismatches += 1
    return num_of_mismatches

def translation_string(DNA_sequence):
    aa_sequence = translate(DNA_sequence.upper())
    
    return_string = ""
    for aa in aa_sequence:
        return_string = return_string + aa + "  "
    return return_string

nucleotide_count_dict = {"A": 1,
                         "C": 1,
                         "G": 1,
                         "T": 1,
                         "N": 0,
                         "B": 1,
                         "D": 1,
                         "H": 1,
                         "V": 1,
                         "W": 1,
                         "S": 1,
                         "R": 1,
                         "Y": 1,
                         "M": 1,
                         "K": 1}

def count_number_of_recognition_nucleotides(IUPAC_pattern):
    count = 0
    for nucleotide in IUPAC_pattern:
        count = count + nucleotide_count_dict[nucleotide]
    return count

def main(argv):
    # Start of the script
    
    arg_n = ""
    arg_i = 0
    arg_f = 0
    arg_p = ""
    arg_c = 0
    arg_o = ""
    arg_help = help_text.format(argv[0], version_number)
    opts, args = getopt.getopt(argv[1:], "hn:i:f:p:c:o:v", ["help", "nucleotide=", "initial=", "final=", "peptide=", "cutter=", "output=", "version"])
    
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(0)
        elif opt in ("-v", "--version"):
            print(version_number)
            sys.exit(0)
        elif opt in ("-n", "--nucleotide"):
            arg_n = arg
            arg_f = len(arg_n)
        elif opt in ("-i", "--initial"):
            arg_i = int(arg)
        elif opt in ("-f", "--final"):
            arg_f = int(arg)
        elif opt in ("-c", "--cutter"):
            arg_c = int(arg)
        elif opt in ("-o", "--output"):
            arg_o = arg
        elif opt in ("-p", "--peptide"):
            arg_p = arg
    
    output_filename = ""
    
    if arg_o != "":
        output_filename = arg_o
    else:
        output_filename = arg_n + ".reencoding.RE.txt"
    
    
    RE_file = "NEB_restriction_enzyme_list_and_recognition_sequences.txt"
    RE_regex_dict, max_RE_name_length = RE_list_to_regex_dict(RE_file)
    RE_IUPAC_dict, max_RE_name_length = RE_list_to_IUPAC_dict(RE_file)
    
    #print(RE_regex_dict)
    
    with open(output_filename, "w") as output_fh:
        output_fh.write("ORIGINAL SEQUENCE\n\n")
        output_fh.write(translation_string(arg_n)+"\n")
        output_fh.write(arg_n+"\n")
        output_fh.write(complement(arg_n)+"\n")
        
        for RE in RE_regex_dict:
            IUPAC_pattern = RE_IUPAC_dict[RE]
            if arg_c == count_number_of_recognition_nucleotides(IUPAC_pattern) or arg_c == 0:
                regex_pattern = RE_regex_dict[RE]
                recognition_string = enzyme_recognition_range(arg_n.upper(), regex_pattern)
                if recognition_string.find("-") >= 0:
                    output_fh.write(recognition_string + "\t" + RE + " (" + IUPAC_pattern + ")\n")
                else:
                    pass
            else:
                pass
        
        for index, new_encoding in enumerate(reverse_translate(arg_n, arg_p)):
            output_fh.write("*********************************************************************************\n")
            output_fh.write(f"RE-ENCODING #{index + 1}\n\n")
            output_fh.write(translation_string(new_encoding)+"\n")
            output_fh.write(new_encoding+"\n")
            output_fh.write(complement(new_encoding)+"\n")

            for RE in RE_regex_dict:
                IUPAC_pattern = RE_IUPAC_dict[RE]
                if arg_c == count_number_of_recognition_nucleotides(IUPAC_pattern) or arg_c == 0:
                    regex_pattern = RE_regex_dict[RE]
                    recognition_string = enzyme_recognition_range(new_encoding.upper(), regex_pattern)
                    if recognition_string.find("-") >= 0:
                        output_fh.write(recognition_string + "\t" + RE + " (" + IUPAC_pattern + ")\n")
                    else:
                        pass
                else:
                    pass
             

if __name__ == "__main__":
    main(sys.argv)