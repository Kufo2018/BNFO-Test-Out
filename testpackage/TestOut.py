import math

listOfSequences = []

# open a (new) file to write
output_file = open("Output.txt", "w")


# Reads FASTA file
def read_file(file_path):
    # Locates and reads FASTA file
    with open(file_path) as file:
        aligned_sequences = file.readlines()

    # Filters raw list of aligned sequences
    for line in aligned_sequences:
        if line[0] != ">":
            listOfSequences.append(line.rstrip('\n'))


# Pairs sequences
def pair_sequences(sequence_list):
    for x in sequence_list:
        for y in sequence_list:
            calculate_jc(x, y)
        output_file.write("\n")


# Calculates JC
def calculate_jc(sequence_a, sequence_b):
    aligns_but_mismatched = 0
    ungapped_sites = 0
    values = []

    # Determines the p value
    for x, y in zip(sequence_a, sequence_b):

        # Aligned but mismatched sites
        if x != y and x != "-" and y != "-":
            aligns_but_mismatched += 1

        # Ungapped sites
        if x != "-" and y != "-":
            ungapped_sites += 1

    p = aligns_but_mismatched / ungapped_sites

    distance = round(-3 / 4 * math.log(1 - ((4 / 3) * p)), 4)

    values.append(distance)

    for z in values:
        output_file.write(str(format(z, '.4f')))
        output_file.write(" ")


read_file("8.1_AlignedDNA.fas")
pair_sequences(listOfSequences)
output_file.close()
