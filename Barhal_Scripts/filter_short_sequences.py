input_file = "/workspace/jmarkwelchlab/P_0003_Neisseriaceae/CONTIG_LENGTH_TEST/03_GENOMES_EDITED/Edited_Concatenated_both.fa"
output_file = "/workspace/jmarkwelchlab/P_0003_Neisseriaceae/CONTIG_LENGTH_TEST/03_GENOMES_EDITED/Edited_Concatenated_both_seq_less_than_300.fa"
sequence = ""
keep_sequence = False

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        if line.startswith(">"):  # Defline
            if len(sequence) < 300 and keep_sequence:
                outfile.write(defline)
                outfile.write(sequence)
            defline = line
            sequence = ""
            keep_sequence = True
        else:  # Sequence
            sequence += line

    # Check the last sequence in the file
    if len(sequence) < 300 and keep_sequence:
        outfile.write(defline)
        outfile.write(sequence)

