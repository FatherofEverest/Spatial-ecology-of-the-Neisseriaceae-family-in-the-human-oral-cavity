# Define the mapping from old names to new names
mapping = {}

# Read the tab-separated text file and create a mapping
with open("/workspace/jmarkwelchlab/P_0003_Neisseriaceae/CONTIG_LENGTH_TEST/03_GENOMES_EDITED/both_decompose.tsv", "r") as mapping_file:
    for line in mapping_file:
        old_name, new_name = line.strip().split("\t")
        mapping[old_name] = new_name

# Process the FASTA file and update deflines
with open("/workspace/jmarkwelchlab/P_0003_Neisseriaceae/CONTIG_LENGTH_TEST/03_GENOMES_EDITED/Edited_Concatenated_both_seq_less_than_300.fa", "r") as fasta_file, open("/workspace/jmarkwelchlab/P_0003_Neisseriaceae/CONTIG_LENGTH_TEST/03_GENOMES_EDITED/Edited_Concatenated_both_seq_less_than_300_renamed.fa", "w") as output_file:
    for line in fasta_file:
        if line.startswith(">"):
            # Extract the old name from the defline
            old_name = line.strip()[1:]
            # Get the corresponding new name from the mapping
            new_name = mapping.get(old_name, new_name)
            # Extract the last four digits from the old name
            last_four_digits = old_name[-4:]
            # Append the last four digits to the new name
            new_name += "_contig_" + last_four_digits
            # Write the updated defline
            output_file.write(">" + new_name + "\n")
        else:
            # Write the sequence lines unchanged
            output_file.write(line)

