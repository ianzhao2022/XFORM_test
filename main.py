from Bio import SeqIO, Align
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def read_fasta(filename):
    # Read the FASTA file and return a list of sequences
    sequences = []
    with open(filename, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            sequences.append(record.seq)
    return sequences

def get_upstream_downstream_sequences(sequence, upstream_location, downstream_location):
    # Extracts and returns the upstream and downstream sequences from the main sequence.
   
    upstream_start = upstream_location.start
    upstream_end = upstream_location.end
    downstream_start = downstream_location.start
    downstream_end = downstream_location.end

    upstream_sequence = sequence[upstream_start:upstream_end]
    downstream_sequence = sequence[downstream_start:downstream_end]

    return upstream_sequence, downstream_sequence

def read_genbank(file_path):
    # Reads a GenBank file and returns a dictionary containing the sequence names as keys and a tuple of
    # (sequence, features, length, upstream_location, downstream_location) as values.

    sequences_with_info = {}

    with open(file_path, 'r') as file:
        for record in SeqIO.parse(file, 'genbank'):
            sequence = str(record.seq)
            features = record.features
            length = len(sequence)
            
            # Find the first and last feature locations and save them as upstream and downstream locations
            upstream_location = features[0].location
            downstream_location = features[-1].location
            
            # Save the information in the dictionary
            sequences_with_info[record.id] = (sequence, features, length, upstream_location, downstream_location)

    return sequences_with_info

def get_complement(dna_sequence):
    # Get the complementary sequence of a DNA sequence.

    # Define the translation table to replace each nucleotide with its complement
    complement_table = str.maketrans("ATCG", "TAGC")
    
    # Use the translate method to get the complementary sequence
    complementary_sequence = dna_sequence.translate(complement_table)
    
    return complementary_sequence

def search_upstream_downstream_in_sequence(main_sequence, upstream_sequence, downstream_sequence, cutsite_sequence):
    # Searches for the location of upstream and downstream sequences in the main sequence.

    def find_sequence_in_main(sequence):
        seq_variants = [sequence, get_complement(sequence), sequence[::-1], get_complement(sequence[::-1])]
        for variant_type, seq_type in enumerate(seq_variants, start=1):
            index = main_sequence.find(seq_type)
            if index != -1:
                print(f"sequence (or its variant) found at index {index}:{index+len(sequence)}")
                return index, variant_type
        print("sequence and its variants not found")
        return -1, 0
        
    print("Upstream ", end='')
    upstream_result = find_sequence_in_main(upstream_sequence)
    print("Downstream ", end='')
    downstream_result = find_sequence_in_main(downstream_sequence)
    print("Cut site ", end='')
    find_sequence_in_main(cutsite_sequence)\
        
    if upstream_result[1] != downstream_result[1]:
        return KeyError("Upstream and downstream sequences are not the same type. Please check the sequences.")
    
    return (upstream_result[0], upstream_result[0] + len(upstream_sequence)), (downstream_result[0], downstream_result[0] + len(downstream_sequence)), upstream_result[1]

def homologous_recombination(target_seq, insert_seq, location_in_sequence):
    # Convert sequences to Python strings
    target_seq_str = str(target_seq)
    insert_seq_str = str(insert_seq)

    # get the positions of upstream and downstream homologies in the location information handle the special case
    upstream_location_start = location_in_sequence[0][0]
    upstream_location_end = location_in_sequence[0][1]
    downstream_location_start = location_in_sequence[1][0]
    downstream_location_end = location_in_sequence[1][1]
    match_type = location_in_sequence[2]
    
    # Perform the replacement (homologous recombination)
    if match_type ==1:
        new_sequence = target_seq_str[:upstream_location_start] + insert_seq_str + target_seq_str[downstream_location_end:]
    elif match_type ==2:
        new_sequence = target_seq_str[:upstream_location_start] + get_complement(insert_seq_str) + target_seq_str[downstream_location_end:]
    elif match_type ==3:
        new_sequence = target_seq_str[:downstream_location_start] + insert_seq_str[::-1] + target_seq_str[upstream_location_end:]
    elif match_type ==4:
        new_sequence = target_seq_str[:downstream_location_start] + get_complement(insert_seq_str[::-1]) + target_seq_str[upstream_location_end:]
    else:
        new_sequence = target_seq_str
    
    return new_sequence

if __name__ == "__main__":
    # Paths to input files
    target_file = "/mnt/d/Github/XFORM_test/input/SD36051-001.fasta"
    insert_genbank_file = "/mnt/d/Github/XFORM_test/input/ADD-4353.gb"
    # Set the cutsite sequence
    cutsite_sequence = "tctggcgcaggtgatatgta"
    
    # Read the target sequences (could be multiple)
    target_seqs = read_fasta(target_file)
    
    # Read the insert sequence and homology information from GenBank file
    genbank_sequences = read_genbank(insert_genbank_file)
    
    # Process and print information about insert sequence(upstream and downstream sequences, features, etc.)
    for sequence_name, data in genbank_sequences.items():
        sequence, features, length, upstream_location, downstream_location = data
        print(f"Sequence name: {sequence_name}")
        print(f"Sequence length: {length} base pairs")
        print(f"Upstream location: {upstream_location}")
        print(f"Downstream location: {downstream_location}")

        # Extract and print upstream and downstream sequences
        upstream_sequence, downstream_sequence = get_upstream_downstream_sequences(sequence, upstream_location, downstream_location)
        print(f"Upstream sequence: {upstream_sequence}")
        print(f"Downstream sequence: {downstream_sequence}")

        # Print Features
        print("Features:")
        for feature in features:
            print(f"\tFeature type: {feature.type}")
            print(f"\tLocation: {feature.location}")
            print(f"\tQualifiers: {feature.qualifiers}")
            
    # Set the insert sequence
    insert_seq = sequence
        
    # Perform homologous recombination for each target sequence
    recombined_sequences = []
    for target_seq in target_seqs:
        # Get location of upstream and downstream sequences in the target sequence
        location_in_sequence= search_upstream_downstream_in_sequence(target_seq, upstream_sequence, downstream_sequence, cutsite_sequence.upper())
        # Perform homologous recombination
        new_sequence = homologous_recombination(target_seq, insert_seq, location_in_sequence)
        recombined_sequences.append(new_sequence)
    
    # Write the new sequences to a new FASTA file
    with open("recombined_sequences.fasta", "w") as output_file:
        for i, sequence in enumerate(recombined_sequences, 1):
            new_record = SeqRecord(Seq(sequence), id=f"Recombined_Sequence_{i}", description="")
            SeqIO.write(new_record, output_file, "fasta")
