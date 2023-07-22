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

    upstream_sequence = sequence[upstream_start+1:upstream_end]
    downstream_sequence = sequence[downstream_start+1:downstream_end]

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
            downstream_location = features[-1].locationWWW
            
            # Save the information in the dictionary
            sequences_with_info[record.id] = (sequence, features, length, upstream_location, downstream_location)
            upstream_sequence, downstream_sequence = get_upstream_downstream_sequences(sequence, upstream_location, downstream_location)
    return sequence, upstream_sequence, downstream_sequence


def homologous_recombination(target_seq, insert_seq, upstream_homology, downstream_homology):
    # Convert sequences to Python strings
    target_seq_str = str(target_seq)
    insert_seq_str = str(insert_seq)

    # Find the positions of upstream and downstream homologies in the target sequence
    up_index = target_seq_str.find(upstream_homology)
    down_index = target_seq_str.find(downstream_homology) + len(downstream_homology)
    print(up_index, down_index)
    
    # Perform the replacement (homologous recombination)
    new_sequence = target_seq_str[:up_index] + insert_seq_str + target_seq_str[down_index:]
    return new_sequence

if __name__ == "__main__":
    # Paths to input files
    target_file = "/mnt/d/Github/XFORM_test/input/SD36051-001.fasta"
    insert_fasta_file = "/mnt/d/Github/XFORM_test/input/ADD-4353.fasta"
    insert_genbank_file = "/mnt/d/Github/XFORM_test/input/ADD-4353.gb"
    
    # Read the target sequences (could be multiple)
    target_seqs = read_fasta(target_file)
    
    # Read the insert sequence and homology information from GenBank file
    insert_seq, upstream_homology, downstream_homology = read_genbank(insert_genbank_file)
    
    # Perform homologous recombination for each target sequence
    recombined_sequences = []
    for target_seq in target_seqs:
        new_sequence = homologous_recombination(target_seq, insert_seq, upstream_homology, downstream_homology)
        recombined_sequences.append(new_sequence)
    
    # Write the new sequences to a new FASTA file
    with open("recombined_sequences.fasta", "w") as output_file:
        for i, sequence in enumerate(recombined_sequences, 1):
            new_record = SeqRecord(Seq(sequence), id=f"Recombined_Sequence_{i}", description="Homologous Recombination Result")
            SeqIO.write(new_record, output_file, "fasta")
