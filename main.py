import os
import numpy as np
from collections import Counter


AMINO_ACIDS = sorted("ACDEFGHIKLMNPQRSTVWY")

def fasta_processor(fasta_dir):
    fasta_file_found = False
    summary_data = []

    aminoacid_result_path = os.path.join(fasta_dir, "analysis_aminoacid.txt")
    folder_summary_path = os.path.join(fasta_dir, "folder_summary.txt")

    with open(aminoacid_result_path, 'w') as file:
        file.write("seq_id," + ",".join(AMINO_ACIDS) + "\n")

    for file_name in os.listdir(fasta_dir):
        if file_name.endswith('.fasta') or file_name.endswith('.fa'):
            fasta_file_found = True
            file_path = os.path.join(fasta_dir, file_name)
            print(f"[INFO] Processing file: {file_name}")
            seq_count = process_fasta_file(file_path, aminoacid_result_path)
            summary_data.append(f"{file_name}, {seq_count}")

    if not fasta_file_found:
        print("[ERROR] No FASTA files found in the folder. Please check file extensions.")
    else:
        with open(folder_summary_path, 'w') as summary_file:
            summary_file.write("Filename, Total Sequences\n")
            for line in summary_data:
                summary_file.write(line + "\n")
        print(f"[DONE] Summary saved to: {folder_summary_path}")

def process_fasta_file(file_path, aminoacid_result_path):
    sequences = {}
    sequence_id = ""
    sequence = ""
    total_sequences = 0

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence_id:
                    sequences[sequence_id] = sequence
                sequence_id = line.split('|')[1][1:]
                sequence = ""
                total_sequences += 1
            else:
                sequence += line
        if sequence_id:
            sequences[sequence_id] = sequence  # save the last sequence

    aminoacid_processor(sequences, aminoacid_result_path)
    return total_sequences

def aminoacid_processor(sequences, result_path):
    vectors = {} 

    with open(result_path, 'a') as file:
        for seq_id, seq in sequences.items():
            length = len(seq)
            char_freq = Counter(seq)

            
            freq_vector = np.array(
                [char_freq.get(aa, 0) / length for aa in AMINO_ACIDS],
                dtype=np.float32
            )

           
            file.write(seq_id)
            for freq in freq_vector:
                file.write(f",{freq:.4f}")
            file.write("\n")

            vectors[seq_id] = freq_vector

    return vectors



if __name__ == "__main__":
    fasta_dir = input("Enter directory of your FASTA files:\n").strip()
    if os.path.isdir(fasta_dir):
        fasta_processor(fasta_dir)
    else:
        print("[ERROR] Provided path is not a valid directory.")
