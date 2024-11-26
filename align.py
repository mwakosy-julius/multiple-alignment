import streamlit as st
import numpy as np

MATCH = 1
MISMATCH = -1
GAP = -2

def parse_fasta_sequences(fasta_input):
    sequences = []
    current_sequence = ""
    
    lines = fasta_input.strip().splitlines()
    
    for line in lines:
        if line.startswith(">"):
            if current_sequence:
                sequences.append(current_sequence)
                current_sequence = ""
        else:
            current_sequence += line.strip()
    
    if current_sequence:
        sequences.append(current_sequence)
    
    return sequences

def needleman_wunsch(seq1, seq2):
    n, m = len(seq1), len(seq2)
    score = np.zeros((n + 1, m + 1))
    for i in range(n + 1):
        score[i][0] = i * GAP
    for j in range(m + 1):
        score[0][j] = j * GAP

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = score[i - 1][j - 1] + (MATCH if seq1[i - 1] == seq2[j - 1] else MISMATCH)
            delete = score[i - 1][j] + GAP
            insert = score[i][j - 1] + GAP
            score[i][j] = max(match, delete, insert)

    align1, align2 = "", ""
    i, j = n, m
    while i > 0 and j > 0:
        if score[i][j] == score[i - 1][j - 1] + (MATCH if seq1[i - 1] == seq2[j - 1] else MISMATCH):
            align1 += seq1[i - 1]
            align2 += seq2[j - 1]
            i -= 1
            j -= 1
        elif score[i][j] == score[i - 1][j] + GAP:
            align1 += seq1[i - 1]
            align2 += "-"
            i -= 1
        else:
            align1 += "-"
            align2 += seq2[j - 1]
            j -= 1
    while i > 0:
        align1 += seq1[i - 1]
        align2 += "-"
        i -= 1
    while j > 0:
        align1 += "-"
        align2 += seq2[j - 1]
        j -= 1
    return align1[::-1], align2[::-1]

def calculate_distance_matrix(sequences):
    n = len(sequences)
    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            align1, align2 = needleman_wunsch(sequences[i], sequences[j])
            if len(align1) != len(align2):
                raise ValueError("Aligned sequences must be of the same length.")
            score = sum(1 if a == b else -1 for a, b in zip(align1, align2))
            matrix[i][j] = -score
            matrix[j][i] = matrix[i][j]  # Ensure symmetry
    return matrix

def guide_tree(distance_matrix):
    n = len(distance_matrix)
    clusters = {i: [i] for i in range(n)}
    distances = {(i, j): distance_matrix[i, j] for i in range(n) for j in range(i + 1, n)}
    new_key = n  # Start generating unique keys after initial indices

    while len(clusters) > 1:
        # Find closest pair
        pair = min(distances, key=distances.get)
        i, j = pair
        new_cluster = clusters[i] + clusters[j]
        clusters[new_key] = new_cluster
        del clusters[i]
        del clusters[j]

        # Update distances
        for k in list(clusters.keys()):
            if k != new_key:
                dist = sum(distances.get((min(x, y), max(x, y)), 0) for x in new_cluster for y in clusters[k])
                dist /= len(new_cluster) * len(clusters[k])
                distances[(min(new_key, k), max(new_key, k))] = dist

        # Remove old distances
        distances = {k: v for k, v in distances.items() if i not in k and j not in k}
        new_key += 1

    return list(clusters.values())

def progressive_alignment(sequences, guide):
    alignments = {i: seq for i, seq in enumerate(sequences)}
    for group in guide:
        if len(group) == 2:
            align1, align2 = needleman_wunsch(alignments[group[0]], alignments[group[1]])
            alignments[group[0]] = align1
            alignments[group[1]] = align2
        else:
            # Handle more than two sequences in a group
            base_align = alignments[group[0]]
            for idx in group[1:]:
                base_align, new_align = needleman_wunsch(base_align, alignments[idx])
                alignments[group[0]] = base_align
                alignments[idx] = new_align
    return list(alignments.values())

st.title("Multiple Sequence Alignment App")

fasta_input = st.text_area(
        "Enter FASTA-formatted sequences:",
        height=200,
        placeholder=">seq1\nACTGCTAGCTAG\n>seq2\nACT-CTAGC-AG\n>seq3\nAC-GCT-GCT-G"
    )

if fasta_input.strip():
    sequences = parse_fasta_sequences(fasta_input)
    if len(sequences) < 2:
        st.error("Please enter at least two sequences in FASTA format.")

    with st.spinner("Performing multiple sequence alignment..."):
        distance_matrix = calculate_distance_matrix(sequences)
        guide = guide_tree(distance_matrix)
        aligned_sequences = progressive_alignment(sequences, guide)

        st.success("Alignment complete!")
        st.write(aligned_sequences)

