#!/usr/bin/env python3
"""
find_tmn_cluster_medoids.py
Find the most representative (medoid) sequence from each Tmn cluster

Usage:
    python3 find_tmn_cluster_medoids.py
"""

from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from pathlib import Path

def load_cluster_ids(cluster_file):
    """Load sequence IDs from cluster file"""
    ids = []
    with open(cluster_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:  # Skip empty lines
                ids.append(line)
    return ids

def extract_cluster_sequences(fasta_file, cluster_ids):
    """Extract sequences matching cluster IDs from multifasta"""
    cluster_seqs = []
    print(f"  Looking for {len(cluster_ids)} sequences...")
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in cluster_ids:
            cluster_seqs.append(record)
    
    print(f"  Found {len(cluster_seqs)} sequences")
    return cluster_seqs

def calculate_pairwise_identity(seq1, seq2):
    """Calculate percentage identity between two sequences"""
    # Use global alignment
    alignments = pairwise2.align.globalxx(str(seq1), str(seq2), 
                                          one_alignment_only=True)
    
    if not alignments:
        return 0.0
    
    alignment = alignments[0]
    matches = alignment[2]  # Number of matches
    max_len = max(len(seq1), len(seq2))
    
    identity = (matches / max_len) * 100
    return identity

def find_medoid(sequences):
    """
    Find the medoid sequence - the one with minimum average distance to all others
    
    The medoid is the most "central" sequence in the cluster
    """
    n = len(sequences)
    print(f"  Calculating pairwise distances for {n} sequences...")
    
    # Calculate distance matrix (distance = 100 - identity)
    distances = np.zeros((n, n))
    
    for i in range(n):
        if i % 10 == 0:
            print(f"    Processed {i}/{n} sequences...")
        
        for j in range(i+1, n):
            identity = calculate_pairwise_identity(
                sequences[i].seq, 
                sequences[j].seq
            )
            distance = 100 - identity  # Convert identity to distance
            distances[i][j] = distance
            distances[j][i] = distance
    
    # Find sequence with minimum average distance to all others
    avg_distances = distances.mean(axis=1)
    medoid_idx = np.argmin(avg_distances)
    
    medoid_seq = sequences[medoid_idx]
    avg_identity = 100 - avg_distances[medoid_idx]
    
    print(f"  Medoid: {medoid_seq.id}")
    print(f"  Average identity to cluster: {avg_identity:.1f}%")
    
    return medoid_seq, avg_identity

def compare_medoids(medoid1, medoid2):
    """Compare two medoid sequences"""
    print("\n" + "="*60)
    print("COMPARING MEDOIDS")
    print("="*60)
    
    identity = calculate_pairwise_identity(medoid1.seq, medoid2.seq)
    
    print(f"\nLarge cluster medoid: {medoid1.id}")
    print(f"  Length: {len(medoid1.seq)} aa")
    print(f"\nSmall cluster medoid: {medoid2.id}")
    print(f"  Length: {len(medoid2.seq)} aa")
    print(f"\nSequence identity: {identity:.1f}%")
    print(f"Sequence divergence: {100-identity:.1f}%")
    
    # Perform alignment for visualization
    print("\nPerforming pairwise alignment...")
    alignments = pairwise2.align.globalxx(str(medoid1.seq), str(medoid2.seq))
    
    if alignments:
        print("\nAlignment preview (first 200 positions):")
        alignment = alignments[0]
        print(alignment[0][:200])  # Sequence 1
        print(alignment[1][:200])  # Sequence 2
        print(f"\nFull alignment score: {alignment[2]}")

def main():
    # File paths
    large_cluster_ids_file = "LargeCluster_Tmngenes.txt"
    small_cluster_ids_file = "SmallCluster_Tmngenes.txt"
    all_sequences_file = "all_tmn_protein.fasta"  # Your multifasta file
    
    print("="*60)
    print("FINDING TMN CLUSTER MEDOIDS")
    print("="*60)
    
    # Load cluster IDs
    print("\n1. Loading cluster IDs...")
    large_ids = load_cluster_ids(large_cluster_ids_file)
    small_ids = load_cluster_ids(small_cluster_ids_file)
    
    print(f"  Large cluster: {len(large_ids)} sequences")
    print(f"  Small cluster: {len(small_ids)} sequences")
    
    # Extract sequences for each cluster
    print("\n2. Extracting large cluster sequences...")
    large_sequences = extract_cluster_sequences(all_sequences_file, large_ids)
    
    print("\n3. Extracting small cluster sequences...")
    small_sequences = extract_cluster_sequences(all_sequences_file, small_ids)
    
    # Find medoids
    print("\n4. Finding large cluster medoid...")
    large_medoid, large_avg_id = find_medoid(large_sequences)
    
    print("\n5. Finding small cluster medoid...")
    small_medoid, small_avg_id = find_medoid(small_sequences)
    
    # Compare medoids
    compare_medoids(large_medoid, small_medoid)
    
    # Save medoid sequences
    print("\n6. Saving medoid sequences...")
    
    # Save individual medoids
    SeqIO.write(large_medoid, "large_cluster_medoid.fasta", "fasta")
    SeqIO.write(small_medoid, "small_cluster_medoid.fasta", "fasta")
    print("  Saved: large_cluster_medoid.fasta")
    print("  Saved: small_cluster_medoid.fasta")
    
    # Save both together for easy alignment
    both_medoids = [large_medoid, small_medoid]
    SeqIO.write(both_medoids, "both_cluster_medoids.fasta", "fasta")
    print("  Saved: both_cluster_medoids.fasta")
    
    # Save summary
    with open("medoid_summary.txt", 'w') as f:
        f.write("TMN CLUSTER MEDOID ANALYSIS\n")
        f.write("="*60 + "\n\n")
        f.write(f"Large cluster medoid: {large_medoid.id}\n")
        f.write(f"  Length: {len(large_medoid.seq)} aa\n")
        f.write(f"  Average identity to large cluster: {large_avg_id:.1f}%\n\n")
        f.write(f"Small cluster medoid: {small_medoid.id}\n")
        f.write(f"  Length: {len(small_medoid.seq)} aa\n")
        f.write(f"  Average identity to small cluster: {small_avg_id:.1f}%\n\n")
        f.write(f"Identity between medoids: {calculate_pairwise_identity(large_medoid.seq, small_medoid.seq):.1f}%\n")
    
    print("  Saved: medoid_summary.txt")
    
    print("\n" + "="*60)
    print("DONE!")
    print("="*60)
    print("\nNext steps:")
    print("1. Align medoids: mafft both_cluster_medoids.fasta > medoids_aligned.fasta")
    print("2. Visualize: Use Jalview, Snapgene, or any alignment viewer")
    print("3. Analyze: Look for key differences in domains, active sites, etc.")

if __name__ == "__main__":
    main()
