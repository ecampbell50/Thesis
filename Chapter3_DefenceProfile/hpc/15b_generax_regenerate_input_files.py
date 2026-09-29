# ======================================================================
# Thesis Chapter 3, Figure 3.7
# Rebuilds the GeneRax family/mapping input files.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/11_DefenseProfileofSsuis/1_SsuisGenomes/ssuis_defence_project/GeneReconcilliation_of_tmn/Phylogeny_alignments/regenerate_files.py
# ======================================================================
import re
from ete3 import Tree

t = Tree("rooted_tree.nwk")
species_tree_labels = t.get_leaf_names()

accession_to_species = {}
for label in species_tree_labels:
    match = re.search(r'(GC[AF]_\d+_\d+)$', label)
    if match:
        accession_to_species[match.group(1)] = label

mapping_lines = []
new_alignment = []

with open("Tmn_MAFFT_aligned.faa") as f:
    for line in f:
        if line.startswith(">"):
            original_header = line.strip().lstrip(">")
            after_pipe = original_header.split("|")[1]
            accession_with_dot = after_pipe.split("_")[0] + "_" + after_pipe.split("_")[1]
            accession_normalised = accession_with_dot.replace(".", "_")

            if accession_normalised in accession_to_species:
                species_label = accession_to_species[accession_normalised]
                new_alignment.append(">" + species_label + "\n")
                mapping_lines.append(f"{species_label}\t{species_label}")
                print(f"OK: {original_header} -> {species_label}")
            else:
                print(f"WARNING: Could not map {original_header}")
                new_alignment.append(line)
        else:
            new_alignment.append(line)

with open("Tmn_MAFFT_aligned_cleaned.faa", "w") as f:
    f.writelines(new_alignment)

with open("Tmn_mapping.link", "w") as f:
    f.writelines(line + "\n" for line in mapping_lines)

print(f"\nDone. Mapped {len(mapping_lines)} sequences.")
