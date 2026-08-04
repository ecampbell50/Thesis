#!/usr/bin/env python3
"""
create_itol_serotype_dataset.py
Generate iToL COLORSTRIP dataset from serotype dataframe with custom colours

Usage in Jupyter:
    from create_itol_serotype_dataset import create_serotype_itol_dataset
    
    create_serotype_itol_dataset(
        serotype_df=serotype_df,
        serotype_colours_file='serotype_colours.csv',
        output_file='serotype_dataset_itol.txt'
    )
"""

import pandas as pd

def create_serotype_itol_dataset(serotype_df, serotype_colours_file, output_file='serotype_dataset_itol.txt'):
    """
    Create iToL COLORSTRIP dataset from serotype dataframe
    
    Parameters:
    -----------
    serotype_df : pandas.DataFrame
        DataFrame with 'genome_id' and 'Serotype_Consensus' columns
    serotype_colours_file : str
        Path to CSV file with columns: serotype, colour, name
    output_file : str
        Output filename for iToL dataset
    """
    
    # Load serotype colours (no header in CSV)
    print(f"Loading serotype colours from: {serotype_colours_file}")
    serotype_colours = pd.read_csv(serotype_colours_file, header=None, names=['serotype', 'colour', 'name'])
    print(f"  Loaded {len(serotype_colours)} serotype colours")
    
    # Create colour mapping dictionary
    colour_map = {}
    for _, row in serotype_colours.iterrows():
        serotype = str(row['serotype']).strip().strip('"')  # Remove quotes if present
        colour = str(row['colour']).strip().strip('"')
        colour_map[serotype] = colour
    
    print(f"\nColour mapping:")
    for serotype, colour in colour_map.items():
        print(f"  {serotype}: {colour}")
    
    # Extract genome_id and serotype from dataframe
    print(f"\nProcessing {len(serotype_df)} genomes...")
    
    # Handle different possible column names for genome_id
    if 'genome_id' in serotype_df.columns:
        genome_col = 'genome_id'
    elif 'Strain' in serotype_df.columns:
        genome_col = 'Strain'
    else:
        genome_col = serotype_df.index.name or serotype_df.columns[0]
        print(f"Using column: {genome_col}")
    
    # Get genome IDs and serotypes
    genome_data = []
    serotype_counts = {}
    
    for idx, row in serotype_df.iterrows():
        if genome_col == 'Strain':
            genome_id = str(row[genome_col]).replace('.gff', '').strip()
        elif 'genome_id' in serotype_df.columns:
            genome_id = str(row['genome_id']).strip()
        else:
            genome_id = str(idx).strip()
        
        serotype = str(row['Serotype_Consensus']).strip()
        
        # Handle missing/unknown serotypes
        if serotype == '-' or serotype == 'nan' or serotype == '' or pd.isna(serotype):
            serotype = 'Unknown'
        
        # Handle 'group s' type entries (convert to Unknown)
        if 'group' in serotype.lower():
            serotype = 'Unknown'
        
        # Get colour for this serotype
        if serotype in colour_map:
            colour = colour_map[serotype]
        else:
            # Assign greyscale for unlisted serotypes
            colour = '#999999'
            print(f"  Warning: No colour defined for serotype '{serotype}', using grey")
        
        genome_data.append({
            'genome_id': genome_id,
            'serotype': serotype,
            'colour': colour
        })
        
        # Count serotypes
        serotype_counts[serotype] = serotype_counts.get(serotype, 0) + 1
    
    print(f"\nSerotype distribution:")
    for serotype in sorted(serotype_counts.keys()):
        count = serotype_counts[serotype]
        print(f"  {serotype}: {count} genomes")
    
    # Create iToL dataset file
    print(f"\nWriting iToL dataset to: {output_file}")
    
    with open(output_file, 'w') as f:
        # Header
        f.write("DATASET_COLORSTRIP\n")
        f.write("#Serotype dataset - selected serotypes in colour, others in greyscale\n")
        f.write("\n")
        f.write("SEPARATOR SPACE\n")
        f.write("\n")
        f.write("DATASET_LABEL Serotype\n")
        f.write("\n")
        f.write("COLOR #ff0000\n")
        f.write("\n")
        f.write("COLOR_BRANCHES 0\n")
        f.write("\n")
        
        # Legend - get unique serotypes in order
        unique_serotypes = sorted(set([x['serotype'] for x in genome_data]))
        unique_colours = [colour_map.get(s, '#999999') for s in unique_serotypes]
        
        f.write("LEGEND_TITLE Serotype\n")
        f.write("LEGEND_POSITION_X 100\n")
        f.write("LEGEND_POSITION_Y 100\n")
        f.write("LEGEND_HORIZONTAL 1\n")
        f.write("LEGEND_SHAPES " + " ".join(["2"] * len(unique_serotypes)) + "\n")
        f.write("LEGEND_COLORS " + " ".join(unique_colours) + "\n")
        f.write("LEGEND_LABELS " + " ".join(unique_serotypes) + "\n")
        f.write("LEGEND_SHAPE_SCALES " + " ".join(["1"] * len(unique_serotypes)) + "\n")
        f.write("\n")
        
        # Data section
        f.write("DATA\n")
        
        for item in genome_data:
            f.write(f"{item['genome_id']} {item['colour']} {item['serotype']}\n")
    
    print(f"\nDone! Created iToL dataset with {len(genome_data)} genomes")
    print(f"\nTo use in iToL:")
    print(f"1. Open your tree in iToL")
    print(f"2. Go to 'Datasets' tab")
    print(f"3. Drag and drop '{output_file}'")
    print(f"4. Your serotypes will be displayed!")
    
    return output_file


create_serotype_itol_dataset(
    serotype_df=serotype_df,
    serotype_colours_file='serotype_colours.csv',
    output_file='serotype_dataset_itol.txt'
)

# Standalone execution example
if __name__ == "__main__":
    print("This script is meant to be imported in Jupyter")
    print("\nExample usage:")
    print("  from create_itol_serotype_dataset import create_serotype_itol_dataset")
    print("  create_serotype_itol_dataset(serotype_df, 'serotype_colours.csv', 'output.txt')")
