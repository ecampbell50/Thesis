#!/usr/bin/env python3
# iTOL presence/absence shapes for selected defence subtypes (Fig 3.6).


import pandas as pd
import numpy as np

def create_shapes_dataset_presence_absence(df, subtypes, output_file='shapes_dataset.txt',
                                            dataset_label='Defense Systems', colors=None,
                                            shape_type=2, unique_colors=True):
    """
    Create an iToL shapes dataset showing only presence/absence (not counts)
    Each defense system gets a different colour!
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Dataframe with genome_id as index and subtypes as columns
        Values should be counts (0, 1, 2, etc.) - will be converted to presence/absence
    subtypes : list
        List of column names (subtypes) to include in the dataset
    output_file : str
        Path to output file
    dataset_label : str
        Label for the dataset in iToL
    colors : list or None
        List of hex colors for each subtype. If None, uses default colorblind-friendly palette
    shape_type : int
        Shape type (1=square, 2=circle, 3=star, 4=right triangle, 5=left triangle)
    unique_colors : bool
        If True, each defense system gets its own color (default: True)
    
    Returns:
    --------
    None (writes file)
    """
    
    # Check that all subtypes exist in the dataframe
    missing_subtypes = [st for st in subtypes if st not in df.columns]
    if missing_subtypes:
        print(f"Warning: The following subtypes are not in the dataframe: {missing_subtypes}")
        subtypes = [st for st in subtypes if st in df.columns]
        if not subtypes:
            print("Error: No valid subtypes found!")
            return
    
    print(f"Creating iToL shapes dataset for {len(subtypes)} subtypes:")
    for st in subtypes:
        print(f"  - {st}")
    
    # Generate colorblind-friendly colors if not provided
    if colors is None:
        # Expanded colorblind-friendly palette with distinct colors
        default_colors = [
            '#e41a1c',  # Red
            '#377eb8',  # Blue
            '#4daf4a',  # Green
            '#984ea3',  # Purple
            '#ff7f00',  # Orange
            '#ffff33',  # Yellow
            '#a65628',  # Brown
            '#f781bf',  # Pink
            '#999999',  # Grey
            '#66c2a5',  # Teal
            '#fc8d62',  # Salmon
            '#8da0cb',  # Lavender
            '#e78ac3',  # Magenta
            '#a6d854',  # Lime
            '#ffd92f',  # Gold
            '#e5c494',  # Tan
            '#b3b3b3',  # Silver
            '#8dd3c7',  # Aqua
            '#bebada',  # Light purple
            '#fb8072',  # Light red
        ]
        colors = (default_colors * ((len(subtypes) // len(default_colors)) + 1))[:len(subtypes)]
    
    # Convert to presence/absence (0 or 1)
    data_subset = df[subtypes].copy()
    for st in subtypes:
        data_subset[st] = (data_subset[st] > 0).astype(int)
    
    # Count genomes with each subtype
    print(f"\nPresence/absence statistics:")
    for st in subtypes:
        count = data_subset[st].sum()
        percent = (count / len(data_subset)) * 100
        print(f"  {st}: {count}/{len(data_subset)} genomes ({percent:.1f}%)")
    
    print(f"\nWriting iToL dataset to: {output_file}")
    
    with open(output_file, 'w') as f:
        # Header
        f.write("DATASET_EXTERNALSHAPE\n")
        f.write("#External shapes dataset showing presence/absence of defense systems\n")
        f.write("#Each defense system has a unique color\n")
        f.write("#Shape is shown only when the system is present (value=1)\n\n")
        
        # Mandatory settings
        f.write("SEPARATOR COMMA\n\n")
        f.write(f"DATASET_LABEL,{dataset_label}\n\n")
        f.write("COLOR,#ff0000\n\n")
        
        # Field colors - each defense system gets its own color!
        f.write("FIELD_COLORS," + ",".join(colors) + "\n\n")
        
        # Field labels (subtype names)
        f.write("FIELD_LABELS," + ",".join(subtypes) + "\n\n")
        
        # Optional settings
        f.write("#Optional settings\n")
        f.write("DASHED_LINES,1\n")
        f.write(f"SHAPE_TYPE,{shape_type}\n")
        f.write("COLOR_FILL,1\n")
        f.write("SHOW_VALUES,0\n")  # Don't show numbers, just presence/absence
        f.write("SHOW_LABELS,1\n")
        f.write("HEIGHT_FACTOR,1\n")
        f.write("SHAPE_SPACING,10\n")
        f.write("MARGIN,15\n\n")
        
        # Legend
        f.write("#Legend\n")
        f.write(f"LEGEND_TITLE,{dataset_label}\n")
        f.write("LEGEND_SHAPES," + ",".join([str(shape_type)] * len(subtypes)) + "\n")
        f.write("LEGEND_COLORS," + ",".join(colors) + "\n")
        f.write("LEGEND_LABELS," + ",".join(subtypes) + "\n\n")
        
        # Data section
        f.write("DATA\n")
        
        # Write data for each genome (0 or 1 for each system)
        for genome_id in data_subset.index:
            values = data_subset.loc[genome_id, subtypes].values
            values_str = ",".join(str(int(v)) for v in values)
            f.write(f"{genome_id},{values_str}\n")
    
    print("Done!")
    print(f"\nDataset statistics:")
    print(f"  Total genomes: {len(data_subset)}")
    print(f"  Total defense systems: {len(subtypes)}")
    print(f"  Each system has a unique color!")
    print(f"\nTo use in iToL:")
    print(f"1. Open your tree in iToL")
    print(f"2. Go to 'Datasets' tab")
    print(f"3. Drag and drop '{output_file}'")
    print(f"4. Colored shapes will appear for present systems!")


def create_shapes_dataset_with_counts(df, subtypes, output_file='shapes_dataset.txt',
                                       dataset_label='Defense Systems', colors=None,
                                       shape_type=2):
    """
    Create an iToL shapes dataset showing actual counts (if you want numbers)
    
    Same parameters as create_shapes_dataset_presence_absence
    """
    
    # Check that all subtypes exist in the dataframe
    missing_subtypes = [st for st in subtypes if st not in df.columns]
    if missing_subtypes:
        print(f"Warning: The following subtypes are not in the dataframe: {missing_subtypes}")
        subtypes = [st for st in subtypes if st in df.columns]
        if not subtypes:
            print("Error: No valid subtypes found!")
            return
    
    print(f"Creating iToL shapes dataset for {len(subtypes)} subtypes:")
    for st in subtypes:
        print(f"  - {st}")
    
    # Generate colors if not provided
    if colors is None:
        default_colors = [
            '#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00',
            '#ffff33', '#a65628', '#f781bf', '#999999', '#66c2a5',
            '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f'
        ]
        colors = (default_colors * ((len(subtypes) // len(default_colors)) + 1))[:len(subtypes)]
    
    # Extract data for the selected subtypes
    data_subset = df[subtypes].copy()
    
    # Count genomes with each subtype
    for st in subtypes:
        count = (data_subset[st] > 0).sum()
        print(f"  {st}: {count} genomes")
    
    print(f"\nWriting iToL dataset to: {output_file}")
    
    with open(output_file, 'w') as f:
        # Header
        f.write("DATASET_EXTERNALSHAPE\n")
        f.write("#External shapes dataset showing counts of defense systems\n\n")
        
        # Mandatory settings
        f.write("SEPARATOR COMMA\n\n")
        f.write(f"DATASET_LABEL,{dataset_label}\n\n")
        f.write("COLOR,#ff0000\n\n")
        
        # Field colors
        f.write("FIELD_COLORS," + ",".join(colors) + "\n\n")
        
        # Field labels (subtype names)
        f.write("FIELD_LABELS," + ",".join(subtypes) + "\n\n")
        
        # Optional settings
        f.write("#Optional settings\n")
        f.write("DASHED_LINES,1\n")
        f.write(f"SHAPE_TYPE,{shape_type}\n")
        f.write("COLOR_FILL,1\n")
        f.write("SHOW_VALUES,1\n")  # Show actual counts
        f.write("VALUE_AUTO_COLOR,1\n")
        f.write("SHOW_LABELS,1\n")
        f.write("HEIGHT_FACTOR,1\n")
        f.write("SHAPE_SPACING,10\n\n")
        
        # Legend
        f.write("#Legend\n")
        f.write("LEGEND_TITLE,Defense Systems\n")
        f.write("LEGEND_SHAPES," + ",".join([str(shape_type)] * len(subtypes)) + "\n")
        f.write("LEGEND_COLORS," + ",".join(colors) + "\n")
        f.write("LEGEND_LABELS," + ",".join(subtypes) + "\n\n")
        
        # Data section
        f.write("DATA\n")
        
        # Write data for each genome
        for genome_id in data_subset.index:
            values = data_subset.loc[genome_id, subtypes].values
            values_str = ",".join(str(int(v)) for v in values)
            f.write(f"{genome_id},{values_str}\n")
    
    print("Done!")


# Your defense systems to visualize
subtypes = ['RM_III', 'CAS_Class2-II-A', 'Aditi', 'PDC-S11', 
            'PDC-M22', 'PDC-S30', 'RM_II', 'tmn', 'Ogmios', 'PDC-S04']

# Create the dataset - each system gets a different color!
# Shows only presence/absence (no numbers!)
create_shapes_dataset_presence_absence(
    df=subtype_counts,
    subtypes=subtypes,
    output_file='defense_systems.txt',
    dataset_label='Defense Systems'
)
