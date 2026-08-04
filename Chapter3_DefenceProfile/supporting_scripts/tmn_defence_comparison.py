"""
SIMPLE VERSION - Adapt to your dataframe structure
Copy this into your Jupyter notebook and modify variable names
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# ============================================================================
# SETUP: Define your data
# ============================================================================

# MODIFY THESE LINES based on your actual dataframes:
# defense_df = your defense dataframe (genome_id + defense system columns)
# tmn_large_ids = list of genome IDs with Large Tmn cluster
# tmn_small_ids = list of genome IDs with Small Tmn cluster

# Example if you have the cluster files:
tmn_large_ids = [line.strip() for line in open('LargeCluster_Tmngenes.txt')]
tmn_small_ids = [line.strip() for line in open('SmallCluster_Tmngenes.txt')]

# Extract genome IDs (remove the @protein_id part)
tmn_large_genomes = list(set([id.split('@')[0] for id in tmn_large_ids]))
tmn_small_genomes = list(set([id.split('@')[0] for id in tmn_small_ids]))

print(f"Genomes with Large Tmn: {len(tmn_large_genomes)}")
print(f"Genomes with Small Tmn: {len(tmn_small_genomes)}")

# ============================================================================
# Add Tmn status to your defense dataframe
# ============================================================================

# Assuming defense_df has a 'genome_id' column
defense_df['tmn_status'] = 'No_Tmn'
defense_df.loc[defense_df['genome_id'].isin(tmn_large_genomes), 'tmn_status'] = 'Large'
defense_df.loc[defense_df['genome_id'].isin(tmn_small_genomes), 'tmn_status'] = 'Small'
defense_df['has_tmn'] = defense_df['tmn_status'] != 'No_Tmn'

print(f"\nNo Tmn: {(defense_df['tmn_status'] == 'No_Tmn').sum()}")
print(f"Large Tmn: {(defense_df['tmn_status'] == 'Large').sum()}")
print(f"Small Tmn: {(defense_df['tmn_status'] == 'Small').sum()}")

# ============================================================================
# Calculate total defense genes
# ============================================================================

# MODIFY THIS: List your defense system columns
# Example: defense_columns = ['RM_I', 'RM_II', 'CRISPR', 'Abi', 'CBASS', etc.]
defense_columns = [col for col in defense_df.columns 
                   if col not in ['genome_id', 'Strain', 'tmn_status', 'has_tmn']]

print(f"\nDefense columns: {defense_columns[:10]}...")  # Show first 10

# Calculate totals
defense_df['total_defense'] = defense_df[defense_columns].sum(axis=1)
defense_df['num_systems'] = (defense_df[defense_columns] > 0).sum(axis=1)

# ============================================================================
# STATISTICAL TESTS
# ============================================================================

print("\n" + "="*70)
print("WITH TMN vs WITHOUT TMN")
print("="*70)

with_tmn = defense_df[defense_df['has_tmn']]['total_defense']
without_tmn = defense_df[~defense_df['has_tmn']]['total_defense']

print(f"\nWith Tmn (n={len(with_tmn)}): {with_tmn.mean():.1f} ± {with_tmn.std():.1f}")
print(f"Without Tmn (n={len(without_tmn)}): {without_tmn.mean():.1f} ± {without_tmn.std():.1f}")
print(f"Difference: {with_tmn.mean() - without_tmn.mean():.1f} genes")

stat, pval = stats.mannwhitneyu(with_tmn, without_tmn)
print(f"\nMann-Whitney U: p = {pval:.4e} {'***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'}")

print("\n" + "="*70)
print("LARGE vs SMALL vs NO TMN")
print("="*70)

for status in ['No_Tmn', 'Large', 'Small']:
    subset = defense_df[defense_df['tmn_status'] == status]['total_defense']
    print(f"{status:10} (n={len(subset):4}): {subset.mean():5.1f} ± {subset.std():4.1f}")

# Kruskal-Wallis test
groups = [defense_df[defense_df['tmn_status'] == s]['total_defense'] 
          for s in ['No_Tmn', 'Large', 'Small']]
stat, pval = stats.kruskal(*groups)
print(f"\nKruskal-Wallis: p = {pval:.4e} {'***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'}")

# Pairwise comparisons
print("\nPairwise comparisons:")
large = defense_df[defense_df['tmn_status'] == 'Large']['total_defense']
small = defense_df[defense_df['tmn_status'] == 'Small']['total_defense']
stat, pval = stats.mannwhitneyu(large, small)
print(f"  Large vs Small: p = {pval:.4e}")

stat, pval = stats.mannwhitneyu(large, without_tmn)
print(f"  Large vs No_Tmn: p = {pval:.4e}")

stat, pval = stats.mannwhitneyu(small, without_tmn)
print(f"  Small vs No_Tmn: p = {pval:.4e}")

# ============================================================================
# VISUALIZATIONS
# ============================================================================

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Plot 1: Boxplot
ax = axes[0]
defense_df.boxplot(column='total_defense', by='tmn_status', ax=ax)
ax.set_title('Total Defense Genes by Tmn Status')
ax.set_xlabel('Tmn Status')
ax.set_ylabel('Total Defense Genes')
plt.sca(ax)
plt.xticks(rotation=0)

# Plot 2: Violin plot
ax = axes[1]
sns.violinplot(data=defense_df, x='tmn_status', y='total_defense', ax=ax, 
               order=['No_Tmn', 'Large', 'Small'])
ax.set_title('Distribution of Defense Genes')
ax.set_xlabel('Tmn Status')
ax.set_ylabel('Total Defense Genes')

# Plot 3: Bar plot with individual points
ax = axes[2]
order = ['No_Tmn', 'Large', 'Small']
positions = range(len(order))

means = [defense_df[defense_df['tmn_status']==s]['total_defense'].mean() for s in order]
sems = [defense_df[defense_df['tmn_status']==s]['total_defense'].sem() for s in order]

ax.bar(positions, means, yerr=sems, capsize=5, alpha=0.7, 
       color=['gray', 'red', 'blue'])
ax.set_xticks(positions)
ax.set_xticklabels(order)
ax.set_xlabel('Tmn Status')
ax.set_ylabel('Mean Total Defense Genes')
ax.set_title('Mean ± SEM')

plt.tight_layout()
plt.savefig('tmn_defense_comparison.png', dpi=300, bbox_inches='tight')
plt.show()

# ============================================================================
# DEFENSE SYSTEM ENRICHMENT
# ============================================================================

print("\n" + "="*70)
print("DEFENSE SYSTEM ENRICHMENT")
print("="*70)

enrichment = []
for col in defense_columns:
    no_tmn_prev = (defense_df[defense_df['tmn_status']=='No_Tmn'][col] > 0).mean() * 100
    large_prev = (defense_df[defense_df['tmn_status']=='Large'][col] > 0).mean() * 100
    small_prev = (defense_df[defense_df['tmn_status']=='Small'][col] > 0).mean() * 100
    
    enrichment.append({
        'System': col,
        'No_Tmn': no_tmn_prev,
        'Large': large_prev,
        'Small': small_prev,
        'Large_diff': large_prev - no_tmn_prev,
        'Small_diff': small_prev - no_tmn_prev
    })

enrich_df = pd.DataFrame(enrichment).sort_values('Large_diff', ascending=False)

print("\nTop 10 systems enriched in Large Tmn:")
print(enrich_df.head(10)[['System', 'No_Tmn', 'Large', 'Large_diff']])

print("\nTop 10 systems enriched in Small Tmn:")
enrich_df_small = enrich_df.sort_values('Small_diff', ascending=False)
print(enrich_df_small.head(10)[['System', 'No_Tmn', 'Small', 'Small_diff']])

# Save results
defense_df[['genome_id', 'tmn_status', 'total_defense', 'num_systems']].to_csv('tmn_defense_summary.csv', index=False)
enrich_df.to_csv('defense_enrichment.csv', index=False)

print("\n✓ Saved: tmn_defense_summary.csv")
print("✓ Saved: defense_enrichment.csv")
print("✓ Saved: tmn_defense_comparison.png")
