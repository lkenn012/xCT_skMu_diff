"""
Code for generating enrichment plots of ShinyGO data for enrichments in genes correlated to Slc7a11 in C2C12 differentiation datasets as reported in Kanaan et al (2023)
	- Loads a data file contain ShinyGO ouputs for both positive and negative correlated gene sets
	- Transforms values (i.e. FDR values) for better visualization
	- Plots and formats the figures
	- Save figure

"""

# import modules
import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt


"""
load our enrichment data, modified from ShinyGO (such that negative correlations are denoted by negative fold enrichments)
	Biological process: corrs_enrichBioProc.csv
	Cellular Component: corrs_enrichCellCom.csv
	Molecular Function: corrs_enrichMolFunc.csv
	KEGG: corrs_enrichKEGG.csv
"""
enrich = pd.read_csv(rf'data\corrs_enrichKEGG.csv') 	# REPLACE with your directory


"""
Some pathways are identical or redundant, remove these for plotting.
Redundant pathways are identical in terms of total genes (Pathway Genes), genes in our set (nGenes), specific genes in our set that are found in the pathway (Genes), and similar/redundant pathway information
	
	BioProc: ['RNA splicing, via transesterification reactions ', 'RNA splicing, via transesterification reactions with bulged adenosine as nucleop', 'Process utilizing autophagic mechanism ']
	CellCom: ['Fascia adherens ', 'Organellar ribosome ', 'Envelope ']
	MolFunc: ['N-acylsphingosine amidohydrolase activity ']
	KEGG: None, comment out code lines 34-36

Can replace 'redundant_paths' with the arrays above for corresponding enrichments to replicate published enrichment figures

"""
# redundant_paths = ['RNA splicing, via transesterification reactions ', 'RNA splicing, via transesterification reactions with bulged adenosine as nucleop', 'Process utilizing autophagic mechanism ']
# for name in redundant_paths:
# 	enrich.drop(enrich.loc[enrich['Pathway']==name].index, inplace=True)

enrich.sort_values('Fold Enrichment', inplace=True)
enrich['-log10(FDR)*sign(Pearson)'] = -np.log(enrich.loc[:,'Enrichment FDR']) * np.sign(enrich.loc[:,'Fold Enrichment'])


"""
Actual plotting stuff to plot the enrichment
"""
sns.set_style('whitegrid')
sns.set_context("paper", rc={"font.size":8,"axes.titlesize":16})
plt.rcParams['font.family'] = 'Arial'
fig, ax = plt.subplots(figsize=(6, 6))

sns.scatterplot(data=enrich, 
	 x='Fold Enrichment',
	 y='Pathway',
	 size='nGenes', 
	 sizes=(min(enrich.loc[:,'nGenes'])*5, max(enrich.loc[:,'nGenes'])*0.75), 	# Range of size to map {size} to values in our data, adjust scaling so larger values do not overlap. KEGG: [min*5, max*0.75]
	 hue='-log10(FDR)*sign(Pearson)', 
	 palette='RdBu_r', 
	 edgecolor='black', 	# edge color of points
	 ax=ax
	 )

plt.xlim(-max(abs(enrich.loc[:,'Fold Enrichment']))-0.5, max(abs(enrich.loc[:,'Fold Enrichment']))+0.5) 	# specify x-axis range
plt.title('KEGG Enrichment', fontsize=14) 	# NOTE: replace with title of enrichment

# Now formatting stuff
# By default, legend will display points for size (nGenes) and hue (significance), however we will be using a color bar for hue so only want to include size values
h, l = ax.get_legend_handles_labels()
leg = ax.legend(h[8:], l[8:], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title='Genes', frameon=False) 	# Note: For KEGG plots, use h[8:], l[8:]

ax.axvline(0, linewidth=2, color='gray') 	# Gray line at 0 for reference
ax.xaxis.grid(False)
plt.tick_params(axis='x', bottom=True, color='gray') 	# Add ticks at x-axis
ax.invert_yaxis()

"""
For continuous FDR values, need to create a colorbar for the legend
"""

# Define limits of the color bar based on largest absolute values
if enrich.loc[:,'-log10(FDR)*sign(Pearson)'].abs().max() >= enrich.loc[:,'-log10(FDR)*sign(Pearson)'].abs().min():
	color_lim = enrich.loc[:,'-log10(FDR)*sign(Pearson)'].abs().max()
else:
	color_lim = enrich.loc[:,'-log10(FDR)*sign(Pearson)'].abs().min()

# Create colorbar map
norm = plt.Normalize(color_lim, -color_lim)
sm = plt.cm.ScalarMappable(cmap='RdBu_r', norm=norm)
sm.set_array([])

# Create axis object for colorbar and format
cax = fig.add_axes([ax.get_position().x1+0.055, ax.get_position().y0+0.22, 0.075, ax.get_position().height/3])
cb = ax.figure.colorbar(sm, fraction=0.04, cax=cax) #, drawedges=True)
cb.outline.set_color('black')
cax.tick_params(size=0)
# cax.set_color('black')
cax.set_xlabel('-log10(FDR)*sign(Pearson)', horizontalalignment='left', x=-0.45)
cax.xaxis.set_label_position('top') 

# plt.xticks(ax.get_xticklabels()) 	# Add x-axis ticks
ax.spines['bottom'].set_color('gray') 	# Adjust bottom border to anchor the figure
ax.spines['bottom'].set_linewidth(2)

# Save figure
plt.savefig(f'enrichment_KEGG_new.tiff', bbox_inches='tight', pad_inches=0.9, dpi=300)
