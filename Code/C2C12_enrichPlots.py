"""
Code for generating enrichment plots of ShinyGO data for enrichments in genes correlated to Slc7a11 in C2C12 differentiation datasets
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
enrich = pd.read_csv(rf'data\corrs_enrichMolFunc.csv') 	# REPLACE with your directory


"""
Some pathways are identical or redundant, remove these for plotting.
Redundant pathways are identical in terms of total genes (Pathway Genes), genes in our set (nGenes), specific genes in our set that are found in the pathway (Genes), and similar/redundant pathway information
	
	BioProc: ['RNA splicing, via transesterification reactions ', 'RNA splicing, via transesterification reactions with bulged adenosine as nucleop', 'Process utilizing autophagic mechanism ']
	CellCom: ['Fascia adherens ', 'Organellar ribosome ', 'Envelope ']
	MolFunc: ['N-acylsphingosine amidohydrolase activity ']
	KEGG: None, comment out code lines 34-36

Can replace 'redundant_paths' with the arrays above for corresponding enrichments to replicate published enrichment figures

"""
redundant_paths = ['N-acylsphingosine amidohydrolase activity ']
for name in redundant_paths:
	enrich.drop(enrich.loc[enrich['Pathway']==name].index, inplace=True)

enrich.sort_values('Fold Enrichment', inplace=True)
enrich['-log10(FDR)*sign(Pearson)'] = -np.log(enrich.loc[:,'Enrichment FDR']) * np.sign(enrich.loc[:,'Fold Enrichment'])


"""
Actual plotting stuff to plot the enrichment
"""
sns.set_style('whitegrid')
fig, ax = plt.subplots(figsize=(6, 6))
sns.scatterplot(data=enrich, x='Fold Enrichment', y='Pathway', size='nGenes', hue='-log10(FDR)*sign(Pearson)', palette='vlag_r', ax=ax)
plt.xlim(-max(abs(enrich.loc[:,'Fold Enrichment']))-0.5, max(abs(enrich.loc[:,'Fold Enrichment']))+0.5) 	# specify x-axis range

# Now formatting stuff
# By default, legend will display points for size (nGenes) and hue (significance), however we will be using a color bar for hue so only want to include size values
h, l = ax.get_legend_handles_labels()
leg = ax.legend(h[7:], l[7:], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title='Genes', frameon=False) 	# Note: For KEGG plots, use h[8:], l[8:]

ax.axvline(0, linewidth=1.5, color='gray') 	# Grey line at 0 for reference
ax.xaxis.grid(False)


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
sm = plt.cm.ScalarMappable(cmap='vlag_r', norm=norm)
sm.set_array([])

# Create axis object for colorbar and format
cax = fig.add_axes([ax.get_position().x1+0.055, ax.get_position().y0+0.22, 0.075, ax.get_position().height/3])
ax.figure.colorbar(sm, fraction=0.04, cax=cax)
cax.set_xlabel('-log10(FDR)*sign(Pearson)', horizontalalignment='left', x=-0.45)
cax.xaxis.set_label_position('top') 

# Save figure
plt.savefig(f'enrichment_MF.tiff', bbox_inches='tight', pad_inches=0.9, dpi=1200)
