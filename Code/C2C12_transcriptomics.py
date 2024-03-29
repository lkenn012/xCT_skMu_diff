"""
Code for producing analysis and figures of C2C12 transcriptomics data as reported in Kanaan et al (2023)
"""

# import modules
import pandas as pd
import numpy as np
import collections
import re
from sklearn.preprocessing import StandardScaler

import matplotlib.pyplot as plt
import seaborn as sns

import scipy.cluster.hierarchy as hierarch
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr, t
from statsmodels.stats.dist_dependence_measures import distance_correlation

import glob
from datetime import datetime


# define function to access all our files of interest in a given data folder
def get_fNames(name_match, dir_path):

	"""
	name_match: common file name pattern for expression data files
	dir_path: path directory to the data files
	"""

	# this should access all files in our data folder which have {*Some wildcards}+name_match file name
	if dir_path == "" or dir_path is None:
		fs = [f for f in glob.glob(rf'\*{name_match}')] 	# if no directory is given, check current
	else:
		fs = [f for f in glob.glob(rf'{dir_path}\*{name_match}')] 	# else, check the specified directory

	return fs


# define function to get z-score normalization
def z_norm(df):

	"""
	df: dataframe to normalize, assumes rows are features (genes) and columns are observations 
	"""

	scaler = StandardScaler() 	# z-score normalization

	# need to pop and then re-insert first row, since that contains time data
	df_t = df.T
	times = df_t.pop('Time') # remove 'Time' data for normalization
	z_df = pd.DataFrame(data=scaler.fit_transform(df_t), index=df_t.index, columns=df_t.columns) 	# scaled df

	z = z_df.T 	# rows as genes
	z.loc['Time'] = times

	return z

# define function to calculate log2 fold change of gene expression from t-0
def fold_change(expression_df):

	"""
	expression_df: gene expression data for normalization, assumes rows are genes
	"""

	# To get the fold change, need to get the numer of replicates and then divide each replicates' data by its t0 value
	min_Time = expression_df.loc['Time'].min()
	n = len(expression_df.loc[:,expression_df.loc['Time'] == min_Time].columns)
	n_times = len(set(expression_df.loc['Time'].tolist())) 	# get number of times for iterating over

	# now, for each replicate must calculate the fold change
	minT_exprDF = expression_df.loc[:, expression_df.loc['Time'] == min_Time]
	expr_df = expression_df
	times = expr_df.loc['Time'] # remove 'Time' data for normalization
	expr_df.drop('Time', axis=0, inplace=True)
	minT_exprDF = minT_exprDF.drop('Time', axis=0) 

	minT_exprDF = minT_exprDF + 0.000001 	# add small constants to avoid 0/0
	expr_df = expr_df + 0.000001

	# normalize each replicate data by its respective value
	for i in range(n):
		idx = i
		t0 = minT_exprDF.iloc[:,i]
		for j in range(n_times):
			try:
				expr_df.iloc[:,idx] = expr_df.iloc[:,idx] / t0 	# set column values to expression / t0 expression
			except IndexError as e:
				print(f'{str(e)}')
				print(f'out of ranges, possibly inconsistent replicates:\n{expr_df.columns}') 	# need to catch any inconsistent replicates (i.e. less replicates for a time points)
			idx += n 	# update the idx

	FC_df = pd.DataFrame(data=np.log2(expr_df), index=expr_df.index, columns=expr_df.columns) 	# log2 df 
	FC_df.loc['Time'] = times 	# time row	

	return FC_df


# define function for calculating p-values for dcorr values according to a two-sided P-value for 
# the t-distribution with n-2 degrees of freedom (https://doi.org/10.1093/bioinformatics/btad210)
def p_val(x, n):

	"""
	x = the sample value, a dcorr correlation value
	n = our sample size, the number of dcorr values in our matrix
	"""

	# Compute T statistic and corresponding p-value from cumulative distribution function
	try:
		t_stat = x*(n-2)**0.5/(1-x**2)**0.5
		p = 2*(1-t.cdf(t_stat, n-2)) 	# two-sided p-value
	except ZeroDivisionError:
		p = 1

	return p

# define recursive function for flattening a list which contains lists and int values such that the return is a list of all elements
def flatten(xs):
	"""
	xs = a list of lists and values
	"""
	if isinstance(xs, collections.Iterable):
		return [a for x in xs for a in flatten(x)]
	else:
		return [xs]


# recursive function to call in get_nearest to get original values from clusters in a linkage matrix
def get_child(Z, i, n, idxs=False):

	"""
	Z = a scipy (or sklearn?) linkage matrix where Z[i] = [node_j, node_k, distance, n_values]
	i = index where we need to get child leaves
	n = number of original values, will inform what is a new cluster vs. original value
	idxs = Optional, list of linkage matrix idxs to append child idxs to
	"""
	# Get our index values to return
	val_0 = Z[int(i)][0]
	val_1 = Z[int(i)][1]
	idx_0 = int(i)
	idx_1 = int(i)

	# Check if these values are original values or clusters (where val > n)
	if val_0 > n:
		idx_0, val_0 = get_child(Z, val_0-n, n, idxs)
	if val_1 > n:
		idx_1, val_1 = get_child(Z, val_1-n, n, idxs)
	else:
		pass

	# return the indexes obtained over recursively getting clusters until original values
	if idxs:
		idxs = flatten([idx_0, idx_1])
		return idxs, [val_0, val_1]

	else:
		return [val_0, val_1]


# Define function to get the nearest neighbours of a given observation in a hierarchical clustering linkage map
def get_nearest(Z, node_index, max_steps=False, max_dist=False, get_links=False):

	"""
	Z = a scipy (or sklearn?) linkage matrix where Z[i] = [node_j, node_k, distance, n_values]
		is the new linkage formed at step_i in the clustering algorithm, where a node may be an original value or a new cluster if node_index > n (original samples) and node index = n+i
	node_index = value index of either the original value we are interested in, or a new cluster formed by the algorithm if node_index > n
	max_steps = thresholding parameter for the number of clustering steps to extend our nearest neighbour search for
	max_dist = thresholding parameter for a distance cut-off to threshold the number of clusters we obtain
	sub_links = whether to get and return the sub linkage matrix, 
	"""

	n = np.shape(Z)[0] + 1 	# number of original observations/values, corresponds to rows + 1 in Z
	links_idx = [] 	# list of all linkage matrix indexes for generating a sub matrix if get_links=True
	steps = 0 	# Thresholding checks
	last_i = 0

	# Find the first clustering of our node of interest and get all original values in this clustering (using get_child())
	for i, l in enumerate(Z):
		# First, iterate over the matrix and work "up" the linkage matrix while checking for our thresholding parameters
		if l[0] == node_index or l[1] == node_index:
			links_idx.append(i)

			if max_steps and steps < max_steps: 	# Check that we have not reached max steps
				steps += 1
				node_index = i + n

			elif max_dist and l[2] < max_dist: 	# Check that we have not surpassed max distance
				last_i = i
				node_index = i + n

			else:
				if max_dist:
					l = Z[last_i] 	# For max_dist threshold, we check that our cluster distance is less than thershold, when False we must proceed with the last cluster under the distance threshold
					del links_idx[-1]
				# once steps = max_steps we can work back through the linkage matrix to get all the original values contained within clusters
				if l[0] > n:
					# print(f'[l[0]] > n')
					idxs_0, val_0 = get_child(Z, l[0]-n, n, links_idx)
				else:
					val_0 = l[0]

				if l[1] > n:
					# print(f'[l[1]] > n')
					idxs_1, val_1 = get_child(Z, l[1]-n, n, links_idx)
				else:
					val_1 = l[1]

	# convert val_0 and val_1 from lists of lists to flattened list for easy indexing
	vals = [x for xs in [val_0, val_1] for x in xs]

	# now can return the final list of values
	if get_links:
		idxs = [x for xs in [links_idx, idxs_0, idxs_1] for x in xs]
		idxs = sorted(list(set(idxs)))
		return Z[idxs], vals
	else:
		return vals


"""
Define main function which will:
	load GSE data files
	pre-process data:
		convert nan expression values to 0
		remove any genes with 0 expression across t
		remove any missing or duplicate genes (based on ID)
	Processing for correlation analysis:
		Standardize expression values
		Select datasets
	Correlate standardized gene expression
	Plot expression data for genes of interest
"""
def main():

	"""
	parameters to define the analysis, i.e., how to pre-process data and whether to compute correlations
	Note: correlations & clustering are computed with different normalizations than values used in expression plot
	Thus, visualization and correlations will differ from reported results if only one normalization is used for both 
	Commenting out either section for different normalization will avoid this.
	"""

	remove_nan = True 	# replace NAN gene expression with very small value (necessary for some plotting and analysis)
	scale = True 	# whether to scale the data, used in correlation & clustering analysis
	FC = False		# Specify fold change
	compute_corrs = False 	# Specify computing correlations
	compute_clust = True  	# Specify computing clustering
	expr_threshold = (-8,8)		# values to threshold expression Fold Change data for plotting for where extreme values arise

	# Define variables for accessing our data 
	f_match = 'allExpression_data.csv'  	# REPLACE with file suffix of choice (e.g. "{GEO_id}_allExpression_data.csv")
	f_names = get_fNames(name_match=f_match, dir_path='data') 	# specifiy the file names we want to use REPLACE dir_path with path string, or empty string if same directory
	date = datetime.today().strftime('%d-%m-%Y')
	out_name = f'{date}_log2FC' 	# name for output files


	# Load data files as dfs
	dfs = []
	data_names = []
	for name in f_names:

		df = pd.read_csv(rf'{name}', index_col=0, header=None)
		print(f'initial df data (file: {name}):\n{df}')
		data_names.append(re.search(r'GSE[0-9]+', name).group()) 	# name of dataset for naming outputs later


		# Some formatting
		df.loc['Time'] = df.iloc[0] 	# creates a named row of the times
		df = df.iloc[1:] 	# drop the unnamed 'Time' row
		df = df.loc[:, df.loc['Time'] >= 0] 	# Remove any timepoints before differentiation, due to variations in treatments

		df.rename({np.nan: 'missing_genes'}, inplace=True) 	# rename row which contains missing gene symbols 

		# Now need to clean the data
		try:
			df.drop('missing_genes', axis=0, inplace=True) 	# remove any 'missing_gene' rows

		except KeyError:
			pass

		# Check for number of time points, if we have only one time point then cannot conduct any further analyses
		if len(set(df.loc["Time"])) < 2:
			print(f'Only two time points, skipped')
			continue 


		# Remove missing values, replace with 0 expression
		# Optionally, can imput missing values
		if remove_nan:
			n_missing = df.isnull().sum().sum()
			print(f'expression data missing values: {n_missing}')
			df.fillna(0, inplace=True)

		# Scale the data using z-score transformation, if desired
		if scale:
			df = z_norm(df)

		# May want to, instead use Fold change
		elif FC:
			df = fold_change(df) 	# convert data to foldchange (vs. avg. t0 expressions). When log2 is True, instead calculate log2 foldchange

		df = df.loc[~(df==0).all(axis=1)] 	# remove any rows =0 (i.e., no expression)
		df = df.loc[~df.index.duplicated()] 	# remove any dulpicate genes in data

		dfs.append(df)

	if compute_corrs:
		corr_dfs = []
		for i, df in enumerate(dfs):

			# Generating correlations between 2 time points, not relevant
			if len(set(df.loc["Time"])) <= 2:
				print(f'Only two time points, skipped')
				continue

			else:
				slc7a11_expr = df.loc['Slc7a11']
				df.drop('Time', inplace=True)

				dcorrs = [distance_correlation(x=slc7a11_expr, y=df.iloc[i]) for i in range(len(df))]
				p_vals = [p_val(x=x, n=len(df.columns)) for x in dcorrs]
				pears = [pearsonr(x=slc7a11_expr, y=df.iloc[i])[0] for i in range(len(df))]

				corrs_df = pd.DataFrame(data=[pears,dcorrs, p_vals], index=[f'{data_names[i]} Pears', f'{data_names[i]} Dcorr', f'{data_names[i]} p-value'], columns=df.index)
				corr_dfs.append(corrs_df)

		# Save all correlation results
		out_df = pd.concat(corr_dfs, axis=0).T
		print(f'out_df:\n{out_df}')


		out_df.to_csv(rf'C2C12_correlations.csv')

	if compute_clust:
		clust_dfs = []
		for i, df in enumerate(dfs):

			# Generating correlations between 2 time points, not relevant
			if len(set(df.loc["Time"])) <= 2:
				print(f'Only two time points, skipped')
				continue

			else:
				# Calculate distance matrix, default metric = euclidean, calculates distance between rows
				sample_dist = pdist(df, 'cosine') 	
				sample_links = hierarch.average(sample_dist) 	# Calculate linkages between samples using ward metric

				# Interested in obtain the nearest neighbours to Slc7a11
				sub_links, nearest_genes = get_nearest(Z=sample_links, node_index=df.index.get_loc('Slc7a11'), max_dist=0.98, get_links=True)
				nearest_genes = flatten(nearest_genes)

				sample_dist = pdist(df.iloc[nearest_genes], 'cosine') 	# Calculate distance matrix, default metric = euclidean, calculates distance between rows
				sample_links = hierarch.average(sample_dist) 	# Calculate linkages between samples using ward metric

				# get nearest genes distances
				sq_dist = squareform(sample_dist + sample_dist.T)
				slc7a11_dist = pd.DataFrame(data=sq_dist, index=df.iloc[nearest_genes].index, columns=df.iloc[nearest_genes].index)
				slc7a11_dist = slc7a11_dist.loc[:,'Slc7a11']
				slc7a11_dist.to_csv(rf'{data_names[i]}_clusters.csv')
				clust_dfs.append([nearest_genes, slc7a11_dist, sample_links])  	# Append to list of cluster data for any further analysis


	###
	## 	Code for generating line plot of Genes with Slc7a11 & average lines
	###

	# now combine our data so we have all samples as rows
	expr_df = pd.concat(dfs, axis=1) 
	expr_df.columns = pd.RangeIndex(expr_df.columns.size) # Column headers repeat, so simply reset these

	expr_df.rename({np.nan: 'Time'}, inplace=True) 	# rename row which contains time points

	# May need to expr_threshold data to remove extreme outliers
	if expr_threshold:
		temp = expr_df.loc[expr_df.index != 'Time'] 	# do not want to threshold time values

		# Threshold our expression values such that log2FC are within our threshold bounds
		temp = temp.applymap(lambda x: expr_threshold[0] if x < expr_threshold[0] else (expr_threshold[1] if x > expr_threshold[1] else x))
		temp.loc['Time'] = expr_df.loc['Time']
		expr_df = temp

	### Want to visualize the expression of Slc7a11 and relevant genes for differentation, glutathione metabolism over differentiation
	gene_sets = {'Slc7a11': [r'Slc7a11$', r'Slc7a11.\d$'], 	# regex patterns for exact probe matches
		'Myh1': [r'Myh1$', r'Myh1.\d$']
		}

	sub_expr = []
	group_ids = []
	for key, item in gene_sets.items():
		for i in item:
			temp = expr_df.loc[expr_df.index.str.contains(i, regex=True)]
			sub_expr.append(temp) 	# get rows which contain our gene symbol (could be GENE.1, GENE.2, etc. due to probes)

			group_ids.extend([key]*len(temp))

	subExpr_df = pd.concat(sub_expr)
	
	# Add time column, convert to hours and round up for plotting
	times = expr_df.loc['Time'] * 24
	times = times.round(1)
	times.replace([0.5, 1.5], [1.0, 2.0], inplace=True)
	subExpr_df.columns = times

	# Only interested in time points within [112 hours, 8 days (192hours)], remove all beyond this
	subExpr_df = subExpr_df.drop(subExpr_df.loc[:,(subExpr_df.columns < 12) & (subExpr_df.columns > 0)].columns, axis=1)
	subExpr_df = subExpr_df.drop(subExpr_df.loc[:,(subExpr_df.columns > 192)].columns, axis=1)

	subExpr_df.loc[:, 'Gene'] = group_ids


	"""
	Define variables for plotting expression data
	"""
	plt.figure(figsize=(8,6))
	sns.set_style('white')
	fig,ax = plt.subplots()

	melt = subExpr_df.melt(id_vars='Gene', var_name='Time', value_name='expr') 	# Reshape data for plotting

	sns.set_palette(['#0F99B2', '#808080'])
	sns.lineplot(data=melt, x='Time', y='expr', markers=['o', 'o'], style='Gene', dashes=False, hue='Gene', estimator='mean', errorbar=('ci', 95), ax=ax)

	# Some formatting
	ax.axhline(0, ls='--', linewidth=1, color='black')
	sns.despine(offset=0, trim=True) 	# adjusts the surrounding plot box
	plt.xticks(np.arange(min(melt['Time']), max(melt['Time'])+1, step=24))
	plt.xlabel('Time (hours)')
	plt.ylabel('Gene expression (log2 Fold change)')


	# Save figure
	plt.tight_layout()
	plt.savefig(f'C2C12_geneExpression_linePlot.tiff', dpi=1200)
	plt.close()

# run main code
main()