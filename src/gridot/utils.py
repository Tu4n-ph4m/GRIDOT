"""adapted from GrID-Net
"""
from .train import *

import numpy as np
import pandas as pd
import os
from anndata import AnnData
from scipy.stats import f
from scipy.sparse import csr_matrix
import scanpy as sc
import scanpy.external as sce
from scipy.stats import spearmanr
from schema import SchemaQP
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.neighbors import KNeighborsClassifier
import torch
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import issparse




def gene_info(x):
    # Extract gene names
#     print(x)
    gene_name = list(filter(lambda x: 'gene_name' in x,  x.split(";")))[0].split("=")[1]
    gene_type = list(filter(lambda x: 'gene_type' in x,  x.split(";")))[0].split("=")[1]
    level = int(list(filter(lambda x: 'level' in x,  x.split(";")))[0].split("=")[1])
    return (gene_name, gene_type, level)

def update_anndata_obs(rna_adata,gencode_genes):
    merged_df = rna_adata.var.merge(gencode_genes, left_on='features', right_on='gene_name',how = 'left')
    index = list(merged_df.features)
    merged_df.index =index
#     print(merged_df)
    merged_df.rename(columns={
    'start': 'txstart',
    'seqname': 'chr_no',
}, inplace=True)

    merged_df = merged_df.loc[:,['features','chr_no','txstart']]
    return merged_df

def gencode_reading(gff_file,gene_list):
    gencode = pd.read_table(gff_file, comment="#",
                    sep = "\t", names = ['seqname', 'source', 'feature', 'start' , 'end', 'score',
                                         'strand', 'frame', 'attribute'])
    gencode_genes = gencode[(gencode.feature == "gene")][['seqname', 'start', 'end', 'attribute']].copy().reset_index().drop('index', axis=1) # Extract genes

    gencode_genes["gene_name"], gencode_genes["gene_type"], gencode_genes["gene_level"] = zip(*gencode_genes.attribute.apply(lambda x: gene_info(x)))
    gencode_genes = gencode_genes.drop_duplicates(subset='gene_name', keep='first')
#### remove chM
    gencode_genes = gencode_genes.loc[gencode_genes.gene_name.isin(gene_list)]
    gencode_genes = gencode_genes[~gencode_genes['seqname'].str.contains('M')]
    gencode_genes = gencode_genes.loc[gencode_genes.gene_name.isin(gene_list)]
    # gencode_genes['seqname'] = gencode_genes['seqname'].str.replace('chr', '')

    return gencode_genes,gencode_genes['gene_name'].to_list()

def preprocessing_adata_celltype(atac_anndata, rna_anndata,gene_file,celltype):
    print(f"Pre-processing scATAC-seq...")
    atac_adata = atac_anndata.copy()
    rna_adata = rna_anndata.copy()

    # Find the intersection of indices
    common_indices = atac_adata.obs.index.intersection(rna_adata.obs.index)

    # Subset both AnnData objects to keep only the common indices
    atac_adata = atac_adata[common_indices]
    rna_adata = rna_adata[common_indices]

    atac_adata = atac_adata[atac_adata.obs['label']==celltype]
    sc.pp.filter_cells(atac_adata, min_genes=100)
    sc.pp.filter_genes(atac_adata, min_cells=3)
    print(atac_adata)
    
    print(f"Pre-processing scRNA-seq...")
    rna_adata.var['features'] = rna_adata.var.index
    rna_adata = rna_adata[rna_adata.obs['label']==celltype]
    
    sc.pp.filter_cells(rna_adata, min_genes=100)
    sc.pp.filter_genes(rna_adata, min_cells=3)
    print(rna_adata)
    
    print('Adding start site to meta file...')
    gencode_genes,filtered_genes = gencode_reading(gene_file,rna_adata.var['features'])
    rna_adata = rna_adata[:, filtered_genes]
    rna_adata.var = update_anndata_obs(rna_adata,gencode_genes)
    
    return atac_adata, rna_adata


def preprocessing_adata(atac_anndata, rna_anndata,gene_file):
    print(f"Pre-processing scATAC-seq...")
    atac_adata = atac_anndata.copy()
    rna_adata = rna_anndata.copy()

    # Find the intersection of indices
    common_indices = atac_adata.obs.index.intersection(rna_adata.obs.index)

    # Subset both AnnData objects to keep only the common indices
    atac_adata = atac_adata[common_indices]
    rna_adata = rna_adata[common_indices]

#     atac_adata = atac_adata[atac_adata.obs['label']==celltype]
    sc.pp.filter_cells(atac_adata, min_genes=100)
    sc.pp.filter_genes(atac_adata, min_cells=3)
    print(atac_adata)
    
    print(f"Pre-processing scRNA-seq...")
    rna_adata.var['features'] = rna_adata.var.index
#     rna_adata = rna_adata[rna_adata.obs['label']==celltype]
    
    sc.pp.filter_cells(rna_adata, min_genes=100)
    sc.pp.filter_genes(rna_adata, min_cells=3)
    print(rna_adata)
    
    print('Adding start site to meta file...')
    gencode_genes,filtered_genes = gencode_reading(gene_file,rna_adata.var['features'])
    rna_adata = rna_adata[:, filtered_genes]
    rna_adata.var = update_anndata_obs(rna_adata,gencode_genes)
    
    return atac_adata, rna_adata


def n_hop_adjacency(adj_matrix, n):
    if n == 1:
        return np.array(adj_matrix)
    else:
        return np.dot(adj_matrix, n_hop_adjacency(adj_matrix, n - 1))


import os

def create_folder(path):
    try:
        os.makedirs(path)
        print(f"Folder created at: {path}")
    except FileExistsError:
        print(f"Folder already exists at: {path}")
    except Exception as e:
        print(f"An error occurred: {e}")


def cell_root(rna_adata,marker):
    if marker in rna_adata.var_names :
    # Get the expression values for TOP2A
        marker_gene_expression = rna_adata[:, marker].X

        # Find the index of the maximum expression value
        max_index = marker_gene_expression.argmax()

        # Retrieve the cell name and the maximum expression value
        max_expression_cell = rna_adata.obs_names[max_index]
        max_expression_value = marker_gene_expression[max_index]

        print(f"Cell with highest {marker} expression: {max_expression_cell,max_index}")
        print(f"Expression level: {max_expression_value}")
    else:
        print(f"Gene {marker} not found in the dataset.")
        max_index = 0
    return max_index

def pseudo_atac_generation(original_atac,coupling_matrix,n,threshold=0.1):
    '''
    original_atac: atac raw counts in sparse matrix form
    coupling_matrix: coupling matrix in tensor form
    threshold: binarization threshold
    n: number of cells that contribute to pseudo-bulking approximation
    '''
    # Sort each row in descending order
    sorted_tensor, _ = torch.sort(coupling_matrix, descending=True, dim=1)
    # Get the threshold (n-th largest value) for each row
    thresholds = sorted_tensor[:, n - 1].unsqueeze(1)
    # Keep values >= threshold, set others to 0
    curated_pi_samp = torch.where(coupling_matrix >= thresholds, coupling_matrix, torch.zeros_like(coupling_matrix))
    unfiltered_atac = curated_pi_samp @ original_atac.todense() / curated_pi_samp.sum(1).reshape(-1, 1)
    threshold = threshold

    mask = unfiltered_atac > threshold
    new_atac = np.where(mask, 1, 0)

    return new_atac


def keep_top_n_per_row(tensor, n):
    # Sort each row in descending order
    sorted_tensor, _ = torch.sort(tensor, descending=True, dim=1)
    # Get the threshold (n-th largest value) for each row
    thresholds = sorted_tensor[:, n - 1].unsqueeze(1)
    # Keep values >= threshold, set others to 0
    return torch.where(tensor >= thresholds, tensor, torch.zeros_like(tensor))

def transfer_accuracy(domain1, domain2, type1, type2,n):
	"""
	Metric from UnionCom: "Label Transfer Accuracy"
	"""
	knn = KNeighborsClassifier(n_neighbors=n)
	knn.fit(domain2, type2)
	type1_predict = knn.predict(domain1)
	count = 0
	for label1, label2 in zip(type1_predict, type1):
		if label1 == label2:
			count += 1
	return count / len(type1)


def preprocess_multimodal(rna_adata,atac_adata,rna_normalize_total=1e5,atac_normalize_total=1e4,
			   rna_filter_gene_percent=0.01,atac_filter_peak_percent=0.001):
	
	rna_adata.raw = rna_adata
	atac_adata.raw = atac_adata
	
	# get number of counts per cell
	rna_adata.obs['n_counts'] = np.array(rna_adata.X.sum(1)).squeeze()
	atac_adata.obs['n_counts'] = np.array(atac_adata.X.sum(1)).squeeze()

	# normalize counts & log transformation
	sc.pp.normalize_total(rna_adata,target_sum=rna_normalize_total)
	sc.pp.normalize_total(atac_adata,target_sum=atac_normalize_total)
	sc.pp.log1p(rna_adata)
	sc.pp.log1p(atac_adata)
	
	# filter genes & peaks
	n_cells = rna_adata.shape[0]
	sc.pp.filter_genes(rna_adata,min_cells=int(rna_filter_gene_percent*n_cells))
	sc.pp.filter_genes(atac_adata,min_cells=int(atac_filter_peak_percent*n_cells))

def set_dpt_root(rna_adata,start_marker_gene):
	
	marker_gene_values = rna_adata[:,start_marker_gene].X.toarray().squeeze()
	rna_adata.uns['iroot'] = np.argsort(-marker_gene_values)[0]

def schema_representations(rna_adata,atac_adata,schema_reference='rna',
						   num_hv_genes=2000,num_components=50):		
		
	rna_adata_copy = rna_adata.copy()
	atac_adata_copy = atac_adata.copy()    
	
	# highly variable genes
	sc.pp.highly_variable_genes(rna_adata_copy,n_top_genes=num_hv_genes,subset=True)
	sc.pp.scale(rna_adata_copy)    
	
	# TF-IDF transformation (use raw peak counts)
	tfidf = TfidfTransformer()
	idf = tfidf.fit_transform(atac_adata.raw[:,atac_adata.var.index.values].X)
	atac_adata_copy.X = idf
	
	# apply PCA
	sc.tl.pca(rna_adata_copy,n_comps=num_components)
	sc.tl.pca(atac_adata_copy,n_comps=num_components)
	rna_adata.obsm['X_pca'] = rna_adata_copy.obsm['X_pca']
	atac_adata.obsm['X_pca'] = atac_adata_copy.obsm['X_pca']

	# remove PCs that correlate strongly with total read counts
	atac_pca_inds2keep = []
	rna_pca_inds2keep = []
	for i in range(num_components):
		rho,p = spearmanr(atac_adata_copy.obsm['X_pca'][:,i],
						  atac_adata_copy.obs['n_counts'])
		if rho < 0.9:
			atac_pca_inds2keep.append(i)

		rho,p = spearmanr(rna_adata_copy.obsm['X_pca'][:,i],
						  rna_adata_copy.obs['n_counts'])
		if rho < 0.9:
			rna_pca_inds2keep.append(i)
			
	params = {"decomposition_model": "pca","num_top_components": num_components}
	model = SchemaQP(min_desired_corr=0.99,mode = 'scale',params=params)

	if schema_reference == 'rna':
		secondary_modalities = [atac_adata.obsm['X_pca'][:,atac_pca_inds2keep]]
		W = model.fit_transform(rna_adata.obsm['X_pca'][:,rna_pca_inds2keep], 
								secondary_modalities,['feature_vector'],[1])
	elif schema_reference == 'atac':
		secondary_modalities = [rna_adata.obsm['X_pca'][:,rna_pca_inds2keep]]
		W = model.fit_transform(atac_adata.obsm['X_pca'][:,atac_pca_inds2keep],
								secondary_modalities,['feature_vector'],[1])

	return W
	
def identify_all_peak_gene_link_candidates(rna_adata,atac_adata,distance_thresh=1e6,
										   rna_chr_key='chr_no',rna_txstart_key='txstart',
										   atac_chr_key='chr_no',atac_start_key='start',
										   atac_end_key='end'):
	
	gene_peak_dict = {}
	gene_peak_dist_dict = {}

	for chr_no in sorted(list(set(atac_adata.var[atac_chr_key].values))):

		filtered_atac_df = atac_adata.var[atac_adata.var[atac_chr_key] == chr_no]
		filtered_gene_df = rna_adata.var[rna_adata.var[rna_chr_key] == chr_no]
		filtered_gene_symbols = filtered_gene_df.index.values

		peak_data = filtered_atac_df[[atac_start_key,atac_end_key]].values
		peak_start = peak_data[:,0]
		peak_end = peak_data[:,1]

		gene_txstart = filtered_gene_df['txstart'].values

		distance_mat = peak_start[:,None].astype(int)-gene_txstart.astype(int)
		for atac_ind,gene_ind in zip(*np.where(abs(distance_mat) <= distance_thresh)):
			gene = filtered_gene_symbols[gene_ind]
			atac_idx = filtered_atac_df.index.values[atac_ind]
			if gene not in gene_peak_dict:
				gene_peak_dict[gene] = []
			gene_peak_dict[gene].append(atac_idx)
			gene_peak_dist_dict[(chr_no,gene,atac_idx)] = distance_mat[atac_ind,gene_ind]

		distance_mat = peak_end[:,None].astype(int)-gene_txstart.astype(int)
		for atac_ind,gene_ind in zip(*np.where(abs(distance_mat) <= distance_thresh)):
			gene = filtered_gene_symbols[gene_ind]
			atac_idx = filtered_atac_df.index.values[atac_ind]
			if gene not in gene_peak_dict:
				gene_peak_dict[gene] = []
			gene_peak_dict[gene].append(atac_idx)

			if (chr_no,gene,atac_idx) in gene_peak_dist_dict:
				if abs(distance_mat[atac_ind,gene_ind]) < abs(gene_peak_dist_dict[(chr_no,gene,atac_idx)]):
					gene_peak_dist_dict[(chr_no,gene,atac_idx)] = distance_mat[atac_ind,gene_ind]

	data_dict = {k: [] for k in ['chr_no','gene','atac_id','dist']}

	for (chr_no,gene,atac_idx),dist in gene_peak_dist_dict.items():
		data_dict['chr_no'].append(chr_no)
		data_dict['gene'].append(gene)
		data_dict['atac_id'].append(atac_idx)
		data_dict['dist'].append(int(dist))

	data_df = pd.DataFrame(data_dict)
	
	return data_df
	
def construct_dag(joint_feature_embeddings,iroot,n_neighbors=15,pseudotime_algo='dpt'):
	
	"""Constructs the adjacency matrix for a DAG.
	Parameters
	----------
	joint_feature_embeddings: 'numpy.ndarray' (default: None)
		Matrix of low dimensional embeddings with rows corresponding
		to observations and columns corresponding to feature embeddings
		for constructing a DAG if a custom DAG is not provided.
	iroot: 'int' (default: None)
		Index of root cell for inferring pseudotime for constructing a DAG 
		if a custom DAG is not provided.
	n_neighbors: 'int' (default: 15)
		Number of nearest neighbors to use in constructing a k-nearest
		neighbor graph for constructing a DAG if a custom DAG is not provided.
	pseudotime_algo: {'dpt','palantir'} 
		Pseudotime algorithm to use for constructing a DAG if a custom DAG 
		is not provided. 'dpt' and 'palantir' perform the diffusion pseudotime
		(Haghverdi et al., 2016) and Palantir (Setty et al., 2019) algorithms, 
		respectively.
	"""

	pseudotime,knn_graph = infer_knngraph_pseudotime(joint_feature_embeddings,iroot,
		n_neighbors=n_neighbors,pseudotime_algo=pseudotime_algo)
	dag_adjacency_matrix = dag_orient_edges(knn_graph,pseudotime)

	return dag_adjacency_matrix

def infer_knngraph_pseudotime(joint_feature_embeddings,iroot,n_neighbors=15,pseudotime_algo='dpt'):

	adata = AnnData(joint_feature_embeddings)
	adata.obsm['X_joint'] = joint_feature_embeddings
	adata.uns['iroot'] = iroot

	if pseudotime_algo == 'dpt':
		sc.pp.neighbors(adata,use_rep='X_joint',n_neighbors=n_neighbors)
		sc.tl.dpt(adata)
		adata.obs['pseudotime'] = adata.obs['dpt_pseudotime'].values
		knn_graph = adata.obsp['distances'].astype(bool).astype(float)
	elif pseudotime_algo == 'palantir':
		sc.pp.neighbors(adata,use_rep='X_joint',n_neighbors=n_neighbors)
		sce.tl.palantir(adata, knn=n_neighbors,use_adjacency_matrix=True,
			distances_key='distances')
		pr_res = sce.tl.palantir_results(adata,
			early_cell=adata.obs.index.values[adata.uns['iroot']],
			ms_data='X_palantir_multiscale')
		adata.obs['pseudotime'] = pr_res.pseudotime
		knn_graph = adata.obsp['distances'].astype(bool).astype(float)

	return adata.obs['pseudotime'].values,knn_graph

def dag_orient_edges(adjacency_matrix,pseudotime):

	A = adjacency_matrix.astype(bool).astype(float)
	D = -1*np.sign(pseudotime[:,None] - pseudotime).T
	D = (D == 1).astype(float)
	D = (A.toarray()*D).astype(bool).astype(float)

	return D

def construct_S0_S1(D):

	D_0 = D.copy()
	D_1 = D.copy() + np.eye(D.shape[0])
	S_0 = D_0.copy()
	S_1 = D_1.copy()

	D_0_sum = D_0.sum(1)
	D_0_sum[D_0_sum == 0] = 1
	S_0 = (S_0.T/D_0_sum)
	S_1 = (S_1.T/D_1.sum(1))

	return S_0,S_1

def load_multiome_data(data_dir,dataset,sampling=None,preprocess=True):

	if sampling == 'geosketch':
		atac_adata = sc.read(os.path.join(data_dir,
			'{}.atac.sketch.h5ad'.format(dataset)))
		rna_adata = sc.read(os.path.join(data_dir,
			'{}.rna.sketch.h5ad'.format(dataset)))
	elif sampling == 'uniform':
		atac_adata = sc.read(os.path.join(data_dir,
			'{}.atac.uniform.h5ad'.format(dataset)))
		rna_adata = sc.read(os.path.join(data_dir,
			'{}.rna.uniform.h5ad'.format(dataset)))
	else:
		atac_adata = sc.read(os.path.join(data_dir,
			'{}.atac.h5ad'.format(dataset)))
		rna_adata = sc.read(os.path.join(data_dir,
			'{}.rna.h5ad'.format(dataset)))

	if preprocess:

		# scale by maximum 
		# (rna already normalized by library size + log-transformed)
		X_max = rna_adata.X.max(0).toarray().squeeze()
		X_max[X_max == 0] = 1
		rna_adata.X = csr_matrix(rna_adata.X / X_max)

		# atac: normalize library size + log transformation
		sc.pp.normalize_total(atac_adata,target_sum=1e4)
		sc.pp.log1p(atac_adata)

	return rna_adata,atac_adata

def retain_desired_indices(rna_adata,atac_adata,rna_idx,atac_idx):

	atac_X = atac_adata.X[:,sorted(list(set(atac_idx)))].toarray()
	rna_X = rna_adata.X[:,sorted(list(set(rna_idx)))].toarray()

	sorted_atac_idx = sorted(list(set(atac_idx)))
	sorted_rna_idx = sorted(list(set(rna_idx)))

	atac_idx_map = {idx1:idx2 for idx1,idx2 in \
		zip(*[sorted(list(set(atac_idx))),range(atac_X.shape[1])])}
	rna_idx_map = {idx1:idx2 for idx1,idx2 in \
		zip(*[sorted(list(set(rna_idx))),range(rna_X.shape[1])])}

	atac_idx = np.array([atac_idx_map[idx] for idx in atac_idx])
	rna_idx = np.array([rna_idx_map[idx] for idx in rna_idx])

	return rna_X,atac_X,rna_idx,atac_idx

def rss_ratio_ftest(rss_ratio,n_layers,n):
	
	from scipy.stats import f

	# number of reduced model parameters 
	# (W,b for each RNA layer)
	p1 = n_layers * 2 
	# number of full model parameters 
	# (W,b for each layer of RNA, ATAC + interaction)
	p2 = n_layers * 2 * 2 + 1 
	
	dfn = 2 * n_layers + 1
	dfd = n - 4 * n_layers - 1

	F = (rss_ratio - 1) * (n - p2)/(p2-p1)

	return 1-f.cdf(F,dfn,dfd)

def calculate_fdr(p):
	return fdrcorrection(p)[1]
