# In[1]:
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn import metrics
import plotly.offline as py
py.init_notebook_mode(connected=False)

import plotly.graph_objs as go
from plotly.graph_objs import XAxis, YAxis, ZAxis, Scene
from sklearn.decomposition import FastICA as ICA
from sklearn.manifold import LocallyLinearEmbedding as LLE
from sklearn.manifold import SpectralEmbedding as LaplacianEigenMaps
from sklearn.manifold import Isomap
# In[2]:
 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.verbosity = 3            
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80)

# In[3]:
adata=sc.read_text("C:/Users/saite/Desktop/Datasets/dataset1.txt")
#adata=sc.read_csv("C:/Users/saite/Desktop/Datasets/wang.csv")
# =============================================================================
#reading PBMC dataset
# adata=sc.read_10x_mtx(
#     'C:/Users/saite/Desktop/Datasets/PBMC/filtered_gene_bc_matrices/hg19',  # the directory with the `.mtx` file
#     var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
#     cache=True)
# =============================================================================
#adata=adata.transpose()
adata.var_names_make_unique()
#adata.obs_names_make_unique()
# In[4]:
print(adata)
sc.pl.highest_expr_genes(adata, n_top=20)
#Computes, for each gene, the fraction of counts assigned to that gene within a cell. 
#top n genes with the highest mean fraction over all cells
#usual suspects”, i.e., mitochondrial genes, actin, ribosomal protein, MALAT1

sc.pp.filter_cells(adata, min_genes=200)#Filter cell outliers based on counts and numbers of genes expressed. 
#cells with less than 200 expressed genes (UMI count>0)
sc.pp.filter_genes(adata, min_cells=3)#Filter genes based on number of cells or counts. genes expressed in less than 3 cells are removed
print(adata)

# In[5]:
#filtering low quality genes
mito_genes = adata.var_names.str.startswith('MT-')
# for each cell compute fraction of counts in mito genes vs. all genes
adata.obs['percent_mito'] = np.ravel(np.sum(
    adata[:, mito_genes].X, axis=1)) / np.ravel(np.sum(adata.X, axis=1))
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = np.ravel(adata.X.sum(axis=1))

sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],jitter=0.4, multi_panel=True)
sc.pl.scatter(adata, x='n_counts', y='percent_mito')
sc.pl.scatter(adata, x='n_counts', y='n_genes')

# In[6]:
#selecting number of genes and MT genes for downstream analysis
adata = adata[adata.obs.n_genes < 2500, :]
adata = adata[adata.obs.percent_mito < 0.05, :]
print(adata)

# In[7]:
#Normalization and log-scaling of data
sc.pp.normalize_total(adata, target_sum=1e4)#this is CPM normalization
#f exclude_highly_expressed=True, very highly expressed genes are excluded from the computation of the normalization factor (size factor) for each cell. 
#This is meaningful as these can strongly influence the resulting normalized values for all other genes
sc.pp.log1p(adata)#Logarithmize the data.
sc.pl.highest_expr_genes(adata, n_top=20)
print(adata)

# In[8]:
#extracting highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
#The normalized dispersion is obtained by scaling with the mean and standard deviation of the dispersions for genes falling into a given bin for mean expression of genes. 
#This means that for each bin of mean expression, highly variable genes are selected.

# In[9]:
print(adata)
sc.pl.highly_variable_genes(adata)
adata = adata[:, adata.var.highly_variable]
print(adata)

# In[10]:
#regressing out technical and biological noise
sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])#Regress out (mostly) unwanted sources of variation.
sc.pp.scale(adata, max_value=10)

# In[11]:
#pre-processed data
print(adata)
ppd=adata.X

# In[12]:
def Kmeans_2D(XT,plt_title,k):
    kmeans = KMeans(n_clusters=k).fit(XT)
# Compute cluster centers and predict cluster indices
    X_clustered = kmeans.labels_

    trace_Kmeans = go.Scatter(x=XT[:, 0], y= XT[:, 1],mode="markers",
                    showlegend=False,
                    text = adata.var_names,
                    marker=dict(
                            size=8,
                            color = X_clustered,
                            colorscale = 'Portland',
                            showscale=False, 
                            line = dict(
            width = 2,
            color = 'rgb(255, 255, 255)'
        )))

    layout = go.Layout(
    title= plt_title,
    plot_bgcolor='rgba(0,0,0,0)',
    paper_bgcolor='rgba(0,0,0,0)',
    hovermode= 'closest',
     xaxis= {'showgrid': False,
             'showline': False,
             'zeroline':False,
             'color':"rgb(255,255,255)"
             },
    yaxis={'showgrid': False,
           'showline': False,
           'zeroline':False,
           'color':"rgb(255,255,255)"
             },
    showlegend= False)

    data = [trace_Kmeans]
    fig1 = dict(data=data, layout= layout)
    py.plot(fig1, filename=plt_title)
    print('Clustering Score',silhouette_score(XT,X_clustered))
    return X_clustered
# In[13]:
def Kmeans_3D(XT,plt_title,k):
    kmeans = KMeans(n_clusters=k)
# Compute cluster centers and predict cluster indices
    X_clustered = kmeans.fit_predict(XT)

    fig = go.Figure(data=[go.Scatter3d(x=XT[:, 0], y= XT[:, 1], z=XT[:, 2], mode="markers",
                    showlegend=False,
                    text = adata.var_names,
                    marker=dict(
                            size=8,
                            color = X_clustered,
                            colorscale = 'Portland',
                            showscale=False, 
                            line = dict(
            width = 2,
            color = 'rgb(255, 255, 255)'
        )))])

    layout = fig.update_layout(
    title= plt_title,
    hovermode= 'closest',
    scene = Scene(bgcolor='rgba(255,255,255,255)',
        xaxis={'showgrid': False,
               'backgroundcolor':'rgba(255,255,255,255)',
           'showline': False,
           'zeroline':False,
           'color':"rgb(255,255,255)"
             },
       yaxis={'showgrid': False,
           'showline': False,
           'zeroline':False,
           'backgroundcolor':'rgba(255,255,255,255)',
           'color':"rgb(255,255,255)"
             },
       zaxis={'showgrid': False,
           'showline': False,
           'zeroline':False,
           'backgroundcolor':'rgba(255,255,255,255)',
           'color':"rgb(255,255,255)"
             }),
    showlegend= False)

    data = fig
    fig1 = dict(data=data, layout= layout)
    py.plot(fig1, filename=plt_title)
    print('Clustering Score',silhouette_score(XT,X_clustered))  
# In[14]:
#ICA on reduced MLLE
def ICA_2D(XT,plt_title,n):
    iso=ICA(n_components=2)#isomap
    XT=iso.fit_transform(XT)
    labels=Kmeans_2D(XT,plt_title,n)
    return XT,labels
# In[15]:
#function to get optimal parameters
def nn_check(ppd):
    for i in range(8,26):
        lle=LLE(n_components=3, n_neighbors=i,method='modified',modified_tol=1e-12)
        XT=lle.fit_transform(ppd)
        print('running')
        validity(XT,i)
    print('done')        
        
# In[16]:
#Validity of indices to evaluate clustering performance
def validity(val_data,i):
    for n in range(4,15):
        km=KMeans(n_clusters=n)
        clustered=km.fit_predict(val_data)#silhouette_score,  calinski_harabasz_score
        sh=metrics.silhouette_score(val_data, clustered)
        ch=metrics.calinski_harabasz_score(val_data, clustered)
        db=metrics.davies_bouldin_score(val_data, clustered)
        df=pd.read_csv('C:/Users/saite/Downloads/validity1.csv')
        df2 = df.append(pd.DataFrame([[i,n,sh,ch,db]], columns=df.columns))
        df2.to_csv (r'C:/Users/saite/Downloads/validity1.csv',index=False)
# In[17]:
#Standard LLE
def SLLE():
    le = LLE(n_components=3,n_neighbors=14)
    slle=le.fit_transform(ppd)
    km_slle2 = Kmeans_2D(slle,"KM Clustering on 2D Standard LLE.html",8)
    km_slle3 = Kmeans_3D(slle,"KM Clustering on 3D Standard LLE.html",8)

# In[18]:
def Iso_map():
    iso=Isomap(n_components=3,n_neighbors=9)
    isomaps=iso.fit_transform(ppd)
    km_iso2 = Kmeans_2D(isomaps,"KM Clustering on 2D Isomaps.html",8)
    km_iso3 = Kmeans_3D(isomaps,"KM Clustering on 3D Isomaps.html",8)

# In[19]:
def Laplacian_eigenmap():
    le=LaplacianEigenMaps(n_components=3, n_neighbors=10)
    lmaps=le.fit_transform(ppd)
    km_lmaps2 = Kmeans_2D(lmaps,"KM Clustering on 2D Laplacian Eigenmaps.html",8)
    km_lmaps3 = Kmeans_3D(lmaps,"KM Clustering on 3D Laplacian Eigenmaps",8)

# In[20]:
SLLE()
Iso_map()
Laplacian_eigenmap()

# In[21]:
#to get the best parameters for dimensionality reduction and clustering
nn_check(ppd)
# In[22]:
#Modified LLE 
lle=LLE(n_components=5, n_neighbors=8,method='modified',modified_tol=1e-12)
middle=lle.fit_transform(ppd)#passing adata.X giving different clustering results

# In[23]:
lle=LLE(n_components=3, n_neighbors=11,method='modified',modified_tol=1e-12)
reduced_lle=lle.fit_transform(ppd)
# In[24]:
km_mds2 = Kmeans_2D(reduced_lle,"KM Clustering on 2D Modified LLE.html",7)
km_mds3 = Kmeans_3D(reduced_lle,"KM Clustering on 3D Modified LLE.html",7)

# In[25]:
#call ICA function
reduced_ica, ica_km = ICA_2D(reduced_lle,"ICA+KM on MLLE.html",7)

# In[26]:
#saving data and clustering labels 
df=pd.DataFrame(reduced_ica)
df.to_csv (r'C:/Users/saite/source/mlle.csv', index = False, header=True)

df=pd.DataFrame(ica_km)
df.columns=['clusters']
df.to_csv (r'C:/Users/saite/source/km.csv', index = False, header=True)

# In[27]:
#reading back to insert into scanpy
mlle=pd.read_csv("C:/Users/saite/source/mlle.csv")
km=pd.read_csv("C:/Users/saite/source/km.csv")

adata.obsm['X_pca']=np.asanyarray(mlle)
adata.obs['km']=list(km['clusters'])
adata.obs['km']=adata.obs['km'].astype('category')
sc.pl.pca(adata,color=['km'])
# In[28]:
#ranking genes from the clusters
sc.tl.rank_genes_groups(adata, 'km', method='wilcoxon')#wilcoxon, t-test, logreg
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)

# In[29]:
#storing top marker genes into a variable    
marker_genes=pd.DataFrame(adata.uns['rank_genes_groups']['names'])
#mg1=marker_genes2.stack().reset_index()
z=marker_genes[:25]
# In[30]:
#renaming the clusters with cell types
new_cluster_names = [
    'SARS-Cov', 'Neutrophil2',
    'Club cell',
    'Th17/iTreg-stimulated CD4+ central memory T cell','Natural killer cell',
    'Stem/Club/Hillock epithelial cell','Club cell lung']
adata.rename_categories('km', new_cluster_names)
sc.pl.pca(adata, color='km', title='Cluster Annotation', frameon=False, save='.pdf')

# In[31]:
#visualizing marker genes 
sc.pl.rank_genes_groups_stacked_violin(adata,n_genes=5)
sc.pl.rank_genes_groups_dotplot(adata,n_genes=5)

# In[32]:
#visualizing top marker genes with its clusters
marker_genes_dict={'SARS-Cov': ['SCoV1-ORF7a','SCoV1-ORF6','TOP2A','GINS2','CKB'],
'Neutrophil': ['ISG15', 'RSAD2', 'IFIT2', 'IFIT3', 'IFIT1', 'XAF1', 'OASL'],
'Club-cell': ['GLUL','SLPI', 'PSCA', 'NEAT1'],
'Th17/iTreg-stimulated_CD4+central_memory_T_cell': ['IFIT3', 'IFIT1', 'OASL', 'IFIT2', 'RSAD2', 'RNF213', 'ISG15', 'MX1'],
'natural_killer_cell': ['IFIT3','IFIT1', 'MX1', 'SAMD9', 'RSAD2', 'RNF213', 'SAMD9L', 'XAF1', 'IFI6', 'OAS1'],
'Stem/Club/Hillock_epithelial_cell': ['PSCA','AQP3','DDIT4','GLUL'],
'Club_cell_Lung': ['MSMB','S100A4','PIGR','SAT1']}

sc.pl.matrixplot(adata, marker_genes_dict, 'km', dendrogram=True,cmap='Blues', standard_scale='var')
sc.pl.dotplot(adata, marker_genes_dict, 'km', dendrogram=True)
sc.pl.stacked_violin(adata, marker_genes_dict, groupby='km', swap_axes=False, dendrogram=True)
sc.pl.heatmap(adata, marker_genes_dict, groupby='km', cmap='viridis', dendrogram=True)

# In[33]:
#Cluster annotation
sc.pl.pca(adata, color='km', add_outline=True, 
           legend_fontsize=12, legend_fontoutline=2,frameon=False,
           title='cluster annotation', palette='Set1')
#clusters with outline
sc.pl.pca(adata, color='km', add_outline=True, 
           legend_fontsize=12, legend_fontoutline=2,frameon=False,
           title='cluster annotation', palette='Set1')