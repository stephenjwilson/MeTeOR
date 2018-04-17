from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans
import numpy as np
import seaborn as sns
import pylab as plt
import matplotlib
from itertools import cycle
from IPython import embed
import copy,itertools
import pyupset as pyu
import pickle
import pandas as pd
import pylab as pl

### Process Data Files #####################################################
def process_network(netfl,entity,threshold):
    '''Extracts all relationships related to an entity from a network with header demarked by #.'''
    genes=[]
    w=[]
    for line in open(netfl).readlines():
        line=line.strip().split('\t')
        if line[0]=='#':
            continue
        if entity in line:
            #print(line)
            line.pop(line.index(entity))# pop the entity of interest out of line
            if line[0] not in genes:
                tmp=float(line[1])
                if tmp>threshold:
                    genes.append(line[0])
                    w.append(tmp)
    return [genes,w]

def getNetGS(netfl,entity,mapping,threshold=0):
    '''Processes a network into a gold standard'''
    genes,_ =process_network(netfl, entity,threshold)
    # Map
    genes=[mapping[gene] for gene in genes if gene in mapping.keys()]
    
    return genes
    
def gen_randlist(size,netfl,resultspath,samples=10):
    '''Generate a Random'''
    genes=[]
    for line in open(netfl).readlines():
        line=line.strip().split('\t')
        if line[0]=='#':
            continue
        for gene in line[:2]:
            if 'G.' in gene:
                genes.append(gene)
    ls=[]
    
    for i in range(samples):  # @UnusedVariable
        randls=np.random.choice(genes,size)
        weights=np.arange(0,len(randls))
        ls.append((randls,weights))
    return ls

def getPathwayGenes(GSEA='../Data/GSEA.tsv',thresh=2./3):
    '''Uses GSEA to get EGFR associated by pathways.'''
    lines=open(GSEA).readlines()
    dic={}
    for line in lines:
        line=line.strip().split('\t')
        conf=1./len(line) # The probability of gene associations is related the the inverse of the number of genes
        for pair in itertools.combinations(line[1:],2):
            try:
                dic[pair]+=(conf) # The confidence is accumulated over the pathways 
            except:
                dic[pair]=(conf)
            
            pair=(pair[1],pair[0]) # ensure symmetrical
            try:
                dic[pair]+=(conf)
            except:
                dic[pair]=(conf)
    items=list(zip(*list(dic.items())))
    thresh=np.percentile(items[1],thresh*100.)

    # Threshold data
    gtgenes=[]

    for pair in dic:
        if 'EGFR' in pair:
            if dic[pair]>thresh:
                pair=list(pair)
                pair.pop(pair.index('EGFR'))
                gtgenes.append(pair[0])
    
    return list(set(gtgenes))
    
    
def getCoReg(fl='EGFR_cor_pancancer.tsv',datadir='../../data',thresh=0.25):
    '''Retrieves the list of genes coregulated by EGFR PanCancer'''
    lines=open('%s/%s'%(datadir,fl)).readlines()
    genes=[]
    for line in lines[1:]:
        line=line.strip().split('\t')
#         if np.abs(float(line[-1]))<thresh: # according to q
#             genes.append(line[0].strip())
        if np.abs(float(line[2]))>thresh: # according to r
            genes.append(line[0].strip())
    return genes

def getCancerGenes(fl='../Data/CancerGenes/CancerGenes.csv',sep=',',col=2):
    '''Get Genes associated with Cancer from COSMIC and MutSig'''
    lines=open(fl).readlines()
    genes=[]
    header=lines[0][col]
    for line in lines:
        genes.append(line.strip().split(sep)[col])
    return (header,genes)

def process_MeTeOR_OUTPUT(fl):
    '''Reads and processes the output from MeTeOR single entity flat files.'''
    lines=open(fl).readlines()
    per_row = []
    titles=lines[0].strip().replace('_',' ').split('\t')[1:]
    mgenes=[]
    predFilt=[[],[]]
    for line in lines[1:]:
        line=line.strip().split('\t')
        mgenes.append(line[0])
        if float(line[1])>2:
            predFilt[0].append(line[0])
            predFilt[1].append(float(line[6]))
        per_row.append([float(x) for x in line[1:]])
    weights=list(zip(*per_row))
    ls=[[list(copy.copy(mgenes)),list(w)] for w in weights]
    
    # add in filtered NMF
    ls.append(predFilt)
    titles.append('NMF_Filt')
    return ls,titles
############################################################################

### Plotting Aids ###
def heatmap(mat,experiments,genes,title,resultspath,clusters=[],genesort=[],control='starvation 16h, no EGF'):
    '''Produces Heatmaps for the matrix'''
    explabels=list(experiments)
    try:
        explabels[explabels.index(control)]='No EGF'
    except:
        pass

    if np.array(genesort).size!=0:
        mat = mat[:, genesort]
        genes=np.array(genes)[genesort]
    
    f,ax = plt.subplots(figsize=(11,9))
    #cmap=sns.diverging_palette(220,10,as_cmap=True)
    cmap=sns.color_palette("RdBu_r", 7)
    #embed()
    #sns.heatmap(mat,yticklabels=explabels,cmap=cmap,ax=ax)
    sns.heatmap(mat,yticklabels=explabels,ax=ax,cmap=cmap)
    ax.set_xticks([])
    f.savefig("{}/heatmaps/heatmap_{}.pdf".format(resultspath,title))

def plot_cluster_signals(mat,labels,n_clusters,plttitle,times,resultspath):
    '''Plots the clustering signals'''
    COLORS=['r', 'b', 'y','c','m','g','pink','navy','darkorange','indigo','lightsteelblue','olive']
    colorcycle=cycle([COLORS[i] for i in range(0,n_clusters)])
   
    # Plot clusters and mean signal
    mid=np.ceil(n_clusters/2.)
    matplotlib.rcParams.update({'font.size': 10})
    fig, axes = plt.subplots(nrows=2,ncols=int(mid),sharex=True,sharey=True)
    errors=[]
    n=0
    clusdev=[]

    for i in range(0,n_clusters):
        color=next(colorcycle)
        sub=mat[:,labels==i]
        w,h=sub.shape  # @UnusedVariable
        if i<mid:
            row=0
            col=int(i)
        else:
            row=1
            col=int(i-mid)
        mean=np.mean(sub, axis=1)
        clusdev.append(np.sum(np.abs(mean)))
        #print mean
        times=[float(time) for time in times]
        try:
            axes[row][col].plot(times,mean,'--r',lw=2)
        except:
            embed()
            exit()
        axes[row][col].set_title('{}'.format(h))

        sumerror=0
        for o in range(0,h): # for each gene sum the error and plot a scatter of its points
            sumerror+=sum((mean-sub[:,o])**2)
            n+=1
            axes[row][col].scatter(times,sub[:,o],alpha=0.2,c=color,edgecolor=color)
        errors.append(sumerror)
        axes[row][col].set_xscale('symlog')
        axes[row][col].set_xlim([1,150])

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    #plt.title("Error={}".format(sumerror/float(n)))
    plt.xlabel("Time (Minutes)")
    plt.ylabel("Protein Change ({})\n".format(plttitle))
    plt.tight_layout()
    plt.savefig("{}/genechanges/{}.pdf".format(resultspath,plttitle))
    plt.close()
    matplotlib.rcParams.update({'font.size': 20})   
    return sum(errors),clusdev

def upset(sets=[[]],resultdir='../../results/validation/EGFR'):
    '''
    sets= list of lists
    '''
    #with open('./test_data_dict.pckl', 'rb') as f:
    #   data_dict = pickle.load(f)
    
    #pyu.plot(data_dict)
    for fl,test in sets:#[['EGFRSetsp.csv','MeTeOR'],['EGFRSetsp2.csv','MeTeORSimple'],['EGFRSetsPredp.csv','MeTeORPred'],['EGFRSetsPredp2.csv','MeTeORPredSimple']]:
        df=pd.read_csv('{}/{}'.format(resultdir,fl),sep='\t')
        mydata={}
        title=fl
        colNames=list(df.columns[2:])
        sNone=[x for x in colNames if 'None' in x][0]
        for col in df.columns[2:]:
                #print(col)
                mydata[col]=pd.DataFrame(df['Genes'][df[col]==1])
        
        for arg in ['size','degree']:
            pyu.plot(mydata, colNames,title,resultdir,sort_by=arg)
            #pyu.plot(mydata, sort_by=arg, query=[('IPMS')])
            pl.savefig('{}/UpsetOverlap_None_{}_{}.pdf'.format(resultdir,arg,test))
    
        mydata.pop(sNone,None)
        colNames.pop(colNames.index(sNone))
        for arg in ['size','degree']:
            pyu.plot(mydata, colNames,title,resultdir,sort_by=arg)
            pl.savefig('{}/UpsetOverlap_{}_{}.pdf'.format(resultdir,arg,test))

### Algorithms ###
def cluster(mat,experiments,genes,times,resultspath,title=''):
    '''Cluster the matrix.'''
    bestlabels=[]
    besterror=''
    bestmethod=''
    r=[9]#range(3,9,1)
    print(r)
    for method in ['kmeans']:#,'agglomerative']:
        for n_clusters in r:
            if method=='agglomerative':
                for metric in ["euclidean", "cityblock"]:#["cosine", "euclidean", "cityblock","jaccard"]):
                    model = AgglomerativeClustering(n_clusters=n_clusters,
                                            linkage="average", affinity=metric)
                    model.fit(mat.T)
                    #print Counter(model.labels_)
                    indexes=np.argsort(model.labels_)
                    plttitle="{}_{}_{}_{}".format(title,method,metric,n_clusters)
                    heatmap(mat,experiments,genes,plttitle,resultspath,model.labels_,indexes)
                    error,clusdev=plot_cluster_signals(mat,model.labels_,n_clusters,plttitle,times,resultspath)
                    if besterror=='':
                        go=1
                    else:
                        if besterror>error:
                            go=1
                        else:
                            go=0
                    if go:
                        besterror=error
                        bestlabels=model.labels_
                        bestmethod=plttitle
                        bestdev=clusdev
            if method=='kmeans':
                model = KMeans(n_clusters=n_clusters, init='k-means++', n_init=100,
                               precompute_distances='auto',copy_x=True)
                model.fit(mat.T)
                indexes=np.argsort(model.labels_)
                #print Counter(model.labels_)
                plttitle="{}_{}_{}".format(title,method,n_clusters)
                heatmap(mat,experiments,genes,plttitle,resultspath,model.labels_,indexes)
                error,clusdev=plot_cluster_signals(mat,model.labels_,n_clusters,plttitle,times,resultspath)

                if besterror=='':
                    go=1
                else:
                    if besterror>error:
                        go=1
                    else:
                        go=0
                if go:
                    besterror=error
                    bestlabels=model.labels_
                    bestmethod=plttitle
                    bestdev=clusdev
    print(besterror,bestmethod)
    return bestlabels,bestdev
