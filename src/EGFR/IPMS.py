import pandas as pd
import numpy as np
import pylab as plt
from sklearn.preprocessing import normalize
from IPython import embed
import sys
from utils import cluster,heatmap

def determineInteractors(mat,entities,threshold=0.05):
    '''Determines interactors by looking at a threshold on the normalized gradient matrix'''
    gtgenes=[]
    for i in range(0,len(entities)):
        gtval = np.max(np.abs(mat), axis=0)[i] # Maximum gradient for gene
        if gtval>threshold:
            gtgenes.append(entities[i])
    return gtgenes
        
def extract_experiments(data,datapath,resultspath,metric='area',method='raw'):#,control='starvation 16h, no EGF'):
    '''Extracts the Mass Spec data, when multiple experiments are stored in a single file.'''
    ############################################
    NORMMETHOD='l1'  # 'l2' destroys signal from low throughput
    ############################################

    # Get data and group by experiment

    df=pd.read_csv('{}/IPMSData.txt'.format(datapath),delimiter='\t')
    #print min(df['_e2g_nGPArea_Sum_cgpAdj']), np.median(df['_e2g_nGPArea_Sum_cgpAdj']),np.mean(df['_e2g_nGPArea_Sum_cgpAdj'])

    ############################################################################################
    min_value=min(df[df['iBAQ']!=0]['iBAQ'])# "min value"
    print(min_value) # 3.02044e-05
    #embed()
    df['iBAQ']=df['iBAQ'].fillna(min_value)

    ########################################################################################
    print(len(df['iBAQ']), "original length")
    #df=df[df['_e2g_nGPArea_Sum_cgpAdj'] >0.03]
    ########################################################################################

    df=df[df['GeneName']!='EGFR']

    #gb=df.groupby('Experiment')

    lines=open('{}/ExperimentsKey.csv'.format(datapath)).read().replace('\r','\n').split('\n')
    expmapping={}
    useexp={}
    experimentset={}
    for line in lines[1:]:
        if line=='':
            continue
        line=line.strip().split(',')
        exp=int(line[0])
        expmapping[exp]=float(line[1])
        useexp[exp]=int(line[3])
        experimentset[exp]=int(line[2])
    df['Set']=[experimentset[x] for x in df['Experiment']]
    df['Use']=[useexp[x] for x in df['Experiment']]
    # Get experiments and all genes in IPMSs
    experiments=sorted(list(set(df['Experiment'])))
    experiments=[ x for x in experiments if useexp[x]]
    allgenes=sorted(list(set(df['GeneName'])))
    print('total genes',len(allgenes))

    #Get unique experiments
    uniqueexperiments=sorted(list(set([expmapping[x] for x in experiments])))
    #replicates=len(experiments)/len(uniqueexperiments)
     
    # Initialize a matrix for each experiment set to store the IPMS data 
    # for each experiment on every gene
    mats=[]
    raw_mats=[]
    for expset in set(experimentset.values()):
        data=np.empty((len(uniqueexperiments),len(allgenes)))
        data.fill(min_value)
        
        rawdata=np.empty((len(uniqueexperiments),len(allgenes)))
        rawdata.fill(min_value)
        labels=[str(x) for x in uniqueexperiments]
        #print labels
        times=uniqueexperiments
        #times=times*60.
                
        toaverage=[]
        toaverageind=[]
        tmpdf=df[df['Set']==expset]
        tmpdf=tmpdf[tmpdf['Use']==1]
        gb=tmpdf.groupby('Experiment')
        tmpexperiments=list(set(tmpdf['Experiment']))

        # Fill matrix
        for i in range(0,len(tmpexperiments)):
            experiment=tmpexperiments[i]
            dfa=gb.get_group(experiment)
    
            genes=list(dfa['GeneName'])
            #print genes
            ind=[allgenes.index(x) for x in genes]
            pair=(labels.index(str(expmapping[experiment])),ind)
            
            if metric=='area':
                if pair in toaverage:
                    data[pair]=data[pair]+dfa['iBAQ']
                    print('HERE')
                    exit()
                else:
                    data[pair]=dfa['iBAQ']
                    rawdata[pair]=dfa['iBAQ']
                    toaverage.append(pair)
                    for p in pair[1]:
                        toaverageind.append((pair[0],p))
            else:
                print('Please choose "area"')
                sys.exit()
            toaverage=list(set(toaverageind))
            
        # Normalization
        data[data==0]=min_value
        rawdata[rawdata==0]=min_value
        if method=='norm':
            data=normalize(data,norm=NORMMETHOD,axis=0) # Normalize each gene
            #data=np.nan_to_num(data/data[labels.index(control),:])
        elif 'gradient' in method:
            if 'norm' in method:
                data=normalize(data,norm=NORMMETHOD,axis=1)
                data=normalize(data,norm=NORMMETHOD,axis=0) # Normalize each gene
            #embed()
            #data=-np.log10(data)
            new=np.zeros_like(data)
            for i in range(1,len(uniqueexperiments)):
                new[i,:]=(data[i,:]-data[i-1,:])/(times[i]-times[i-1])
                #data[i,:]=(data[i,:]-data[0,:])/times[i]
            data=new
            data=data[1:,:]
        else:
            print('defaulting to no normalization')            
        mats.append(data)
        raw_mats.append(rawdata)
    #mbed()
    raw=np.mean(np.array(raw_mats),axis=0)
    data=np.mean(np.array(mats),axis=0)
    if 'gradient' in method:
        labels=labels[1:]
        times=times[1:]

    # Write matrix to file
    f=open('{}/Matrix_{}_{}.txt'.format(resultspath,metric,method),'w')
    w,h=data.shape  # @UnusedVariable
    f.write('Genes\t{}\n'.format('\t'.join(labels)))
    for i in range(0,h):
        s='\t'.join([str(x) for x in list(data[:,i])])
        f.write("{}\t{}\n".format(allgenes[i],s))
    f.close()
    
    
    # Visualize changes as heatmap
    heatmap(data,labels,allgenes,'noclusters_{}_{}'.format(metric,method),resultspath)
    #embed()
    # Plot all gene trends together
    plt.figure()
    plt.plot(times,np.mean(data, axis=1),'--',lw=3)
    #plt.xscale('log')
    for i in range(0,len(allgenes)):
        plt.plot(times,data[:,i],alpha=0.01)
    plt.savefig("{}/allchanges_{}_{}.pdf".format(resultspath,metric,method))
    plt.xscale('linear')
    plt.close()

    # Cluster
    clusters,dev=cluster(data,labels,allgenes,times,resultspath,metric+'_'+method)

    return data,allgenes,labels,(clusters,dev),raw
