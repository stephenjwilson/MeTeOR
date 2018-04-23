### Imports ###
import os
import matplotlib
import pickle
import numpy as np
from IPython import embed

from utils import  process_MeTeOR_OUTPUT, getNetGS, getCoReg, upset#,process_network, getPathwayGenes, getCoReg, getNetGS @UnresolvedImport
from IPMS import extract_experiments, determineInteractors 
from assess_ranking import assess_ranking,formatForUpset 


### Plotting Setting ###
font = {'weight': 'bold',
        'size': 20}
matplotlib.rc('font', **font)
matplotlib.rcParams.update({'figure.autolayout': True})


def run(datadir,resultdir,gene):
    
    outputdir='%s/validation/EGFR' % (resultdir)
    # Make Directories
    try:
        os.mkdir(outputdir)
    except:
        pass
    try:
        os.mkdir('{}/genechanges'.format(outputdir))
    except:
        pass
    try:
        os.mkdir('{}/heatmaps'.format(outputdir))
    except:
        pass
    try:
        os.mkdir("{}/AccuracyPlots".format(outputdir))
    except:
        pass
    
    
    ### Prepare data ###
    print('Prepare Data')
    # ls = zip(*[x[1] for x in networkdata])
    # titles = [x[0] for x in networkdata]
    ls=[[],[]] # Toggle for network input
    titles=[]
    
    # Get MeTeOR output
    nets, labels = process_MeTeOR_OUTPUT('%s/validation/mrr300_MeTeOR%s.txt' %(resultdir,'G'+gene)) # need to replace k automatically
    part1, part2 = zip(*nets)
    ls[0] += part1
    ls[1] += part2
    titles += labels
    ls = [list(x) for x in ls]
     
    ### Set-up IPMS ###
    MS = '%sIPMSData.txt' % (datadir)
    measure = 'area'
    method = 'gradientnorm'
     
    ### Begin Analysis ###
    data = []
    # Get IPMS data
    try:
        f = open("{}/{}_{}.pkl".format(outputdir, method, measure), 'rb')
        msdata, msgenes, experiments, clusters, raw = pickle.load(f)
        f.close()
    except:
        msdata, msgenes, experiments, clusters, raw = extract_experiments(
            MS, datadir, outputdir, metric=measure, method=method)
        # exit()
        #msgenes = mapGenels(msgenes, resultdir)
        f = open("{}/{}_{}.pkl".format(outputdir, method, measure), 'wb')
        pickle.dump((msdata, msgenes, experiments, clusters, raw), f)
        f.close()
         
     
    ### Get GTs ###
    print('Get GTs')
    ### Obtain Mapping ###
    mappingfl=open('{}/mapping/gene/HGNC_2-8-18.txt'.format(datadir))
    mapping={}
    
    for line in mappingfl:
        line=line.strip('\n').split('\t')
        mapping[line[5]]=line[1]
    
    # Compare to STRING EXP
    netfl='%s/networks/STRING10_experimental.txt' % (resultdir)
    threshold=0
    STRINGexp = getNetGS(netfl,gene,mapping,threshold)
    # Compare to CoReg genes
    cancergenes = getCoReg()
    #  Extract IPMS interactions
    IPMSgenes = determineInteractors(msdata, msgenes)
    # Get Pathway Data
    netfl='%s/networks/MSigDBCurated_top_cp.txt' % (resultdir)
    threshold=0
    pathwaygenes = getNetGS(netfl,gene,mapping,threshold)
    # Collect and then assess
    gts = [IPMSgenes, cancergenes, STRINGexp,pathwaygenes]
    gttitles = ['Prospective IPMS', 'Cancer Co-Regulation', 'STRING Experimental','Pathway Knowledge']
    chosengt=3

    ### Assess ranking ###
    print('Assess Rankings')
    outfls=assess_ranking(ls[1], ls[0], titles, gts, gttitles,
                   "{}_{}".format(measure, method), resultdir+'/validation/EGFR',chosengt)
    
    # Format for upset
    formatForUpset(outfls,resultdir,cutoff=20)
    
    # run upset
    upset(sets=[['data__MeTeOR Weight_area_gradientnormupset.txt','MeTeOR'],['data__MeTeOR NMF Pred_area_gradientnormupset.txt','MeTeOR NMF'],['data__NMF_Filt_area_gradientnormupset.txt','MeTeOR NMF Diff']],resultdir=resultdir+'/validation/EGFR')
#     # Get relative rankings and make plots
#     runPlots()
    
if __name__=='__main__':
    ### Set-up ###
    datadir = '../../data'
    gene = '1956'  # EGFR mapped
    resultdir = '../../results'
    run(datadir=datadir,resultdir=resultdir,gene=gene)

