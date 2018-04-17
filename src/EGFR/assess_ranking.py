import matplotlib
import numpy as np
import pylab as plt
from IPython import embed
from itertools import cycle
from operator import itemgetter
from sklearn.metrics import roc_curve, auc, precision_recall_curve
import pandas as pd

def formatForUpset(fls,resultdir,cutoff=100):
    ''' Creates files that are formatted for the the upset library'''   
    cutoff=cutoff-1
    for fl in fls:
        df=pd.read_csv(fl,sep='\t')
        header=list(df.columns.values)
        #header=header[2:]
        # Need to limit to the same genes as the orig file.
        df=df.loc[:cutoff,header]
        # Create None column
        df['None']=np.logical_not(df.loc[:,header[2:]].sum(axis=1)>0)*1
        # put none at the beginning
        inds=list(range(0,len(header)))
        inds.insert(2,len(header))
        df=df.iloc[:,inds]
        # Add Counts to columns
        #df.cols.values
        mapping=df.iloc[:,2:].sum(axis=0)
        mapping=dict((mapping.items()))
        for key in mapping:
            mapping[key]='{} ({:d})'.format(key,int(mapping[key]))
        df.rename(columns=mapping, inplace=True)
        # Put to file
        
        df.to_csv(fl.replace('.txt','upset.txt'),sep='\t', index=False)

def ismember(a, b):
    '''ismember function that returns indicies of a in b'''
    bind = {}
    for i, elt in enumerate(b):
        if elt not in bind:
            bind[elt] = i
    return [bind.get(itm, None) for itm in a]


def truthVector(rgenes, gt):
    '''Creates a truth vector given a the two lists of genes.'''
    truth = np.zeros([len(rgenes)])
    inds = ismember(gt, rgenes)
    inds = [x for x in inds if x != None]
    truth[inds] = 1
    return truth


def write_data(genes, rankings, gts, gttitles, title, atitle, resultspath):
    '''Writes the data to file'''
    # Get data and truths
    tofile = [genes, rankings]
    header = ['Genes', title]
    for gt in gts:
        tofile.append(truthVector(genes, gt))
    header += gttitles
    f = open("{}/data__{}_{}.txt".format(resultspath, title, atitle), 'w')
    f.write('\t'.join(header) + '\n')  # write header
    for data in zip(*tofile):
        f.write('\t'.join([str(x) for x in data]) + '\n')  # write data
    f.close()
    return "{}/data__{}_{}.txt".format(resultspath, title, atitle)

def save_xy_data(X,Y,title,atitle,resultspath,header='FPR\tTPR\n'):
    '''Writes the ROC curve data to file for plotting outside python'''
    f = open("{}/data__{}_{}.txt".format(resultspath, title,atitle), 'w')
    f.write(header)  # write header
    for i in range(0,len(X)):
        f.write('{}\t{}\n'.format(X[i],Y[i]))# write data
    f.close()
     

def assess_ranking(rankings, rgenesets, titles, gts, gttitles, atitle, resultspath,chosengt=2):
    '''Determine if ranking corresponds to ground truth data. '''

    ### Plotting Settings ###
    matplotlib.rcParams.update({'font.size': 18})
    COLORS = ['r', 'b', 'y', 'c', 'm', 'g', 'pink', 'navy',
              'darkorange', 'indigo', 'lightsteelblue', 'olive', 'darkgray',
              'tan', 'maroon', 'dodgerblue', 'thistle']
    colorcycle = cycle([COLORS[i] for i in range(0, len(titles))])
    print(len(titles), 'len of titles')

    tests = ["ROC", "PR"]
    
#     # combined gt
#     combinedgt = list(set([gene for sublist in gts for gene in sublist]))
#     combinedgt.sort()
#     combinedgt = np.array(combinedgt)
    gt=list(set(gts[chosengt]))
    gt.sort()
    gt=np.array(gt)
#     embed()
    for i in range(0, len(gts)):
        tmp = truthVector(gts[i], gt)
        print(('{} genes for {}'.format(tmp.shape[0], gttitles[i])))

    ### Determine true and false positives ###
    numgenes = []
    outfls=[]
    for test in tests: # For ROC and PR analysis
        plt.clf()
        plt.figure(1)
        plt.xscale('linear')

        if test == 'ROC':
            aucs = []
        colors = []

        for o in range(0, len(titles)): # Iterate over networks to analyze
            title = titles[o]
            rgenes = np.array(rgenesets[o])
            ranking = np.array(rankings[o])
#             rgenes=rgenes[ranking>0] # MeTeOR still performs best
#             ranking=ranking[ranking>0]
            
            ### Look at top 1000, because that's what we care about ###
            if 'rank' in title.lower() and 'diff' not in title.lower():
                ranking=1./ranking

            
            if 'diff' not in title.lower():
                rgenes, ranking = [list(x) for x in zip(
                *sorted(zip(rgenes, ranking), key=itemgetter(1), reverse=True))]
            else:
                rankingPred = np.array(rankings[o-2])
                rgenes, ranking, rankingPred = [list(x) for x in zip(
                *sorted(zip(rgenes, ranking,rankingPred), key=itemgetter(2), reverse=True))]

            rgenes = rgenes[:1000]
            ranking = ranking[:1000]
            

            if test == 'ROC':
                numgenes.append(len(rgenes))

            # Compare to combinedGT
            truth = truthVector(rgenes, gt)
            try:
                fpr, tpr, _ = roc_curve(truth, ranking, 1)
            except TypeError:
                print(title)
                embed()
                continue
            
            if test == 'ROC':
                # print titles[o]
                color = next(colorcycle)
                colors.append(color)
                # print titles[o], color
                plt.plot(fpr, tpr, color)
                _auc = auc(fpr, tpr)
                if sum(truth) == len(truth):
                    _auc = 1
                aucs.append(_auc)
                # Saves XY data for the ROC curves
                save_xy_data(fpr, tpr, title+'ROC', atitle, resultspath) 
            else:
                precision, recall, _ = precision_recall_curve(
                    truth, ranking, 1)
                precision=list(precision)
                recall=list(recall)
                precision.reverse()
                recall.reverse()
                # Saves XY data for the PR curves
                save_xy_data(recall, precision, title+'PR', atitle, resultspath,'Recall\tPrecision\n')
                
                color = next(colorcycle)
                colors.append(color)
                plt.plot(recall, precision, color)
            
            # Write all data to file
            outfl=write_data(rgenes, ranking, gts,
                       gttitles, title, atitle, resultspath)
            outfls.append(outfl)
        if 'ROC' in test:
            labels = [x + ' AUC={:0.2f}'.format(y)
                      for x, y in zip(titles, aucs)]
            plt.plot([0, 1], [0, 1], 'k--')
            labels += ['Random AUC=0.50']
            plt.xlabel("FPR", fontsize=22, fontweight='bold')
            plt.ylabel("TPR", fontsize=22, fontweight='bold')
        else:
            labels = titles
            plt.xlabel("Recall", fontsize=22, fontweight='bold')
            plt.ylabel("Precision", fontsize=22, fontweight='bold')
        plt.xlim([0, 1])
        plt.ylim([0, 1])

        plt.savefig("{}/AccuracyPlots/{}_{}.pdf".format(resultspath,
                                                        test, atitle))
        # with legend
        _ = plt.legend(labels, loc='best')
        plt.savefig("{}/AccuracyPlots/{}_{}_leg.pdf".format(resultspath,
                                                            test, atitle))
        plt.close()

    return outfls
