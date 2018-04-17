'''
Created on May 23, 1817

@author: daniel
'''
import os,copy
import csv
import operator
# from sys import argv
import numpy as np
import pandas as pd
from numpy import arange
### Graphviz Libraries ###
from graphviz import Graph
from pygraphviz import *
###                    ###

import seaborn as sns
from seaborn import heatmap
from math import sin, cos, pi
# from itertools import combinations
import matplotlib.pyplot as plt
import  matplotlib.colors as colorsfunc
from IPython import embed
import networkx as nx
import scipy.stats as spstats
import readline # workaround for importr problem
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
from matplotlib.colors import rgb2hex
import itertools

try:
    from networkx.drawing.nx_agraph import graphviz_layout
    layout=graphviz_layout
except ImportError:
    print("PyGraphviz not found; drawing with spring layout; will be slow.")
    layout=nx.spring_layout


def importSupportFile(filePath):
    fileHandle = open(filePath, 'rb')
    reader = csv.reader(fileHandle, delimiter='\t')
    indexDict = {}
    dataDict = {'MeTeOR': {}, 'CoReg': {}, 'IPMS': {},
                'Pathway': {}, 'STRINGEXP': {}}
    stop=52
    for line in reader:
        if(len(indexDict) == 0):
            indexDict['Genes'] = line.index('Genes')
            try:
                indexDict['MeTeOR'] = line.index('MeTeOR Weight')
            except ValueError:
                try:
                    indexDict['MeTeOR'] = line.index('MeTeOR NMF Pred')
                except ValueError:
                    indexDict['MeTeOR'] = line.index('Diff Ranking')
                    stop=0
                #indexDict['MeTeOR'] = line.index('MeTeOR Diff Weight')
            indexDict['CoReg'] = line.index('CoReg')
            indexDict['IPMS'] = line.index('IPMS')
            indexDict['Pathway'] = line.index('Pathway')
            indexDict['STRINGEXP'] = line.index('STRINGEXP')
        elif reader.line_num==stop and stop!=0:#(len(indexDict)==30):
            break
        else:
            if line[indexDict['Genes']] in list(dataDict['MeTeOR'].keys()):
                    continue
            dataDict['MeTeOR'][line[indexDict['Genes']]] = float(
                line[indexDict['MeTeOR']])
            dataDict['CoReg'][line[indexDict['Genes']]] = float(
                line[indexDict['CoReg']])
            dataDict['IPMS'][line[indexDict['Genes']]] = float(
                line[indexDict['IPMS']])
            dataDict['Pathway'][line[indexDict['Genes']]] = float(
                line[indexDict['Pathway']])
            dataDict['STRINGEXP'][line[indexDict['Genes']]] = float(
                line[indexDict['STRINGEXP']])

    return dataDict

def mapHGNC():
    clean={}
    data=open('../../../data/Mappings/Mapping_Ccgdd16Sx91Extended.tsv').readlines()
    GENES=open('../../../data/Mappings/HGNCMap_10-27-16.txt').readlines()
    GENES=[x.strip().split('\t')[1] for x in GENES]
    GENES={x:1 for x in GENES}

    for line in data:
        line=line.strip().split('\t')
        #dic[line[0]]=(line[1],line[2])
        g=None
        m=''
        go=0
        for part in [line[1]]+line[2].replace('protein','').strip().replace(',',';;').split(';;'):
            if 'MESH' in part:
                m=part.strip()
            try:
                GENES[part.strip()]
                g=part.strip()#clean[line[0].strip()]=part.strip()
                go=1
            except:
                pass
        #clean[m]=g
        if go:
            clean[line[0]]=g
    return clean

def plotNetworkGeneric(dataDict, gts, name, styles, colors):
    g = Graph(name, strict=False, engine='neato', format='svg')
    g.attr('node', fontsize='22', fontname='Arial bold', margin='0.025,0.025')
    g.attr('edge', penwidth='3')
    g.node('EGFR', pos="0,0!")
    nodes = [x[0] for x in sorted(list(dataDict['MeTeOR'].items()),
                                  key=operator.itemgetter(1), reverse=True)]
#     nodes = sorted(dataDict['MeTeOR'].keys())
    nodeNum = len(nodes)
    radius = nodeNum / 4.0
    for i in arange(nodeNum):
        nodeName = nodes[i]
        y = cos(((pi * 2) / nodeNum) * i) * radius
        x = sin(((pi * 2) / nodeNum) * i) * radius
        g.node(nodeName, pos="{},{}!".format(x, y))
#         g.edge('EGFR', nodeName, color="#7d7d7d", style='solid')
        g.edge('EGFR', nodeName, color="#000000", style='solid')
        for gt in gts:
            if(dataDict[gt][nodeName] > 0):
                if(styles):
                    gtStyle = styles[gt]
                else:
                    gtStyle = 'solid'
                if(colors):
                    gtColor = colors[gt]
                else:
                    gtColor = 'black'
                if(gtStyle == 'tapered'):
                    g.edge('EGFR', nodeName, color=gtColor, style=gtStyle,
                           dir='both', arrowhead='none', arrowtail='none')
                else:
                    g.edge('EGFR', nodeName, color=gtColor, style=gtStyle)
    g.render(g.name)


def plotTableView(dataDict, gts, name, colors):
    g = Graph(name, format='svg')
    nodes = [x[0] for x in sorted(list(dataDict['MeTeOR'].items()),
                                  key=operator.itemgetter(1))]
    for i in range(len(gts)):
        gt = gts[i]
        with g.subgraph(name='cluster_{}'.format(i)) as sub:
            for n in nodes:
                nodeName = '{}-{}'.format(gt, n)
                if(dataDict[gt][n] > 0):
                    sub.node(name=nodeName, label=n,
                             style='filled', fillcolor=colors[gt])
                else:
                    sub.node(name=nodeName, label=n)
                if(i >= 1):
                    prevGT = gts[i - 1]
                    sub.edge('{}-{}'.format(prevGT, n), nodeName,
                             style='invis')
            sub.attr(label=gt, color='black')
    g.render(g.name)


def plotTableViewTranspose(dataDict, gts, name, colors):
    g = Graph(name, format='svg')
    nodes = [x[0] for x in sorted(list(dataDict['MeTeOR'].items()),
                                  key=operator.itemgetter(1))]
    for i in range(len(nodes)):
        n = nodes[i]
        with g.subgraph(name='cluster_{}'.format(i)) as sub:
            for gt in gts:
                nodeName = '{}-{}'.format(gt, n)
                if(dataDict[gt][n] > 0):
                    sub.node(name=nodeName, label=n,
                             style='filled', fillcolor=colors[gt])
                else:
                    sub.node(name=nodeName, label=n)
                if(i >= 1):
                    prevNode = nodes[i - 1]
                    sub.edge('{}-{}'.format(gt, prevNode), nodeName,
                             style='invis')
            sub.attr(color='white')
    g.render(g.name)


def plotHeatMap(dataDict, gts, name,diff={}):
    nodes = [x[0] for x in sorted(list(dataDict['MeTeOR'].items()),
                                  key=operator.itemgetter(1),reverse=True)]
    
    genes = []
    genes2 = []
    cols = []
    vals = []
    ugenes=[]
    ugenes2=[]
    for i in range(len(nodes)):
        for j in range(len(gts)):
            if diff!={}:
                genes2.append('{} {}'.format(nodes[i],diff['MeTeOR'][nodes[i]]))
            
            genes.append('{}'.format(nodes[i]))
            cols.append(gts[j])
            vals.append((0, j + 1)[dataDict[gts[j]][nodes[i]] > 0])
        ugenes.append(genes[-1])
        if diff!={}:
            ugenes2.append(genes2[-1])
    df = pd.DataFrame({'Genes': np.array(genes),
                       'Ground Truths': np.array(cols),
                       'Supported': np.array(vals)})
    df['Ground Truths'] = pd.Categorical(df['Ground Truths'], gts)
    df['Genes'] = pd.Categorical(df['Genes'], ugenes)
    df = df.pivot("Genes", "Ground Truths", "Supported")
    plt.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Arial'],
                      'weight': 'bold'})
    ax = heatmap(data=df, xticklabels=True, yticklabels=True, cmap='gist_yarg_r',
                 cbar=False, linecolor="black", linewidths=0.5, square=True)
    ax.xaxis.tick_top()
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    plt.savefig(name + '.svg', bbox_inches='tight')
    if diff!={}:
        df = pd.DataFrame({'Genes': np.array(genes2),
                           'Ground Truths': np.array(cols),
                           'Supported': np.array(vals)})
        df['Ground Truths'] = pd.Categorical(df['Ground Truths'], gts)
        df['Genes'] = pd.Categorical(df['Genes'], ugenes2)
        df = df.pivot("Genes", "Ground Truths", "Supported")
        plt.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Arial'],
                          'weight': 'bold'})
        ax = heatmap(data=df, xticklabels=True, yticklabels=True, cmap='gist_yarg_r',
                     cbar=False, linecolor="black", linewidths=0.5, square=True)
        ax.xaxis.tick_top()
        ax.xaxis.set_ticks_position('top')
        ax.xaxis.set_label_position('top')
        plt.yticks(rotation=0)
        plt.xticks(rotation=90)
        plt.savefig(name + 'rank.svg', bbox_inches='tight')


def plotPanel(name,elists,nlists,ecolors,ncolors,edgeweights,nodeshapes,penwidth='10'):
    # Note that nodeshapes has been appropriated for color now
    G=AGraph()

    if 'No' not in name:# don't add EGFR in the no support plot
        G.add_node('EGFR',pos="0,0!",fontsize='150',color='#000000',style='filled',penwidth=penwidth,fillcolor='#000000',fontcolor='#FFFFFF',shape='box')#nodeshapes['EGFR'])

    for nlist,ncolor in zip(nlists,ncolors):
        for node in nlist:
            if node=='EGFR':
                continue
            #print node,nodeshapes[node]
            if nodeshapes[node]!='#FFFFFF':
                fontcolor='#FFFFFF'
            else:
                fontcolor='#000000'
            G.add_node(node, penwidth=penwidth,color='#000000',fillcolor=nodeshapes[node],style='filled',fontsize='100',fontcolor=fontcolor)#,shape=nodeshapes[node]) 

    for edge,ecolor,weight in zip(elists,ecolors,edgeweights):
        if ecolor=='#FFFFFF':
            weight=0
        G.add_edge(edge[0],edge[1], color=ecolor,penwidth="{}".format(weight))
        G.add_edge(edge[1],edge[0], color=ecolor,penwidth="{}".format(weight))

    G.draw('Overlap_Panel{}.svg'.format(name),prog='twopi',
           args='-Goutputorder=edgesfirst -Goverlap=false')

    G.draw('Overlap_Panel{}NEATO.svg'.format(name),prog='neato',
           args='-Goutputorder=edgesfirst -Goverlap=false')
    
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def calcAllPaths(dataDicts,gts,enrichpaths,pathtogene,finalnodes,qs,nodes,nodeshapes):
    ### Map Colors ###
    nlists=[]
    ncolors=[]
    
    elistsgt=[]
    edgeweightsgt=[]
    ecolorsgt1=[]
    ecolorsgt2=[]
    
    ecolorspath=[]
    edgeweightspath=[]
    elistspath=[]
    
    allnodes=[]
    allnodes2=[]
    
    ### Get the ground-truth data ###
    colors=itertools.cycle([rgb2hex(x)                        
                 for x in [(1,0,0),(1,0,0),(1,0,0)]])
    for gt in gts:
        if gt=='Pathway':
            continue
        for dataDict in dataDicts:
            gtnodes = [x[0] for x in sorted(list(dataDict[gt].items()),
                key=operator.itemgetter(1), reverse=True) if int(x[1])==1]
            allnodes+=gtnodes
            newedges=[(x,'EGFR') for x in gtnodes]
            elistsgt+=newedges
            edgeweightsgt+=[10 for x in range(0,len(newedges))]
    ecolorsgt1+=[rgb2hex((1,0,0))  for x in elistsgt]  
    ecolorsgt2+=['#FFFFFF' for x in elistsgt]
    
    ### Add in the enriched pathways data ###
    cmin=0.25
    cmax=0.75
    colors=itertools.cycle([rgb2hex(x)                        
                 for x in plt.cm.gray(np.linspace(cmin,cmax,len(enrichpaths)))])
    # plot colormap
    new_cmap = colorsfunc.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=plt.cm.gray.name, a=cmin, b=cmax),
        plt.cm.gray(np.linspace(cmin, cmax, 100)))
    gradient = np.linspace(0, 0.1, 90)
    gradient = np.vstack((gradient, gradient))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.imshow(gradient, aspect='auto', cmap=new_cmap)
    plt.colorbar(ticks=np.linspace(0, 0.1, 9))
    plt.ylim(0,0.1)
    ax.set_axis_off()
    plt.savefig('Colormap.pdf')
    
    for o in range(0,len(enrichpaths)):
        #genes=pathtogene[enrichpaths[o]]
        ns=finalnodes[o]
        #nodels=list(set(ns).intersection(genes))
        newedges=[x for x in itertools.combinations(ns,2)]
        elistspath+=newedges
        edgeweightspath+=[((1.-qs[o])*100 -90)*.3 for x in range(0,len(newedges))]#[10 for x in range(0,len(newedges))]#
        allnodes2+=ns

    while len(ecolorspath)< len(elistspath):
        color=next(colors)
        ecolorspath.append(color)

    ### Determine node colors ###
    allnodes2=list(set(allnodes2))
    
    allsupport=list(set(allnodes).intersection(allnodes2))
    onlypath=list(set(allnodes2).difference(allnodes))
    onlyexp=list(set(allnodes).difference(allnodes2))
    nosup=list(set(nodes).difference(set(allnodes).union(allnodes2)))

    # Get all genes
    allnodes=list(set(allnodes+allnodes2))
    try:
        allnodes.pop(allnodes.index('EGFR'))
    except:
        pass
    # Color genes
    nlists=allsupport,onlypath+onlyexp
    ncolors1=['' for x in nlists]
    #ncolors2=ncolors+['#7ca1ea','#7ca1ea','#7cea87']#,'#FFFFFF']

    ### Plot networks ###
    plotPanel('allPathways',elistspath+elistsgt,nlists,
              ecolorspath+ecolorsgt2,ncolors1,
              edgeweightspath+edgeweightsgt,nodeshapes)
    
    plotPanel('all',elistspath+elistsgt,nlists,
              ecolorspath+ecolorsgt1,ncolors1,
              edgeweightspath+edgeweightsgt,nodeshapes)
    plotPanel('NoPathwaySupport',[],[nosup],[],['#FFFFFF'],[],nodeshapes)

    return elistspath+elistsgt,nlists,ecolorspath+ecolorsgt1,ncolors1,edgeweightspath+edgeweightsgt,nosup

def calcEachPath(enrichpaths,pathtogene,finalnodes,qs,nodes,nodeshapes,data):
    elists,nlists,ecolors,ncolors,edgeweights,nothit=data # unpack

    ### Map Colors ###
    colors=itertools.cycle([rgb2hex(x)                        
                 for x in plt.cm.gray(np.linspace(0.25,0.75,len(enrichpaths)))])
    for i in range(0,len(enrichpaths)):
        genes=pathtogene[enrichpaths[i]]
        ns=finalnodes[i]
        nlistsi=nlists
        ncolorsi=ncolors
        newedges=[x for x in itertools.combinations(ns,2)]
        elistsi=elists+newedges
        color=next(colors)
        ecolorsi=['#FFFFFF' for x in elists]+[ color for x in newedges]
        edgeweightsi=edgeweights+[10 for x in range(0,len(newedges))]
        plotPanel('{}_{}_{:0.3f}'.format(i,enrichpaths[i],qs[i]),elistsi,nlistsi,ecolorsi,ncolorsi,edgeweightsi,nodeshapes)


def calcGTPaths(dataDicts,gts,enrichpaths,pathtogene,finalnodes,qs,nodes,nodeshapes,data):
    elists,nlists,ecolors,ncolors,edgeweights,nothit=data # unpack
    colors=itertools.cycle([rgb2hex(x)                        
                 for x in [(1,0,0),(1, 0, 0),(1,0,0)]])
    for gt in gts:# Plot per gt
        if gt=='Pathway':
            continue
        nlistsi=[]
        ncolorsi=[]
        elistsi=elists
        ecolorsi=['#FFFFFF' for x in elists]
        edgeweightsi=[0 for x in elists]
        color=next(colors)
        for dataDict in dataDicts:
            gtnodes = [x[0] for x in sorted(list(dataDict[gt].items()),
                    key=operator.itemgetter(1), reverse=True) if int(x[1])==1]
            newnodes=gtnodes#[list(set(gtnodes).difference([x for ls in nlists for x in ls]))]
            nlistsi+=nlists
            ncolorsi+=ncolors
            newedges=[(x,'EGFR') for x in newnodes]
            elistsi+=newedges
            
            ecolorsi+=[color for x in newedges]
            edgeweightsi+=[10 for x in newedges]
            #print len(edgeweightsi),len(ecolorsi),len(elistsi)
        plotPanel('gt{}'.format(gt),elistsi,nlistsi,ecolorsi,ncolorsi,edgeweightsi,nodeshapes)
    
def plotPathwayPanels(dataDicts, gts, name, net,pathways):
    ### Determine the initial nodes of interest
    nodeshapes={}
    allnodes=[]

    shapes=['#FFFFFF','#848484']#'oval','diamond','trapezium','pentagon']
    
    for i in range(0,len(dataDicts)):
        nodes = [x[0] for x in sorted(list(dataDicts[i]['MeTeOR'].items()),
                                  key=operator.itemgetter(1), reverse=True)]
        nodes=sorted(nodes)
        allnodes+=nodes

        for node in nodes:
            if node not in nodeshapes: # do not overwrite shape
                nodeshapes[node]=shapes[i]

    allnodes=sorted(list(set(allnodes)))
    allnodes+=['EGFR']
    #nodeshapes['EGFR']='box'
    ### Determine pathway enrichments relationships ###
    qs,enrichpaths,finalnodes=enrichAll(allnodes,pathways,name,shrinking=1)

    ### Calc All Paths ###
    data=calcAllPaths(dataDicts,gts,enrichpaths,pathways,finalnodes,qs,allnodes,nodeshapes)

    ### Calc Each Path ###
    calcEachPath(enrichpaths,pathways,finalnodes,qs,allnodes,nodeshapes,data)
    
    ### Now plot GT support ###
    calcGTPaths(dataDicts,gts,enrichpaths,pathways,finalnodes,qs,allnodes,nodeshapes,data)

    
def readNetwork(netfl,clean,skipchar='#',sep='\t'):
    '''Requires a path to an edge list of a network information, and outputs a networkx graph instance'''
    lines=open(netfl).readlines()
    keep=[]

    for i in range(0,len(lines)):
        if lines[i][0]==skipchar:
            continue
        line=lines[i]
        line=line.strip().split(sep)
        meta=line[-1].split(';;')
        label1=[x.split('::')[1] for x in meta if 'Entity1Label' in x][0]
        label2=[x.split('::')[1] for x in meta if 'Entity2Label' in x][0]

        # Map Genes
        try:
            line[0]=clean[line[0]]
            line[1]=clean[line[1]]
        except KeyError:
            continue
            
        lines[i]='\t'.join([line[0],line[1],label1,label2,line[2]])
        if float(line[2])>1:
            keep.append(i)
    lines=list(np.array(lines)[keep])
    net=nx.parse_edgelist(lines, nodetype = str, delimiter=sep,data=(('label1',str),('label2',str),('weight',float)))
    return net

def readPathway(pathfl):
    '''Reads a pathway into a pair of dictionaries'''
    lines=open(pathfl).readlines()
    pathtogene={}
    genetopath={}
    for path in lines:
        path=path.strip().split('\t')
        pathtogene[path[0]]=path[1:]
        for gene in path[1:]:
            genetopath[gene]=path
    return genetopath,pathtogene

def enrichAll(mynodes,pathtogene,name,qthresh=0.1,shrinking=1):
    nodes=copy.deepcopy(mynodes)
    enrichpaths=[]
    ps=[]
    genepool=len([x for ls in list(pathtogene.items()) for x in ls])
    stats = importr('stats')

    finalnodes=[]
    for path in pathtogene:
        pathwayg=pathtogene[path]
        inter=len(set(nodes).intersection(pathwayg))
        p=spstats.hypergeom.sf(inter- 1,genepool,len(pathwayg),len(nodes))
        if len(nodes)<1 or inter==0: # don't analyze things when I shouldn't
            continue
        ps.append(p)
        enrichpaths.append(path)
        finalnodes.append(set(nodes).intersection(pathwayg))

    qs = np.array(stats.p_adjust(FloatVector(ps), method = 'BH'))
    ps = np.array(ps)
    enrichpaths=np.array(enrichpaths)
    finalnodes=np.array(finalnodes)

    enrichpaths=enrichpaths[np.where(qs<qthresh)]
    finalnodes=finalnodes[np.where(qs<qthresh)]
    ps=ps[qs<qthresh]
    qs=qs[qs<qthresh]
    
    qs,enrichpaths,finalnodes=list(zip(*sorted(zip(qs,enrichpaths,finalnodes))))
    qs=list(qs)
    ps=list(ps)
    enrichpaths=list(enrichpaths)
    finalnodes=list(finalnodes)

    if shrinking:
        topop=[]
        alreadylinked=[]
        for path in enrichpaths:
            pathwayg=pathtogene[path]
            new=set(nodes).intersection(pathwayg).difference(alreadylinked)
            alreadylinked+=list(new)
            if len(new)==0:
                topop.append(path)

        for path in topop:
            qs.pop(enrichpaths.index(path))
            ps.pop(enrichpaths.index(path))
            finalnodes.pop(enrichpaths.index(path))
            enrichpaths.pop(enrichpaths.index(path))

    ### Write to file ###
    OUTFILE=open('PathwayEnrichment_{}'.format(name),'w')
    OUTFILE.write('Pathway\tq-value\tp-value\tIntersection\t# of MeTeOR\t # of Pathway\tOverlappingGenes\tMeTeORNodes\tPathwayNodes\n')
    mnum=len(mynodes)
    for data in zip(enrichpaths,qs,ps,finalnodes):
        path,q,p,ns=data
        inter=len(ns)
        pathnodes=pathtogene[path]
        pnum=len(pathnodes)
        pathnodes=','.join(pathnodes)
        ns=','.join(ns)
        OUTFILE.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(path,q,p,inter,mnum,pnum,ns,','.join(mynodes),pathnodes))        
    return qs,enrichpaths,finalnodes


def runPlots():
    MeTeORFile='../data__MeTeOR Weight_area_gradientnorm.txt'#,
    MeTeORPredFile='../data__MeTeOR NMF Pred_area_gradientnorm.txt'#,
    MeTeORDiffFile='../data__Diff Ranking_area_gradientnorm.txt'#):
    
    '''Creates the Figure 4 Heatmaps. Currently outputs in SVGs. 
    Various methods exist to convert from SVG.'''
    root='../Results/NetworkPlots/'
    try:
        os.mkdir(root)
    except:
        pass
    if 'NetworkPlots' not in os.getcwd():
        os.chdir(root)
    styles = None
    colors = {'MeTeOR': '#000000', 'CoReg': '#00FF00', 'IPMS': '#FFAA00',
              'Pathway': '#A501FF', 'Combined': '#0246FF',
              'STRINGEXP': '#FF0D00'}
#     gts = ['Pathway', 'CoReg', 'IPMS', 'STRING', 'Combined']
    gts = ['Pathway', 'STRINGEXP', 'CoReg', 'IPMS']

    ### Get Mapping Data ###
    clean=mapHGNC()

    ### Get MeTeOR Gene-Gene Data ###
    netfl='/home/stephen/Documents/LiClipse Workspace/Networks/NetworkData/MeTeOR/MeTeORallGene-Genemapped.txt'
    pathfl='../../Data/GSEA.tsv'
    net=readNetwork(netfl,clean)
    genetopath,pathtogene=readPathway(pathfl)

    ###########################################################################

    MeTeORData = importSupportFile(MeTeORFile)
    MeTeORPredData = importSupportFile(MeTeORPredFile)
    #plotNetworkGeneric(MeTeORData, gts, '-'.join(gts), styles, colors)
    #plotTableView(MeTeORData, gts, '-'.join(gts) + '_Table', colors)
    #plotHeatMap(MeTeORData, gts, '-'.join(gts) + '_Heatmap')
    ###########################################################################


    ### Plot the overlay ###
    #plotNetworkOverlay(MeTeORData, '-'.join(gts), 'Overlay', styles,
    #                colors,net,pathtogene)

    plotPathwayPanels([MeTeORData,MeTeORPredData],gts,'OverlayPanels',net,pathtogene)

    ### Plot other forms ###
    #MeTeORPredData = importSupportFile(MeTeORPredFile)
    #plotNetworkGeneric(MeTeORPredData, gts, '-'.join(gts) + '_Pred', styles,
    #                   colors)
    #plotTableView(MeTeORPredData, gts, '-'.join(gts) + '_Table_Pred', colors)
    #plotHeatMap(MeTeORPredData, gts, '-'.join(gts) + '_Heatmap_Pred',diff=MeTeORDiffData)

    #MeTeORDiffData = importSupportFile(MeTeORDiffFile)
if __name__ == '__main__':
    runPlots()
