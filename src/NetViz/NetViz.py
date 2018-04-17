import networkx as nx
from networkx.drawing.nx_agraph import write_dot
from IPython import embed
import graphviz
import scipy.sparse
import numpy as np
from collections import Counter
from colorutils import Color
import itertools
 
def cmd_interface(program,xargin):
    '''Runs shell commands from python using subprocess.'''
    from subprocess import Popen,PIPE
    #print [program]+xargin
    p=Popen([program]+xargin,stdout=PIPE, stderr=PIPE)
    out,error=p.communicate()
    print(error)
    return out,error

def ReadNetwork(root,fl='coOccurrence.npz', cutoff=100,verbose=1):
    '''Reads the network and outputs a networkx graph instance'''
    mat=scipy.sparse.load_npz('{}/{}'.format(root,fl)).tocsr()
    labels=[]
    lines=open('{}/UIsMapped.txt'.format(root)).read().split('\n')
    for i in range(0,len(lines)):
        line=lines[i].split('\t')
        if len(line)>1:
            labels.append(line[1])
    labels=np.array(labels)    
    # eliminate diagonal as it stores the degree
    mat.setdiag(0)
    mat.eliminate_zeros()
    
    # Eliminate nodes with Unknown labels
    mask=labels!='Unknown'
    mat=mat[:,mask][mask,:]
    labels=labels[mask]
    
    # eliminate edges that are low confidence
    matShape=mat.shape
    nonzero=(mat>cutoff).nonzero()
    vals=mat[nonzero]
    mat=scipy.sparse.coo_matrix((vals.A1,nonzero),shape=matShape).tocsr()

    # Get rid of low degree nodes
    while np.sum(mat.astype(bool).sum(axis=1)<4)>0: # while loop as some <4 degree nodes are linked to degree 5 nodes
        degreeMask=mat.astype(bool).sum(axis=1)>=4
        mat=mat[:,degreeMask.A1][degreeMask.A1,:]
        labels=labels[degreeMask.A1]

    nonzero=scipy.sparse.tril(mat,k=-1).nonzero()# use only the lower triangle minus the diagonal
    
    ind1=nonzero[0]#[mask.A1]
    ind2=nonzero[1]#[mask.A1]
    
    if verbose:
        print('Iterating over {} edges'.format(ind1.shape))

    c=Counter(labels)
    # Make network
    net=nx.Graph()
    #nonNormnet=nx.Graph()
    
    colorLs=['Chemical','Disease','Gene']
    colors=[Color(rgb=[230,159,0]),Color(rgb=[86,180,223]),Color(rgb=[0,158,115]),Color(rgb=[240,228,66]),Color(rgb=[0,114,178]),Color(rgb=[213,94,0])] # color blind safe palette
    colorcomb=list(itertools.combinations_with_replacement(colorLs,2))
    #factors=[0.05,0.2,1,0.2,1,1]
    #print(colorcomb)
    for i,j in zip(ind1,ind2):
        w=mat[i,j]
        l1=labels[i]
        l2=labels[j]
        net.add_node(i,style='invis')
        net.add_node(j,style='invis')
        #nonNormnet.add_node(i,style='invis')
        #nonNormnet.add_node(j,style='invis')

        try:
            color=colors[colorcomb.index((l1,l2))]
            #factor=100/factors[colorcomb.index((l1,l2))]
        except:
            color=colors[colorcomb.index((l2,l1))]
            #factor=100/factors[colorcomb.index((l2,l1))]

        net.add_edge(i,j,weight=w,label1=l1,label2=l2,color='{}'.format(color.hex))
        #nonNormnet.add_edge(i,j,weight=w/100,label1=l1,label2=l2,color='{}'.format(color.hex))

    return net

def main(root,title=''):
    net=ReadNetwork(root) # from Networks folder

    #print "Writing dot..."
    write_dot(net,'{}/MeTeOR{}.dot'.format(root,title))
    #write_dot(nonNormnet,'{}/MeTeOR_nonNorm{}.dot'.format(root,title))
    render(root,title)
    
def render(root,title=''):
    '''Renders the dot'''
    try:
        fl=graphviz.render('sfdp','svg','{}/MeTeOR{}.dot'.format(root,title)) # I used ~/bin/sfdp for a local install and ran from terminal
    except Exception as e: 
        print(e)
        print('Make sure Graphviz installed')
        print('Alternatively, run /path/to/sfdp -Tsvg -o../results/MeTeOR.svg ../results/MeTeOR.dot')
    # terminal command: ~/bin/sfdp -Tsvg -o../results/MeTeOR.svg ../results/MeTeOR.dot
    # Convert to png with: inkscape -z -e MeTeOR.png -w 1024 -h 1024 MeTeOR.svg