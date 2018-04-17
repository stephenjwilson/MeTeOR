import sys,traceback
import time,re,os,scipy,itertools
from operator import itemgetter
from IPython import embed
import pandas as pd
from collections import defaultdict
from scipy.sparse import coo_matrix
import numpy as np
import scipy.sparse
from lxml.etree import XMLSyntaxError

### Internal packages
from utils import setup,parse
from numpy import int32
from _operator import itemgetter

#from lxml.html import soupparser
#import lxml.etree as etree

def fixXML(fl):
    ''' Attempt to fix XML files that cause problems. For example,  <ForeName>M<?ForeName>  is a line from 14620000 and breaks the entire XML parser.'''
    #Delete lines that have symbols between the tags <>
    lines=open(fl).readlines()
    newfl=fl.replace('.xml','')+'fixed.xml'
    OUT=open(newfl,'w')
    fix=[]
    for i in range(0,len(lines)):
        line=lines[i]
        go=1
        valid=True
        while go:
            try:
                ind1=line.index('<')
            except ValueError:
                go=0
            if go:
                ind2=line.index('>')
                tag=line[ind1+1:ind2]
                tag=tag.replace('/','').split(' ')[0]
                valid = re.match('^[\w-]+$', tag) is not None
                if not valid:
                    fix.append(i)
                #print(tag,valid,line)
                line=line[ind2+1:]
        if valid:
            OUT.write(lines[i])
        #if i<2:
        #    OUT.write(lines[i])
    OUT.close()
    return newfl
        
        
def main(datadir='../data',storagedir='/home/stephen/ExtraDrive1/MeSH_V3',resultdir='.',verbose=1):
    '''
    Crawls downloaded MEDLINE data to create matrix of PMID by MeSH/SCR terms as well as the co-occurrence matrix
    '''       
    ### Run Crawl ###
    tic1=time.time()
    addedids=[]
    addedmesh=[]
    missed=0
    allinds=[]
    # Process all PMIDs
    step=5000

    examplegenelinks,enumpmid,enumui,uis=setup(datadir=datadir,resultdir=resultdir,verbose=verbose)

    count=len(enumpmid.keys())
    if verbose:
        print('IDs from file: {}'.format(count))

    # Contains metadata of pmids
    meta=defaultdict(dict)
    indicies=[]
    tic=time.time()
    log=open('Log.txt','w')
    for i in range(0,count,step):
        if i+step>count:
            artnum= '{}:{}'.format(i,count)
            retstart=i
            retend=count
        else:
            artnum= '{}:{}'.format(i,i+step)
            retstart=i
            retend=i+step

        fl='{}/{}.xml'.format(storagedir,i)
        if os.path.exists(fl):
            if os.stat(fl).st_size!=0:
                # Sometimes the tags are malformed. If that's the case, use HTML recovering to fix the tags in the parsing process
                try:
                    indicies,pmids,meshids,miss=parse('{}/{}.xml'.format(storagedir,i),enumpmid,enumui,uis,meta)
                except:
                    try:
                        indicies,pmids,meshids,miss=parse('{}/{}.xml'.format(storagedir,i),enumpmid,enumui,uis,meta,recover=True)
                        if verbose:
                            print('Recovering XML')
                    except:
                        if verbose:
                                print('Removing XML lines that are problematic')
                        try:
                            fl2=fixXML(fl) # Try removing lines that have symbols within tags
                            fl=fl2
                        except:
                            pass
                        try:
                            indicies,pmids,meshids,miss=parse('{}/{}.xml'.format(storagedir,i),enumpmid,enumui,uis,meta,recover=True)
                        except:
                            missed+=5000
                            e = traceback.format_exc()
                            s='XML File:{} has a problem. Skipping...Error:{}\n'.format(i,e)
                            log.write(s)
                            print(s)
                        continue
            else:
                os.remove(fl)
                s='Failed to load:{}..Continuing'.format(fl)
                log.write(s)
                print(s)
                missed+=5000
                continue
        else:
            s='Failed to load:{}..Continuing'.format(fl)
            log.write(s)
            print(s)
            missed+=5000
            continue
        
        # Write out UIs
        OUT=open('{}/UIs.txt'.format(resultdir),'w')
        OUT.write('\n'.join(uis))
        OUT.close()
        
        # Output indices
        OUT=open('{}/{}_IND.txt'.format(storagedir,i),'w')
        for x,y in indicies:
            OUT.write('{}\t{}\n'.format(x,y))
        OUT.close()
        
        allinds+=indicies
        missed+=miss
        addedids+=pmids
        addedmesh+=meshids
    
        if i%100000==0 and verbose and i!=0:
            toc=time.time()
            print(toc-tic, ' for {}'.format(artnum))
            print('Added PMIDs: {}'.format(len(addedids)))
            print('Added MeSH PMIDs: {}'.format(len(addedmesh)))
            print('Missed: {}'.format(missed))
            tic=time.time()
    log.close()
    toc=time.time()
    if verbose:
        print('{} seconds to build'.format((toc-tic1)))
        print('{} minutes to build'.format((toc-tic1)/60))

    ### Output data
    # Write out UIs
    OUT=open('{}/UIs.txt'.format(resultdir),'w')
    OUT.write('\n'.join(uis))
    OUT.close()

    # Write out Meta
    OUT=open('{}/meta.txt'.format(resultdir),'w')
    for key in meta:
        OUT.write('{}\t{}\n'.format(key,''.join('{}:{}'.format(key, val) for key, val in sorted(meta[key].items()))))
    OUT.close()
    
    # Write out all Indices
    OUT=open('{}/AllIND.txt'.format(storagedir),'w')
    for x,y in allinds:
        OUT.write('{}\t{}\n'.format(x,y))
    OUT.close()
    
    # Write out PMIDs
    OUT=open('{}/pmids.txt'.format(resultdir),'w')
    currentPMIDs=np.empty(shape=(len(enumpmid.keys()),1))
    for pmid, ind in enumpmid.items():
        currentPMIDs[ind]=pmid
    for i in range(0,len(currentPMIDs)):
        OUT.write('{}\t{}\n'.format(int(currentPMIDs[i][0]),i))
    OUT.close()
    currentPMIDs=[str(int(x[0])) for x in currentPMIDs]
    total=len(allinds)
    rows=len(currentPMIDs)
    cols=len(uis)
    
    # Write out matrix of indices
    mat=makeTermMat(rows, cols, total, storagedir, resultdir, verbose)
    ### Make CoOccurence Matrix
    makeCliques(mat,resultdir, currentPMIDs=currentPMIDs)
    # Time-stamped
    makeCliques(mat,resultdir,currentPMIDs=currentPMIDs,year=2014,meta=meta)
    
    toc=time.time()
    print('{} seconds to do everything'.format((toc-tic1)))
    print('{} minutes to do everything'.format((toc-tic1)/60))
    
def makeTermMat(rows,cols,total,storagedir,resultdir,verbose=1):
    mat=coo_matrix((rows,cols))
    step=10000000
    allinds=open('{}/AllIND.txt'.format(storagedir))
    for i in range(0,total,step):
        x=[]
        y=[]
        while len(x)<step and i+len(x)!=total:
            line=allinds.readline()
            line=line.strip().split('\t')
            x.append(int(line[0]))
            y.append(int(line[1]))
            if x[-1]>=rows or y[-1]>=cols:
                embed()
                 
        chunk=coo_matrix((np.ones(len(x)),(x,y)),shape=(rows,cols))
        mat+=chunk

        if verbose:
            print('At {}-{} of {}'.format(i,i+step,total))
        
    scipy.sparse.save_npz('{}/pmid_ui'.format(resultdir),mat)
    return mat

def makeCliques(mat,resultdir,currentPMIDs=None,year=None,meta=None,step=100000,verbose=1):
    '''
    Creates the co-occurrence matrix from the PMID by MESH/SCR term matrix
    '''
    
    # if currentPMIDs ==None, load it
    if currentPMIDs==None:
        f=open('{}/pmids.txt'.format(resultdir))
        currentPMIDs=[]
        for line in f:
            line=line.strip('\n').split('\t')
            currentPMIDs.append(line[0])
        
    
    # if subsetting term-pmid matrix by year, do it now
    if year!=None:
        if meta==None:
            # Load meta from file
            meta={}
            f=open('{}/meta.txt'.format(resultdir))
            for line in f:
                line=line.strip('\n').split('\t')
                pmid=line[0]
                d=line[1].split(',')
                for part in d:
                    key,val=part.split(':')
                    try:
                        meta[pmid][key]=val
                    except KeyError:
                        meta[pmid]={key:val}
        # Get PMIDs of all relevant articles
        go=[]
        for pmid in meta.keys():
            d=meta[pmid]
            if 'Year' in d:
                y=int(d['Year'])
                if y<year:
                    go.append(pmid)
        # subset matrix
        go=np.array(go)
        subset=np.in1d(currentPMIDs, go)
        mat=mat[subset,:]

    h,w=mat.shape    
    if verbose:
        print('Starting Dot product to make co-occurrence matrix...')
    tic=time.time()
    outMat=np.dot(mat.transpose(),mat)
    toc=time.time()
    if verbose:
        print('Created co-occurrence matrix with {} elements in {} seconds. Saving...'.format(outMat.nnz, toc-tic))
    if year==None:
        scipy.sparse.save_npz('{}/coOccurrence'.format(resultdir),outMat)
    else:
        scipy.sparse.save_npz('{}/coOccurrence_{}'.format(resultdir,year),outMat)
    return outMat

if __name__=='__main__':
    #sys.settrace(main)
    main(datadir='../data',resultdir='../results',storagedir='/home/stephen/ExtraDrive1/MeSH_V3')
