import scipy.sparse
import numpy as np
from IPython import embed
from collections import defaultdict,Counter
import os
# 
# def getGeneNames(fl='HGNC_2-8-18.txt',datadir='../data'):
#     gtoName={}
#     lines=open('{}/mapping/gene/{}'.format(datadir,fl)).readlines()
#     for line in lines[1:]:
#         line=line.strip('\n').split('\t')
#         gtoName[line[5]]=line[1]
#     return gtoName

def genMapping(datadir='../data',resultdir='.',humanspecific=1,mode='any'):
    '''Maps to the best UIs based on NCBI's gene2pubmed file. Can scan the highest MeSH term related to the gene of pubmed's articles or any term on those articles'''
    # Load Relevant information
    gene2Pubmed=getGene2Pubmed(datadir, humanspecific) # Load NCBI's Entrez gene annotations on PMIDs
    pmidUi,uis,enumpmid,meshNames,nmids=getPmidUi(resultdir,datadir) # Load PMID-MeSH UI annotation
    geneSymbol,geneSyn=getGeneInfo(datadir) # Entrez Symbol and Synonym information
    geneRef=getGene2Refseq(datadir)
    
    ### Get most occurring terms
    print('Loaded Data. Making Mapping')
    terms={}
    meta={}
    
#     gtoName=getGeneNames(datadir=datadir)
    
    for gene in gene2Pubmed.keys():
        # Subset CoOccurence matrix to relevant articles for each gene
        inds=[]
        for pmid in gene2Pubmed[gene]:
            try:
                ind=enumpmid[pmid]
            except:
                continue
            inds.append(ind)
        if inds==[]:
            continue
        articles=np.array(inds)
        subset=pmidUi[articles,:]
        
        # Sum subset for most common terms
        if mode=='max':
            collation=np.sum(subset,axis=0)
            mval=np.max(collation)
            if mval==0: # don't map if no protein MeSH are present
                continue
            terminds=np.where(collation==mval)[1]
        elif mode=='any':
            collation=np.sum(subset,axis=0)
            terminds=collation.nonzero()[1]
        else:
            raise('Please chose a valid mode: any or max')
        tmpuis=uis[terminds]
        names=[meshNames[ui] for ui in tmpuis]
        # Check if names match known synonyms
        synOverlap=np.in1d(names,geneSyn[gene]) # Check exact matches on synonyms
        
        # Check if symbol is in the MeSH Names ( Maybe add entry terms here)
        symbol=geneSymbol[gene]
        if len(symbol)>2:
            symbolOverlap=[]
            for x in names:
                if 'fusion' in x:
                    symbolOverlap.append(False)
                else:
                    symbolOverlap.append(symbol in x)
        else:
            symbolOverlap=[False for x in names]
        
        # Check if there is an nmid
        genenmids=geneRef[gene]
        uinmids=[]
        for ui in tmpuis:
            try:
                uinmids.append(nmids[ui])
            except KeyError:
                uinmids.append('NOMATCH')
        
        nmOverlap=np.in1d(uinmids,genenmids) # Check exact matches on nmids

        # If the name matches the synonyms, the symbol, or the NMID, pick that term
        union=synOverlap+symbolOverlap+nmOverlap
        for i in range(0,len(tmpuis)):
            if union[i]:
                ui=tmpuis[i]
                if ui not in terms.keys():
                    terms[ui]=[]
                    meta[ui]={}

                terms[ui].append(gene)
                # Record how match was determined
                # Setup meta
                if gene not in meta[ui].keys():
                    meta[ui][gene]={}
                # Add Symbol and Weight
                meta[ui][gene]['Weight']='{}'.format(collation[0,terminds[i]])
                meta[ui][gene]['Symbol']='{}'.format(symbol)
                # Extract matching mapping method
                method=[]
                if nmOverlap[i]:
                    method.append('NMMap:{}/{}'.format(geneSymbol[gene],uinmids[i]))
                else:
                    method.append('')
                if symbolOverlap[i]:
                    method.append('SymbolMap:{}/{}'.format(geneSymbol[gene],names[i]))
                else:
                    method.append('')
                if synOverlap[i]:
                    method.append('SynMap:{}'.format(names[i]))
                else:
                    method.append('')
                meta[ui][gene]['Method']=method
    # ensure that extra exists
    try:
        os.mkdir('{}/extra'.format(resultdir))
    except:
        pass
    
    # Output mapping to file
    OUT=open('{}/extra/geneMapping_{}.txt'.format(resultdir,mode),'w')
    OUT.write('MeSHUI\tMeSHType\tEntrezID\tEntrezSymbol\tPubMedCoOccurences(Confidence)\tNM_Map\tSymbolMap\tSynonymMap\n')
    for ui in terms.keys():
        methods=[]
        symbols=[]
        weights=[]
        for gene in terms[ui]:
            symbols.append(meta[ui][gene]['Symbol'])
            weights.append(meta[ui][gene]['Weight'])
            methods.append(meta[ui][gene]['Method']) 
        methods=list(zip(*methods))
        if 'D' in ui:
            meshtype='MeSH'
        else:
            meshtype='SCR'
        OUT.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ui,meshtype,'|'.join(terms[ui]),'|'.join(symbols),'|'.join(weights),'|'.join(methods[0]),'|'.join(methods[1]),'|'.join(methods[2])))
    OUT.close() 
        
    # Get gene mapping of highest confidence
    d={}
    OUT=open('{}/extra/geneMapping_Detailed.txt'.format(resultdir),'w')
    OUT2=open('{}/geneMapping.txt'.format(resultdir),'w')
    OUT.write('MeSHUI\tEntrezID\tPubMedCoOccurences(Confidence)\tEntrezSymbol\tMeSHTpye\tMappingMethod\n')
    OUT2.write('MeSHUI\tEntrezID\tEntrezSymbol\n')
    for ui in terms.keys():

        methods=[]
        symbols=[]
        weights=[]
        genes=[]
        symbol=None
        for tmpgene in terms[ui]: # is this dict a problem, because it will override things
            symbols.append(meta[ui][tmpgene]['Symbol'])
            weights.append(meta[ui][tmpgene]['Weight'])
            methods.append(meta[ui][tmpgene]['Method'])
            genes.append(tmpgene)
            if methods[-1][0]!='': # if the gene is mapped via NM id, choose it
                symbol=symbols[-1]
                weight=weights[-1]
                method='|'.join(methods[-1])
                gene=tmpgene
         
        # select what gene to care about       
        if symbol==None:
            c=Counter(symbols)
            mostcommon=c.most_common() # List of tuples ( len 1 with only one most common)
            mostcommon=mostcommon[0] # subset to tuple
            symbol=mostcommon[0]
            ind=symbols.index(symbol)
            gene=genes[ind]
            method='|'.join(methods[ind])
            weight=weights[ind]
        d[ui]=gene
        if symbol!=None:
            pass
            #d[symbol]=gene
        if 'D' in ui:
            meshtype='MeSH'
        else:
            meshtype='SCR'

        OUT.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(ui,gene,weight,symbol,meshtype,method))
        OUT2.write('{}\t{}\t{}\n'.format(ui,gene,symbol))
    OUT.close()
    OUT2.close()

    return d
    
def getGene2Pubmed(datadir,humanspecific=1):
    mapping=defaultdict(list)
    ### Open the gene2pubmed file
    f=open('{}/mapping/gene/gene2pubmed'.format(datadir))
    for line in f:
        if line[0]=='#':
            continue
        line=line.strip().split('\t')
        if line[0]!='9606' and humanspecific:
            continue
        mapping[line[1]].append(line[2])
    return mapping

def getPmidUi(resultdir,datadir):
    ### Get pmid_ui Matrix 
    pmidUi=scipy.sparse.load_npz('{}/pmid_ui.npz'.format(resultdir))
    tmp=np.sum(pmidUi,axis=1)
    nonInds=tmp.nonzero()[0]
    
    # Get all PMIDs
    PMIDs=[line.strip().split('\t')[0] for line in open('{}/pmids.txt'.format(resultdir)).readlines()]
    PMIDs=np.array(PMIDs)
    # Subset PMIDs based on what is in pmidUi Matrix
    PMIDs=PMIDs[nonInds]
    pmidUi=pmidUi[nonInds]
    
    # Get indices of PMIDs
    enumpmid={}
    for i,val in enumerate(PMIDs):
        enumpmid[val]=i
    
    ### Get terms
    uis=np.array(open('{}/UIs.txt'.format(resultdir)).read().split('\n'))
    ### Get subset of relevant terms
    lines=open('{}/treemap.txt'.format(resultdir)).readlines()
    protein=[]
    meshNames={}
    nmids={}
    for line in lines:
        line=line.strip().split('\t')
        if line[2]=='1':
            protein.append(line[0])
        meshNames[line[0]]=line[3]
        if line[4]!='None':
            nmids[line[0]]=line[4]

    proteinInds=np.in1d(uis, protein) # get indices of the protein values
    # Subset
    pmidUi=pmidUi[:,proteinInds]
    uis=uis[proteinInds]
    return pmidUi,uis,enumpmid,meshNames,nmids

def getGeneInfo(datadir,humanspecific=1):
    ### Load Gene Info
    geneSyn={}
    geneSymbol={}
    f=open('{}/mapping/gene/gene_info'.format(datadir))
    for line in f:
        if line[0]=='#':
            continue
        line=line.strip().split('\t')
        if line[0]!='9606' and humanspecific:
            continue
        geneSymbol[line[1]]=line[2]
        geneSyn[line[1]]=line[13].split('|')
    return geneSymbol,geneSyn


def getGene2Refseq(datadir):
    ### Load Gene mapping to refseq
    geneRef=defaultdict(list)
    f=open('{}/mapping/gene/HGNC_2-8-18.txt'.format(datadir))
    for line in f:
        if 'HGNC' in line:
            continue
        line=line.strip('\n').split('\t')
        refseq=line[6].split(',')+line[8].split(',')+line[4].split(',')
        ui=line[7]
        for val in refseq:
            val=val.strip()
            if val!='':
                geneRef[ui].append(val)
    return geneRef
if __name__=='__main__':
    genMapping(datadir='../data',resultdir='../results')
    
    
# def getGene2Refseq(datadir,humanspecific=1):
#     #Load Gene mapping to refseq
#     geneRef=defaultdict(list)
#     f=open('{}/gene2refseq'.format(datadir))
#     for line in f:
#         if line[0]=='#':
#             continue
#         line=line.strip().split('\t')
#         if line[0]!='9606' and humanspecific:
#             continue
#         pn=line[6].split('.')[0]#protein nucleotide accession
#         pgi=line[7] #protein gi accession
#         gn=line[8].split('.')[0] #genomic nucleotide accession
#         ggi=line[9] #genomic gi accession
#         for val in [pn,pgi,gn,ggi]:
#             if val!='-':
#                 geneRef[line[1]].append(val)
#     return geneRef
