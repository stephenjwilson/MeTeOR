
import csv,itertools,obonet,os
from IPython import embed
import numpy as np
import scipy.sparse
### internal imports
import geneMapping

def main(datadir='../data',resultdir='../results',verbose=1):
    # MeTeOR and generate gene mapping
    mapMeteor(datadir=datadir, resultdir=resultdir,verbose=verbose)
 
    # Genes
    processBIOGRID(datadir=datadir, resultdir=resultdir,verbose=verbose)
    processMSigDB(datadir=datadir, resultdir=resultdir,verbose=verbose)
    processMSigDB(fl='c2.cp.v6.1.entrez.gmt',datadir=datadir, resultdir=resultdir,verbose=verbose)
    processEVEX(datadir=datadir, resultdir=resultdir,verbose=verbose) # Text-mining
    processSTRING(datadir=datadir, resultdir=resultdir,verbose=verbose)
     

    # Chemicals
    processBIOGRIDChem(datadir=datadir, resultdir=resultdir,verbose=verbose)
    processDGIdb(datadir=datadir, resultdir=resultdir,verbose=verbose)
    processSTITCH(datadir=datadir, resultdir=resultdir,verbose=verbose) # Text-mining
    # Diseases
    processHPO(datadir=datadir, resultdir=resultdir,verbose=verbose)
    processCTD(datadir=datadir, resultdir=resultdir,verbose=verbose)

    
### Process Sources ###
### Gene Sources ###  
def processBIOGRID(fl='BIOGRID-ALL-3.4.142.tab2.txt',datadir='../data',resultdir='../results',verbose=1):
    # Check to see if the file already exists
    networkdir='{}/networks'.format(resultdir)
        
    lines=open('{}/preprocess/BIOGRID/{}'.format(datadir,fl)).readlines()
    OUT=open('{}/BIOGRIDLow.txt'.format(networkdir),'w')
    OUT.write('#EntrezID\tEntrezID\tConfidence\tMetaData\n')
    OUT2=open('{}/BIOGRIDHigh.txt'.format(networkdir),'w')
    OUT2.write('#EntrezID\tEntrezID\tConfidence\tMetaData\n')
    experiments=[]
    for line in lines[1:]:
        line=line.strip().split('\t')
        ent1=line[1]# Entrez ID 1
        ent2=line[2]# Entrez ID 2
        experiment=line[11] # Experiment used
        experiments.append(experiment)
        #systemtype=line[12] # physical or genetic
        if 'Low Throughput' in line[17]:
            OUT.write('{}\t{}\t{}\t{}\n'.format(ent1,ent2,1,'Source::{}'.format(line[-1])))
        elif line[17]=='High Throughput':
            OUT2.write('{}\t{}\t{}\t{}\n'.format(ent1,ent2,1,'Source::{}'.format(line[-1])))
        else:
            embed()
            exit()
    OUT.close()
    OUT2.close()
    
    experiments=sorted(list(set(experiments)))
    experiments=[x.replace(' ','').replace('-','') for x in experiments]
    flhs=[open('{}/BIOGRID{}.txt'.format(networkdir,exp),'w') for exp in experiments]
    for flh in flhs:    
        flh.write('#EntrezID\tEntrezID\tConfidence\tMetaData\n')
    for line in lines[1:]:
        line=line.strip().split('\t')
        ent1=line[1]# Entrez ID 1
        ent2=line[2]# Entrez ID 2
        experiment=line[11].replace(' ','').replace('-','')# Experiment used
        ind=experiments.index(experiment)
        flhs[ind].write('{}\t{}\t{}\t{}\n'.format(ent1,ent2,1,'Source::{}'.format(line[-1])))
    
    if verbose:
        print('Processed BIOGRID')

def getGeneNames(fl='HGNC_2-8-18.txt',datadir='../data'):
    gtoName={}
    lines=open('{}/mapping/gene/{}'.format(datadir,fl)).readlines()
    for line in lines[1:]:
        line=line.strip('\n').split('\t')
        gtoName[line[5]]=line[1]
    return gtoName

def getMeshNames(resultdir='../results'):
    toName={}
    lines=open('{}/treemap.txt'.format(resultdir)).readlines()
    for line in lines:
        line=line.strip('\n').split('\t')
        toName[line[0]]=line[3]
    return toName              

def mapMeteor(fls=['coOccurrence.npz','coOccurrence_2014.npz'],datadir='../data',resultdir='../results',verbose=1):
    '''
    Maps the UIs for the MeTeOR Co-Occurrence Matrix
    '''
    uis=open('{}/UIs.txt'.format(resultdir)).read().split('\n')
    # Get mapping
    gmap=getGeneMap(datadir=datadir,resultdir=resultdir,verbose=verbose) # maps 
    cmap=getChemMap(datadir=datadir, resultdir=resultdir, verbose=verbose) # maps MeSH to Chem
    dmap=getMeshDieases(resultdir=resultdir, verbose=verbose)
    
    # get name mappings
    gNames=getGeneNames(datadir=datadir)
    names=getMeshNames(resultdir=resultdir)
    
    OUT=open('{}/UIsMapped.txt'.format(resultdir),'w')
    miss=0
    diseases=[]
    uisMapped=[]
    genes=[]
    chemicals=[]
    for i in range(0,len(uis)):
        ui=uis[i]
        if ui in dmap:
            OUT.write('{}\t{}\t{}\n'.format(ui,'Disease',names[ui]))
            diseases.append(i)
            uisMapped.append(ui)
        else:
            try:
                OUT.write('{}\t{}\t{}\n'.format(gmap[ui],'Gene',gNames[gmap[ui]]))
                genes.append(i)
                uisMapped.append(gmap[ui])
            except KeyError:
                # Not a disease or a gene
                try:
                    OUT.write('{}\t{}\t{}\n'.format(cmap[ui],'Chemical',names[ui]))
                    chemicals.append(i)
                    uisMapped.append(cmap[ui])
                except KeyError:
                    try:
                        OUT.write('{}\t{}\t{}\n'.format(ui,'Unknown',names[ui]))
                    except:
                        OUT.write('{}\t{}\t{}\n'.format(ui,'Unknown','BadID'))
                    miss+=1
                    uisMapped.append(ui)
    OUT.close()

    if verbose:
        print('Unlabeled {}/{} or {:.2f}%'.format(miss,len(uis),miss/len(uis)*100))
    # Push To flat file
    try:
        os.mkdir('{}/networks'.format(resultdir))
    except:
        pass
    for fl in fls:
        if '_' in fl:
            retro='_'+fl.split('_')[-1].replace('.npz','')
        else:
            retro=''
        mat=scipy.sparse.load_npz('{}/{}'.format(resultdir,fl)).tocsr()# load co-occurrence 
        uisMapped=np.array(uisMapped)
        headerMap={'chemical':'CID','gene':'EntrezID','disease':'MeSHID'} # Help make the header
        for mode1,mode2 in itertools.combinations_with_replacement([(diseases,'disease'),(genes,'gene'),(chemicals,'chemical')],2):
            OUT=open('{}/networks/MeTeOR{}{}{}.txt'.format(resultdir,mode1[1],mode2[1],retro),'w')
            # Make the Header
            header=[]
            header.append(headerMap[mode1[1]])
            header.append(headerMap[mode2[1]])
            header='#'+'\t'.join(header)+'\tConfidence\tMetaData\n'
            OUT.write(header)
            #subset mat
            sub=mat[mode1[0],:]
            sub=sub[:,mode2[0]]
            ui1=uisMapped[mode1[0]]
            ui2=uisMapped[mode2[0]]
            # Iterate over edges
            for i,j in np.transpose(sub.nonzero()):
                if i!=j: # don't write out i==j, as this is the degree, not an edge
                    OUT.write('%s\t%s\t%i\t%s\n' % (ui1[i],ui2[j],sub[i,j],'Mode:%s-%s' % (mode1[1],mode2[1]))) # faster than .format
            OUT.close()

def processEVEX(fl='EVEX_relations_9606.tab',datadir='../data',resultdir='../results',verbose=1):
    networkdir='{}/networks'.format(resultdir)
    
    lines=open('{}/preprocess/EVEX/{}'.format(datadir,fl)).readlines()
    OUT=open('{}/EVEX.txt'.format(networkdir),'w')
    OUT.write('#EntrezID\tEntrezID\tConfidence\tMetaData\n')

    for line in lines[1:]:
        line=line.strip().split('\t')
        ent1=line[1]# Entrez ID 1
        ent2=line[2]# Entrez ID 2
        conf=line[3]
        OUT.write('{}\t{}\t{}\t{}\n'.format(ent1,ent2,conf,'Source::{}'.format('EVEX')))
    OUT.close()

    if verbose:
        print('Processed EVEX')

        
def processMSigDB(fl='c2.all.v6.1.entrez.gmt',datadir='../data',resultdir='../results',verbose=1):
    from scipy import stats
    networkdir='{}/networks'.format(resultdir)
    lines=open('{}/preprocess/MsigDBCurated/{}'.format(datadir,fl)).readlines()
    OUT=open('{}/MSigDBCurated_{}.txt'.format(networkdir,fl.split('.')[1]),'w')
    OUT.write('#EntrezID\tEntrezID\tConfidence\tMetaData\n')
    OUT2=open('{}/MSigDBCurated_top_{}.txt'.format(networkdir,fl.split('.')[1]),'w')
    OUT2.write('#EntrezID\tEntrezID\tConfidence\tMetaData\n')
    genes=[]
    for line in lines:
        line=line.strip('\n').split('\t')
        genes+=line[2:]
    numpath=len(lines)
    genes=sorted(list(set(genes))) # eliminate duplicates, order
    mat=np.zeros((numpath,len(genes)))#lil_matrix((numpath,len(genes)))  # Can do in dense as it is only 4738 x 21091. Probably not good for GO
    for i in range(0,len(lines)):
        line=lines[i].strip('\n').split('\t')
        members=line[2:]
        subset=np.in1d(genes, members,assume_unique=True).astype(int)/len(members)
        mat[i,:]=subset
    mat[mat>0]=stats.zscore(mat[mat>0], axis=None) # Don't include the zeros in the normalization!
    if verbose:
        print('Retrieved MSigDB data. Calculating co-occurrence...')
    coOccurrence=np.dot(mat.T,mat)

    # Threshold matrix to 2
    top=np.percentile(coOccurrence,99.9)
    coOccurrence[coOccurrence < 2 ]=0

    links=np.transpose(np.nonzero(coOccurrence))

    for i,j in links:
        if i==j:
            continue
        val=coOccurrence[i,j]
        g1=genes[i]
        g2=genes[j]
        OUT.write('{}\t{}\t{}\t{}\n'.format(g1,g2,val,'Source::{}'.format('MsigDBCurated')))
        if val>top:
                OUT2.write('{}\t{}\t{}\t{}\n'.format(g1,g2,val,'Source::{}'.format('MsigDBCurated')))
        else: # below the threshold
            pass
    OUT.close()
    OUT2.close()
    if verbose:
        print('Processed MSigDBCurated')
                
def processSTRING(fl='9606.protein.links.detailed.v10.5.txt',resultdir='../results',datadir='../data',verbose=1):
    networkdir='{}/networks'.format(resultdir)
    lines=open('{}/preprocess/STRING/{}'.format(datadir,fl)).readlines()
    
    # Get GeneMapping
    s9,s10=getENSEMBL(datadir=datadir,resultdir=resultdir,verbose=verbose)
    
    miss=0
    s9.update(s10) # use both mappings
    s10=s9
    header=lines[0].strip('\n').split(' ')
    for col in header[2:]:
        ind=header.index(col)
        OUT=open('{}/STRING10_{}.txt'.format(networkdir,col),'w')
        OUT.write('#EntrezID\tEntrezID\tConfidence\tMetaData\n')
        for line in lines[1:]:
            line=line.strip('\n').split(' ')
            ent1=line[0]
            ent2=line[1]
            
            score=line[ind]
           
            try:
                ent1=s10[ent1]
            except KeyError:
                miss+=1
                continue
            try:
                ent2=s10[ent2]
            except KeyError:
                miss+=1
                continue
            
            if int(score)>0:
                OUT.write('{}\t{}\t{}\t{}\n'.format(ent1,ent2,score,'Source::{}'.format(col)))
        OUT.close()
    if verbose:
        n=len(lines)-1
        print('Processed STRING. Missed {}/{} or {:.2f}%.'.format(miss,n,miss/n*100))

    
### Chemical Sources ###
def processBIOGRIDChem(fl='BIOGRID-CHEMICALS-3.4.142.chemtab.txt',datadir='../data',resultdir='../results',verbose=1):
    networkdir='{}/networks'.format(resultdir)
    lines=open('{}/preprocess/BIOGRID/{}'.format(datadir,fl)).readlines()
    # Get ChemicalMapping
    mapping=getChemMap(datadir=datadir, resultdir=resultdir, verbose=verbose)
    OUT=open('{}/BIOGRIDchem.txt'.format(networkdir),'w')
    OUT.write('#EntrezID\tCID\tConfidence\tMetaData\n')
    miss=0
    for line in lines[1:]:
        line=line.strip().split('\t')
        ent1=line[2] # Entrez ID
        ent2=line[18] # DrugBank ID
        # Map Drugbank to CID
        try:
            ent2=mapping[ent2]
        except KeyError:
            miss+=1
            continue
        OUT.write('{}\t{}\t{}\t{}\n'.format(ent1,ent2,1,'Source::{}'.format(line[-1])))
    OUT.close()
    if verbose:
        n=len(lines)-1
        print('Processed BIOGRID Chem. Missed {}/{} or {:.2f}%'.format(miss,n,miss/n*100))

def processDGIdb(fl='interactions.tsv',datadir='../data',resultdir='../results',verbose=1):
    networkdir='{}/networks'.format(resultdir)
    lines=open('{}/preprocess/DGIdb/{}'.format(datadir,fl)).readlines()
    # Get ChemicalMapping
    mapping=getChemMap(datadir=datadir, resultdir=resultdir, verbose=verbose)
    fls={}
    miss=0
    nomap=0
    OUT=open('{}/DGIdb.txt'.format(networkdir),'w')
    OUT.write('#EntrezID\tCID\tConfidence\tMetaData\n')
    for line in lines[1:]:
        line=line.strip('\n').split('\t')
        ent1=line[2] # Entrez ID
        ent2=line[8] # ChEMBL ID
        source=line[3] # Source compiled from DGIdb
        if ent2=='':
            nomap+=1
        # Map ChEMBL to CID
        try:
            ent2=mapping[ent2]
        except KeyError:
            miss+=1
            continue
        try:
            fls[source]
        except:
            fls[source]=open('{}/DGIdb_{}.txt'.format(networkdir,source),'w')
            fls[source].write('#EntrezID\tCID\tConfidence\tMetaData\n')
        fls[source].write('{}\t{}\t{}\t{}\n'.format(ent1,ent2,1,'Source::{}'.format(source)))
        OUT.write('{}\t{}\t{}\t{}\n'.format(ent1,ent2,1,'Source::{}'.format(source)))
    for source in fls.keys():
        fls[source].close()
    if verbose:
        n=len(lines)-1
        print('Processed DGIdb. Missed {}/{} or {:.2f}%. No ChEMBL id for:{} or {:.2f}. Adjusted miss: {:.2f}'.format(miss,n,miss/n*100,nomap,nomap/n*100,(miss-nomap)/n*100))

def processSTITCH(fl='9606.protein_chemical.links.detailed.v5.0.tsv',resultdir='../results',datadir='../data',verbose=1):
    networkdir='{}/networks'.format(resultdir)
    lines=open('{}/preprocess/STITCH/{}'.format(datadir,fl)).readlines()
    
    # Get GeneMapping
    _,s10=getENSEMBL(datadir=datadir,resultdir=resultdir,verbose=verbose)
    
    miss=0
    header=lines[0].strip('\n').split('\t')
    for col in header[2:]:
        ind=header.index(col)
        OUT=open('{}/STITCH5_{}.txt'.format(networkdir,col),'w')
        OUT.write('#EntrezID\tCID\tConfidence\tMetaData\n')
        for line in lines[1:]:
            line=line.strip('\n').split('\t')
            cid=line[0].replace('CIDm','').replace('CIDs','') # m is a flat compound, s is stereo-specific
            try:
                int(cid) # if not, it is a SID. Ignore
            except:
                continue
            ensembl=line[1]
            score=line[ind]
            try:
                gene=s10[ensembl]
            except KeyError:
                miss+=1
                continue
            if int(score)>0:
                OUT.write('{}\t{}\t{}\t{}\n'.format(gene,int(cid),score,'Source::{}'.format(col)))
        OUT.close()
    if verbose:
        n=len(lines)-1
        print('Processed STITCH. Missed {}/{} or {:.2f}%.'.format(miss,n,miss/n*100))
            
### Disease Sources ###
def processDisgenet(fl1='curated_gene_disease_associations.tsv',fl2='befree_gene_disease_associations.tsv',fl3='befree_gene_disease_associations_V4.tsv',resultdir='../results',datadir='../data',verbose=1):
    networkdir='{}/networks'.format(resultdir)

    ### Get curated associations
    f=open('{}/preprocess/Disgenet/{}'.format(datadir,fl1))
    lines=f.readlines()
    f.close()

    data={}
    d=getDisgenetMapping(datadir=datadir,resultdir=resultdir,verbose=verbose)
    miss=0

    for line in lines[1:]:
        line=line.strip('\n').split('\t')
        gene=line[0] # Entrez id
        disease=line[2] # UMLS id
        score=line[4] # Confidence score
        try:
            disease=d[disease]
        except KeyError:
            miss+=1
            continue
        sources=line[7].split(';') # sources include many things
        for source in sources:
            if source in data.keys():
                data[source].append('{}\t{}\t{}\t{}\n'.format(gene,disease,score,'Source::{}'.format(line[-1])))
            else:
                data[source]=['{}\t{}\t{}\t{}\n'.format(gene,disease,score,'Source::{}'.format(line[-1]))]
    for source in data.keys():
        if 1:#source in ['HPO','CTD_human']: # only output the sources we use
            OUT=open('{}/{}_fromDisgenet.txt'.format(networkdir,source),'w')
            OUT.write('#EntrezID\tMeSHID\tConfidence\tMetaData\n')        
            for line in data[source]:
                OUT.write(line)
            OUT.close()
    if verbose:
        n=len(lines)-1
        print('Processed DisGeNET. Missed {}/{} or {:.2f}%.'.format(miss,n,miss/n*100))
    
    # do BeFree as well
    f=open('{}/preprocess/Disgenet/{}'.format(datadir,fl2),encoding='iso-8859-1')
    lines=f.readlines()
    f.close()
    OUT=open('{}/BeFree.txt'.format(networkdir),'w')
    OUT.write('#EntrezID\tMeSHID\tConfidence\tMetaData\n')
    miss=0
    for line in lines[1:]:
        line=line.strip('\n').split('\t')
        gene=line[0] # Entrez id
        disease=line[2] # UMLS id
        score=line[4] # Confidence score
        numofpmids=line[5] # can't use their confidence score, as this factors in other sources 
        try:
            disease=d[disease]
        except KeyError:
            miss+=1
            continue
        OUT.write('{}\t{}\t{}\t{}\n'.format(gene,disease,numofpmids,'Source::{}'.format(line[-1])))
        
    if verbose:
        n=len(lines)-1
        print('Processed BeFree. Missed {}/{} or {:.2f}%.'.format(miss,n,miss/n*100))
    
    # do BeFreeV4 as well
    f=open('{}/preprocess/Disgenet/{}'.format(datadir,fl3),encoding='iso-8859-1')
    lines=f.readlines()
    f.close()
    OUT=open('{}/BeFreeV4.txt'.format(networkdir),'w')
    OUT.write('#EntrezID\tMeSHID\tConfidence\tMetaData\n')
    miss=0
    comments=0
    for line in lines:
        if line[0]=='#' or 'umls' not in line[:10]: 
            comments+=1
            continue
        
        line=line.strip('\n').split('\t')
        gene=line[1] # Entrez id
        disease=line[0].replace('umls:','') # UMLS id
        score=line[2] # Confidence score
        numpmids=line[7] # can't use their confidence score, as this factors in other sources 
        try:
            disease=d[disease]
        except KeyError:
            miss+=1
            continue
        OUT.write('{}\t{}\t{}\t{}\n'.format(gene,disease,numpmids,'Source::{}'.format(line[6])))
        
    if verbose:
        n=len(lines)-comments
        print('Processed BeFreeV4. Missed {}/{} or {:.2f}%.'.format(miss,n,miss/n*100))
        
def processHPO(fl='phenotype_annotation.tab',resultdir='../results',datadir='../data',verbose=1):
    networkdir='{}/networks'.format(resultdir)
    f=open('{}/preprocess/HPO/{}'.format(datadir,fl))
    lines=f.readlines()
    f.close()
    d=getHPOtoMeSH(datadir=datadir,resultdir=resultdir,verbose=verbose)
    OUT=open('{}/HPO.txt'.format(networkdir),'w')
    OUT.write('#EntrezID\tMeSHID\tConfidence\tMetaData\n')
    miss=0
    total=0
    for line in lines:
        line=line.strip('\n').split('\t')
        gene=line[1] # Entrez id
        disease=line[4] # HPO id
        # Filter to not IEA
        if line[6]=='IEA':
            continue
        total+=1
        try:
            disease=d[disease]
        except KeyError:
            miss+=1
            continue
        OUT.write('{}\t{}\t{}\t{}\n'.format(gene,disease,1,'Source::{}'.format(line[6])))
    OUT.close()
    if verbose:
        print('Processed HPO. Missed {}/{} or {:.2f}%.'.format(miss,total,miss/total*100))
        
def processCTD(fl='CTD_genes_diseases.tsv',fl2='CTD_chem_gene_ixns.tsv',resultdir='../results',datadir='../data',verbose=1):
    networkdir='{}/networks'.format(resultdir)
    ### Gene Disease
    f=open('{}/preprocess/CTD/{}'.format(datadir,fl))
    OUT=open('{}/CTD_GD.txt'.format(networkdir),'w')
    OUT.write('#EntrezID\tMeSHID\tConfidence\tMetaData\n')
    miss=0
    total=0
    fls={}
    for line in f:
        if line[0]=='#':
            continue
        
        line=line.strip('\n').split('\t')
        

        gene=line[1] # Entrez id
        disease=line[3] # Disease id, sometimes MeSH
        directEvs=line[4] # DirectEvidence code
        if directEvs=='': # means indirect
            continue
        total+=1
            
        if 'MESH' not in disease:
            miss+=1
            continue
        disease=disease.replace('MESH:','')
        # Send out the evidence specific type
        for directEv in directEvs.split('|'):
            directEv=directEv.replace('/','_')    
            try:
                fls[directEv]
            except:
                fls[directEv]=open('{}/CTD_GD_{}.txt'.format(networkdir,directEv),'w')
                fls[directEv].write('#EntrezID\tMeSHID\tConfidence\tMetaData\n')
            fls[directEv].write('{}\t{}\t{}\t{}\n'.format(gene,disease,1,'Source::{}'.format(line[-1])))
        OUT.write('{}\t{}\t{}\t{}\n'.format(gene,disease,1,'Source::{}'.format(line[-1])))
        
    OUT.close()
    f.close()
    if verbose:
        print('Processed CTD Gene-Disease. Missed {}/{} or {:.2f}%.'.format(miss,total,miss/total*100))
        
    ### Chem Gene
    f=open('{}/preprocess/CTD/{}'.format(datadir,fl2))
    cmap=getChemMap(datadir=datadir, resultdir=resultdir, verbose=verbose) # maps MeSH to Chem
    fls={}
    OUT=open('{}/CTD_CG.txt'.format(networkdir),'w')
    OUT.write('#EntrezID\tCID\tConfidence\tMetaData\n') 
    miss=0
    total=0
    for line in f:
        if line[0]=='#':
            continue
        
        line=line.strip('\n').split('\t')

        total+=1
        gene=line[4] # Entrez id
        chem=line[1] # MeSH id for chem
        actions=line[9] # Interaction Actions
        geneForms=line[5] # the gene form
        conf=len(line[10].split('|')) # the number of pmids
        if chem not in cmap:
            miss+=1
            continue
        chem=cmap[chem]
        for action in actions.split('|'):
            action=action.replace(' ' ,'').replace("'",'').replace('-','_').split('^')
            try:
                fls[action[1]]
            except:
                fls[action[1]]=open('{}/CTD_CG_{}.txt'.format(networkdir,action[1]),'w')
                fls[action[1]].write('#EntrezID\tCID\tConfidence\tMetaData\n')               
            fls[action[1]].write('{}\t{}\t{}\t{}\n'.format(gene,chem,conf,'Source::{}|Effect::{}'.format(line[-1],action[0])))
        for geneForm in geneForms.split('|'):
            geneForm=geneForm.replace(' ' ,'').replace("'",'').replace('-','_')
            try:
                fls[geneForm]
            except:
                fls[geneForm]=open('{}/CTD_CG_{}.txt'.format(networkdir,geneForm),'w')
                fls[geneForm].write('#EntrezID\tCID\tConfidence\tMetaData\n')               
            fls[geneForm].write('{}\t{}\t{}\t{}\n'.format(gene,chem,conf,'Source::{}'.format(line[-1])))            
        OUT.write('{}\t{}\t{}\t{}\n'.format(gene,chem,conf,'Source::{}'.format(line[-1])))
          
    for key in fls.keys():
        fls[key].close()
    OUT.close()
    f.close()
    if verbose:
        print('Processed CTD Gene_chem. Missed {}/{} or {:.2f}%.'.format(miss,total,miss/total*100))
    
    

### Make Mapping ###
### Gene Mapping ###
def getGeneMap(datadir='../data',resultdir='../results',verbose=1):
    alltoEntrez={}
    # Generate Base mapping with MeSH and SCR information
    d=geneMapping.genMapping(datadir=datadir,resultdir=resultdir, mode='any') # If this is taking too long, you can switch to mode='max' for mapping only the most likely terms.
    alltoEntrez.update(d)
    # Add on ENSEMBL Data
    s9,s10=getENSEMBL(datadir=datadir,resultdir=resultdir,verbose=verbose)
    # Build Mapping
    alltoEntrez.update(s9)
    alltoEntrez.update(s10)
    if verbose:
        print('Loaded Gene Mapping')
    return alltoEntrez

def getENSEMBL(s9='entrez_gene_id.vs.string.v9.05.28122012.txt',s10='entrez_gene_id.vs.string.v10.28042015.tsv',datadir='../data',resultdir='../results',verbose=1):
    def processSTRING(fl,resultdir):
        d={}
        f=open('{}/mapping/gene/{}'.format(datadir,fl))
        for line in f:
            if line[0]=='#':
                continue
            line=line.strip('\n').split('\t')
            ent1=line[0] # Entrez ID
            ent2=line[1] # STRING ID
            d[ent1]=ent2 # Maps Entrez to STRING
        return d
    s9=processSTRING(s9,resultdir)
    s10=processSTRING(s10,resultdir)
    
    # Update geneMapping
    lines=open('{}/geneMapping.txt'.format(resultdir)).readlines()
    if len(lines[0].split('\t'))<4: # only write it out if it hasn't been done before
        OUT=open('{}/geneMapping.txt'.format(resultdir),'w')
        OUT.write(lines[0].strip('\n')+'\tString9\tString10\n') # update header
        for line in lines[1:]: # Skip header
            line=line.strip('\n').split('\t')
            for d in [s9,s10]:
                try:
                    line.append(d[line[1]])
                except KeyError:
                    line.append('')
            OUT.write('{}\n'.format('\t'.join(line)))
                                
    #reverse dictionary to map string to entrez
    s9={v:k for k,v in s9.items()}
    s10={v:k for k,v in s10.items()}
    return s9,s10

### Chemical Mapping ###  
def getChemMap(datadir='../data',resultdir='../results',verbose=1):
    links=mapAlltoCid(datadir=datadir, resultdir=resultdir)
    alltoCID={}
    for link in links:
        for part in [link[0]]+link[1:]:
            if part!='' and link[1]!='':
                alltoCID[part]=link[1]
    if verbose:
        print('Finished Chemical Mapping...')
    return alltoCID

def mapAlltoCid(fl='drug links.csv',datadir='../data',resultdir='../results',verbose=1):
    # Attempt to Load from file
    lines=[]
    data=None
    try:
        data=open('{}/chemMapping.txt'.format(resultdir))
    except:
        pass
    if data!=None:
        for line in data:
            lines.append(line.strip('\n').split('\t'))
        if verbose:
            print('Loaded Chemical Mapping...')
        return lines
    
    # Calculate if mapping doesn't exist
    if verbose:
        print('Calculating Chemical Mapping...')
    # Load DrugBank mapping
    lines=open('{}/mapping/chemical/{}'.format(datadir,fl)).readlines()
    links=[]
    
    SIDsLackingCID=[]
    for line in  csv.reader(lines[1:], quotechar='"', delimiter=',',
                     quoting=csv.QUOTE_ALL, skipinitialspace=True):
        DBID=line[0] # Drugbank ID
        TTD=line[-1] # TTD ID
        cid=line[6] # PubChem Compound ID
        sid=line[7] # PubChem Substance ID
        if cid=='':
            SIDsLackingCID.append(sid)
        syns=[cid,sid,DBID,TTD]
        links.append(syns)
    if verbose:
        print('Loaded DrugBank Mapping...')
    # Retrieve CID-SID and CID-MeSH mappings(from NCBI)
    CIDtoSIDmap=CIDtoSID(SIDsLackingCID,datadir=datadir,verbose=verbose)
    MESHtoCIDmap=MESHtoCID(datadir=datadir, resultdir=resultdir,verbose=verbose)
    chembltoCIDmap=chembltoCID(datadir=datadir,verbose=verbose)
    
    used=[]
    for i in range(0,len(links)):
        # Replace any CID that doesn't exist based on SID
        if links[i][0]=='':
            try:
                links[i][0]=CIDtoSIDmap[links[i][1]]
            except KeyError:
                continue
        # Add MeSH
        try:
            links[i].insert(0,MESHtoCIDmap[links[i][0]])
            used.append(links[i][0])
        except KeyError:
            links[i].insert(0,'')
    for key in MESHtoCIDmap.keys():
        if key not in used:
            links.append([MESHtoCIDmap[key],key,'','',''])
    # Add in CHEMBL
    for i in range(0,len(links)):
        cid=links[i][1]
        try:
            links[i].append(chembltoCIDmap[cid])
        except KeyError:
            links[i].append('')
    # Write out mapping to fl
    OUT=open('{}/chemMapping.txt'.format(resultdir),'w')
    OUT.write('MeSHID\tCID\tSID\tDBID\tTTD\tChEMBL\n')
    for link in links:
        OUT.write('\t'.join(link)+'\n')
    OUT.close()
    return links

def chembltoCID(fl='src1src22.txt',datadir='../data',verbose=1):
    '''Maps the ChEMBL IDs to CID''' 
    d={}
    f=open('{}/mapping/chemical/{}'.format(datadir,fl))
    for line in f:
        line=line.strip('\n').split('\t')
        d[line[1]]=line[0] # maps CID to Chembl
    if verbose:
        print('Loaded ChEMBL-CID...')
    return d

def CIDtoSID(sids,fl='CID-SID',datadir='../data',verbose=1):
    '''Determine mappings of CID to SID for CIDs that are missing'''
    d={sid:'' for sid in sids}
    
    # Try to load from file
    data=None
    try:
        data=open('{}/mapping/chemical/CIDtoSIDsmall'.format(datadir))
    except:
        pass
    if data!=None:
        for line in data:
            line=line.strip('\n').split('\t')
            d[line[0]]=line[1]
        if verbose:
            print('Loaded CID-SID Mapping...')
        return d
    if verbose:
        print('Calculating CID-SID Mapping...')
    f=open('{}/mapping/chemical/{}'.format(datadir,fl))
    c=0
    for line in f:
        line=line.strip().split('\t')
        # Focus only on sids in file
        try:
            d[line[1]] #checks if the sid is present
            d[line[1]]=line[0] # maps SIDs to CIDs
            c+=1
        except:
            pass
        if c==len(sids):
            break
    OUT=open('{}/mapping/chemical/CIDtoSIDsmall'.format(datadir),'w')
    for key in d.keys():
        OUT.write('{}\t{}\n'.format(key,d[key]))
    OUT.close()
    if verbose:
        print('Loaded CID-SID...')
    return d

def MESHtoCID(fl='CID-MeSHID',datadir='../data',resultdir='../results',verbose=1):
    '''Processes the CID-MeSH file from NCBI into a mapping from MeSH to CID'''
    d={}
    # try to open the CID-MeSHID file (processed form of CID-MeSH which uses MeSH ID's not names)
    try:
        f=open('{}/mapping/chemical/{}'.format(datadir,fl))
    except:
        processMESHtoCID(datadir=datadir,resultdir=resultdir,verbose=verbose)
        f=open('{}/mapping/chemical/{}'.format(datadir,fl))
    for line in f:
        line=line.strip().split('\t')
        d[line[0]]=line[1] # maps CIDs to MeSH
    if verbose:
        print('Loaded CID-MeSHID...')
    return d

def processMESHtoCID(fl='CID-MeSH',resultdir='../results',datadir='../data',verbose=1):
    '''Processes the CID-MeSH file from NCBI into a CID-MeSHID file which uses MeSH IDs not names)'''
    # Open TreeMap, which contains IDs to Names
    f=open('{}/treemap.txt'.format(resultdir))
    nametoid={}
    for line in f:
        line=line.strip('\n').split('\t')
        nametoid[line[3]]=line[0]
    # Open CID-MeSH
    f=open('{}/mapping/chemical/{}'.format(datadir,fl))
    miss=0
    count=0
    OUT=open('{}/mapping/chemical/{}ID'.format(datadir,fl),'w')
    for line in f:
        count+=1
        line=line.strip('\n').split('\t')
        try:
            ui=nametoid[line[1]]
        except KeyError:
            miss+=1
            continue
        OUT.write('{}\t{}\n'.format(line[0],ui))
    if verbose:
        print('Processed CID-MeSH...{} / {} missed chemicals'.format(miss,count))
    return
### Disease Mappings ###
def getDisgenetMapping(fl='disease_mappings.tsv',datadir='../data',resultdir='../results',verbose=1):
    '''Processes the DisGeNET mapping file into a mapping from UMLS to MeSH'''
    d={}
    f=open('{}/mapping/disease/{}'.format(datadir,fl))
    for line in f:
        line=line.strip('\n').split('\t')
        if line[2]=='MSH':
            d[line[0]]=line[3] # maps UMLS to MeSH
    if verbose:
        print('Loaded UMLS to MeSH...')
    return d
def getHPOtoMeSH(fl='hp.obo',datadir='../data',resultdir='../results',verbose=1):
    '''Processes the HPO obo file into a mapping from HPO to MeSH'''

    f=open('{}/mapping/disease/{}'.format(datadir,fl))
    graph=obonet.read_obo(f)
    d={}
    for id_, data in graph.nodes(data=True):
        if 'xref' in data.keys():
            xrefs=data['xref']
            for xref in xrefs:
                if 'MSH' in xref:
                    d[id_]=xref.replace('MSH:','')
    if verbose:
        print('Loaded HPO to MeSH...')
    return d
def getMeshDieases(fl='treemap.txt',resultdir='../results',verbose=1):
    f=open('{}/{}'.format(resultdir,fl))
    diseaseUIs=[]
    for line in f:
        line=line.strip('\n').split('\t')
        ui=line[0]
        if line[1]=='':
            continue
        treenums=[True for num in line[1].split(',') if num[0]=='C']
        if sum(treenums)>0:
            diseaseUIs.append(ui)
    if verbose:
        print('Found {} Disease MeSH terms'.format(len(diseaseUIs)))
    return diseaseUIs

if __name__=='__main__':
    main(resultdir='../results/')
