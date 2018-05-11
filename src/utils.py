from lxml import etree  # @UnresolvedImport
#import xml.etree.cElementTree as etree
from IPython import embed
import re,os
from collections import defaultdict

def setup(resultdir='../results',datadir='../data',IDsfl='',verbose=1):
    '''
    Loads PMIDs if nee
    '''
    if IDsfl=='' and os.path.exists('{}/EBVA.txt'.format(datadir)):
        IDsfl='{}/EBVA.txt'.format(datadir)
        
    # extract all possible uis and if they are in the protein section
    uis,_=extractUIs(resultdir=resultdir,datadir=datadir)
    if IDsfl!='':
        # Load PMIDs
        f=open(IDsfl,'r')
        storedIDs=f.read()
        f.close()
        storedIDs=storedIDs.strip().split('\n')
        if verbose:
            print('Loaded PMIDs...')
    else:
        storedIDs=[]
        
    ### Make a dictionary of PMIDs to their gene ids
    examplegenelink=defaultdict(list)
    lines=open('{}/mapping/gene/gene2pubmed'.format(datadir)).readlines()
    for line in lines:
        line=line.strip().split('\t')
        if '9606'==line[0]: # remove if not human only
            examplegenelink[line[2]].append(line[1])
    if verbose:
        print('Loaded gene2pubmed...')
    # Get indicies
    enumpmid={}
    for i,val in enumerate(storedIDs):
        enumpmid[val]=i
    
    enumui={}
    for i,val in enumerate(uis):
        enumui[val]=i
 
    return examplegenelink,enumpmid,enumui,uis

def parse(infl,enumpmid,enumui,uis,meta,recover=False):

    pmids=[]
    meshids=[]
    miss=0
    indicies=[]
    addeduis=[]
        
    def getTerms(elem):
        '''Extracts MeSH and SCR terms'''
        #mesh
        mesh=[]
        meshel=elem.find('MedlineCitation').find('MeshHeadingList')
        if meshel!=None: # if it exists
            meshel=meshel.findall('MeshHeading')
            for m in meshel:
                mesh.append(m.find('DescriptorName').attrib['UI'])

        #chem
        chemel=elem.find('MedlineCitation').find('ChemicalList')
        if chemel!=None: # if it exists
            chemel=chemel.findall('Chemical')
            for m in chemel:
                mesh.append(m.find('NameOfSubstance').attrib['UI'])

        return mesh
    def getYears(elem):
        dates=[]
        try:
            dates=elem.find('PubmedData').find('History').findall('PubMedPubDate')
        except:
            try:
                year=elem.find('PubmedData').find('MedlineCitation').find('Article').find('ArticleDate').find('Year').text
                if int(year)<1500: # make sure it is accurate
                    return []
            except:
                return []
            return [year]
        years=[]
        for child in dates:
            years.append(child.find('Year').text)
        return years
        
    context = etree.iterparse(infl, events=('end',), tag='PubmedArticle',encoding='utf-8',recover=recover)
    currentPMID=len(enumpmid.keys())
    for _,elem in context:
        pmid= elem.find('MedlineCitation').find('PMID').text # First PMID is the right one
        if pmid==None or pmid=='':
            miss+=1
            continue
        try:
            pmidInd=enumpmid[pmid]
        except KeyError:
            pmidInd=currentPMID
            currentPMID+=1
            enumpmid[pmid]=pmidInd
        
        mesh = getTerms(elem)
        for ui in mesh:
            try:
                indicies.append((pmidInd,enumui[ui]))
            except:
                if 'D' in ui or 'C' in ui:
                    addeduis.append(ui)
                    enumui[ui]=len(uis)
                    uis.append(ui)
        
        years=getYears(elem)
        if years!=[]:
            year=min(years)
            meta[pmid]['Year']=year
            

        if mesh!=[]:
            meshids.append(pmid)
        pmids.append(pmid)
        elem.clear()
    #print('Added {} uis'.format(len(addeduis)))
     
    return indicies,pmids,meshids,miss        


def extractUIs(dbin='d2018.bin',cbin='c2018.bin',datadir='../data/',resultdir='../results'):
    '''Extracts UIs/treenum from MeSH descriptor ASCII file and UIs from supplemental chemicals'''
    treemap=defaultdict(list)
    uis=[]
    protein={}
    names={}
    nmids={}
    for fl in [dbin,cbin]:
        lines=open('{}/{}'.format(datadir,fl)).readlines()
        treenums=[]
        nmid=None
        if 'c20' in fl: # if supplemental records
            go=0
        else:
            go=1
        for line in lines:
            if go:
                if line[0:2]=='MN':
                    treenums.append(line.split('=')[1].strip())
                if line[0:3]=='MH ':
                    name=line.split('=')[1].strip()
            else:
                if line[0:3]=='NM ':
                    name=line.split('=')[1].strip()
                if line[0:3]=='NO ': # this presents the NM_ID mapping
                    nmid=extractGeneIDs(line)
                    if nmid=='':
                        nmid=None
            if line[0:2]=='UI':
                ui=line.split('=')[1].strip()
                if 'D12' in ''.join(treenums) or 'protein' in name.lower():
                    protein[ui]=1
                else:
                    protein[ui]=0
                if go:
                    treemap[ui]+=treenums
                    treenums=[]
                names[ui]=name
                nmids[ui]=nmid
                nmid=None
                uis.append(ui)

    # write out UIs
    OUT=open('{}/UIs.txt'.format(resultdir),'w')
    OUT.write('\n'.join(uis))
    OUT.close()
    
    # write out treemap
    OUT=open('{}/treemap.txt'.format(resultdir),'w')
    for ui in uis:
        if ui in treemap:
            treenum=','.join(treemap[ui])
        else:
            treenum=''
        OUT.write('{}\t{}\t{}\t{}\t{}\n'.format(ui,treenum,protein[ui],names[ui],str(nmids[ui])))
    OUT.close()
    
    return uis,protein

def extractGeneIDs(line,nmid=None):
    '''Attempts to extract a refseq id from a line of the c201*.bin file in the NO or Notes section of an entry'''
    pat = re.compile('.*(.._[0-9]+)')
    try:
        # Specific Match
        nmid = re.match(pat,
                        line.replace(' ', '').replace('-', '_')).groups()[0]
        return nmid
    except AttributeError:
        pass
    pat = re.compile('.*(Gene?Bank|RefSeq);?([0-9A-z]+)')
    try:
        nmid = re.match(pat, line.replace(' ', '').replace('accession', '').replace(
            'no.', '').replace('Ac', '').replace('#', '').replace(':', '')).groups()[1]
        if len(nmid) < 5:
            # print line[1], 'MATCH FAIL'
            # print nmid
            # NMID=''
            pass
        elif len(nmid) > 9:
            nmid = nmid.split('and')[0].strip().split('through')[0].strip().split('homologous')[0].strip().split('donotconfus')[0].strip().split('Donotconfu')[0].strip().split('seealso')[
                0].strip().replace('mRNA', '').split('protein')[0].strip().split('in')[0].strip().split('to')[0].strip().split('isoform')[0].strip().split('-')[0].strip().split('for')[0].strip()
            ind = re.match('.+([0-9])[^0-9]*$', nmid)
            nmid = nmid[:ind.start(1) + 1]
            if len(nmid) > 10:
                # print line[1]
                # print '    '+nmid
                # There is only one that fits it, and that is
                # the id.  AABR03074672
                pass
    except:
        # All cases of failing, the record is a gene
        # nmid=line[1]
        # print line[1], 'FAIL'
        pass
    # print line[1], 'total fail'
    return nmid

def mapName(fl,resultdir='../results',cols=[0,1],hasPrefix=1,delim='\t'):
    '''
    Uses a mapping file to map ids to names
    '''
    mapfl=open('{}/UIsMapped.txt'.format(resultdir)).readlines()
    mapping={}
    prefix={'Gene':'G.','Chemical':'C.','Disease':'D.','Unknown':'U.'}
    for line in mapfl:
        line=line.strip('\n').split('\t')
        if 'badid' not in line[2].lower():
            if hasPrefix:
                mapping[prefix[line[1]]+line[0]]=line[2]
            else:
                mapping[line[0]]=line[2]
    f=open(fl)
    OUT=open('.'.join(fl.split('.')[:-1])+'_Names.txt','w')
    for line in f:
        line=line.strip('\n').split(delim)
        for col in cols:
            try:
                line[col]=mapping[line[col]]
            except KeyError:
                pass
        line='\t'.join(line)+'\n'
        OUT.write(line)
    f.close()
    OUT.close()
    


