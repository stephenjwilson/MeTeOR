import re,urllib.request,urllib.parse,os
from http.client import IncompleteRead
from IPython import embed

def downloadUIs(fl='../data/EBVA.txt',searchterm='',step=5000,verbose=1):
    if searchterm=='':
        searchterm='"eukaryota"[MeSH Terms] OR "bacteria"[MeSH Terms] OR"viruses"[MeSH Terms] OR "archaea"[MeSH Terms]'
        
    efetch='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    esearch='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    db='pubmed'
    webpat=re.compile('<WebEnv>(.*)</WebEnv>')
    qpat=re.compile('<QueryKey>(.*)</QueryKey>')
    count=re.compile('<Count>(.*)</Count>')
    # Search query
    values = { 'db' : db,
               'term' : searchterm,
               'usehistory':'y'}
    valdata=urllib.parse.urlencode(values).encode('utf-8')
    req = urllib.request.Request(esearch, valdata)
    response = urllib.request.urlopen(req)
    data=response.read().decode('utf-8')
    # Get search history
    webenv=re.search(webpat,data).groups()[0]
    qkey=re.search(qpat,data).groups()[0]
    count=int(re.search(count,data).groups()[0])
    
    ids=''
    data=None
    for i in range(0,count,step):
        if i%100000==0 and verbose:
            print('{}/{}'.format(i,count))
        retmode='text'
        rettype='uilist'
        # Fetch 
        values = { 'db' : db,
                   'webenv' : webenv,
                   'query_key' : qkey,
                   'rettype':rettype,
                   'retmode':retmode,
                   'retstart':i,
                   'retmax':step}
        valdata=urllib.parse.urlencode(values).encode('utf-8')
        req=urllib.request.Request(efetch, valdata)
        response = urllib.request.urlopen(req)

        try:
            data=response.read().decode('utf-8')
        except IncompleteRead:
            try:
                data=response.read().decode('utf-8')
            except IncompleteRead:
                print('Failed to retreive')
        if 'Unable to obtain query' in data:
            try:
                data=response.read().decode('utf-8')
            except IncompleteRead:
                print('Failed to retreive')
            data=None

        ids+=data

    OUT=open(fl,'w')
    OUT.write(ids)
    OUT.close()
    return 1
def getData(storagedir='',fl='../data/EVBA.txt',searchterm='',step=5000,store='IDs'):
    '''Three primary functions: 
    1) store='IDs': gets the PMIDs for search and stores them (legacy)
    2) store='XML': gets the XML for the search and stores them
    3) store='None': gets XML and processes it'''
    options=['ids','xml','none']
    if store.lower() not in options:
        print('Bad option for store. Valid options include:{}.'.format(', '.join(options)))
        return 0
    if storagedir=='' and store.lower()=='xml':
        print('Please provide a storage directory as storagedir')
        return 0
    if searchterm=='':
        searchterm='"eukaryota"[MeSH Terms] OR "bacteria"[MeSH Terms] OR"viruses"[MeSH Terms] OR "archaea"[MeSH Terms]'
        
    efetch='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    epost='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi'
    esearch='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    db='pubmed'
    webpat=re.compile('<WebEnv>(.*)</WebEnv>')
    qpat=re.compile('<QueryKey>(.*)</QueryKey>')
    count=re.compile('<Count>(.*)</Count>')
    # Search query
    values = { 'db' : db,
               'term' : searchterm,
               'usehistory':'y'}
    valdata=urllib.parse.urlencode(values).encode('utf-8')
    req = urllib.request.Request(esearch, valdata)
    response = urllib.request.urlopen(req)
    data=response.read().decode('utf-8')
    # Get search history
    webenv=re.search(webpat,data).groups()[0]
    qkey=re.search(qpat,data).groups()[0]
    count=int(re.search(count,data).groups()[0])
    
    ids=''
    data=None
    for i in range(0,count,step):
        fl='{}/{}.xml'.format(storagedir,i)
        if os.path.exists(fl): # Checks if it has already been downloaded
            if os.stat(fl).st_size==0: # # If it exists, but is empty, remove it
                os.remove(fl)
            else:
                continue

        if i%100000==0:
            print(i)
        if store.lower()=='xml' or store.lower()=='none':
            retmode='XML'
            rettype='text'
        else:
            retmode='text'
            rettype='uilist'
        # Fetch 
        values = { 'db' : db,
                   'webenv' : webenv,
                   'query_key' : qkey,
                   'rettype':rettype,
                   'retmode':retmode,
                   'retstart':i,
                   'retmax':step}
        valdata=urllib.parse.urlencode(values).encode('utf-8')
        req=urllib.request.Request(efetch, valdata)
        response = urllib.request.urlopen(req)

        try:
            data=response.read().decode('utf-8')
        except IncompleteRead:
            try:
                data=response.read().decode('utf-8')
            except IncompleteRead:
                print('Failed to retreive')
        if 'Unable to obtain query' in data:
            try:
                data=response.read().decode('utf-8')
            except IncompleteRead:
                print('Failed to retreive')
            data=None
        if store.lower()=='ids':
            ids+=data
        elif store.lower()=='xml':
            if data==None:
                continue
            f=open('{}/{}.xml'.format(storagedir,i),'w')
            f.write(data)
            f.close()
        elif store.lower()=='none':
            pass
        else:
            print('This Error should not occur. store option is bad.')
            return 0
            
    if store.lower()=='ids': 
        OUT=open(fl,'w')
        OUT.write(ids)
        OUT.close()
    return 1


def downloadAll(storagedir):
    downloadUIs()
    getData(storagedir=storagedir,store='XML')


if __name__=='__main__':
    DIR='/home/stephen/ExtraDrive1/MeSH_V3'
    downloadAll(DIR)
