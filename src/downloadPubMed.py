import re,urllib.request,urllib.parse,os, urllib.error
from http.client import IncompleteRead
from IPython import embed
import gzip

def getData(storagedir='',fl='../data/EBVA.txt',searchterm='',step=5000):
    '''Gets the XML for the search and stores them'''
    if storagedir=='':
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
    
    retmode='XML'
    rettype='text'
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

        if data==None:
            continue
        ids+=data
        f=open('{}/{}.xml'.format(storagedir,i),'w')
        f.write(data)
        f.close()        

    OUT=open(fl,'w')
    OUT.write(ids)
    OUT.close()
    return 1


def downloadAll(storagedir,searchterm=''):
    getData(storagedir=storagedir,searchterm=searchterm)

def downloadBulk(storagedir):
    address='ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline'
    address2='ftp://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles'
    fl='{0}/pubmed18n{1:04d}.xml.gz'
    for i in range(1,1105):
        if os.path.exists(fl.format(storagedir,i)):
            continue
        if i <929:
            url=address
        else:
            url=address2
        try:
            open(fl.format(storagediri))
        except:
            print('Downloading {}'.format(i))
            req=urllib.request.Request(fl.format(url,i))
            response = urllib.request.urlopen(req)
            the_page = response.read()
            f=open(fl.format(storagedir,i),'wb')
            f.write(the_page)
            f.close()
            with gzip.open(fl.format(storagedir,i), 'rb') as fh:
                file_content = fh.read()
                f=open(fl.format(storagedir,i).replace('.gz',''),'wb')
                f.write(file_content)
                f.close()
