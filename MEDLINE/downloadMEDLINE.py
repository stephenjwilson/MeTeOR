import urllib.request, urllib.error, urllib.parse
import gzip
address='ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline/'
address2='ftp://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/'
fl='pubmed18n{0:04d}.xml.gz'


for i in range(1,1105):
    if i <929:
        url=address
    else:
        url=address2
    try:
        open(fl.format(i))
    except:
        print('Downloading {}'.format(i))
        req=urllib.request.Request(url+fl.format(i))
        response = urllib.request.urlopen(req)
        the_page = response.read()
        f=open(fl.format(i),'wb')
        f.write(the_page)
        f.close()
        with gzip.open(fl.format(i), 'rb') as fh:
            file_content = fh.read()
            f=open(fl.format(i).replace('.gz',''),'wb')
            f.write(file_content)
            f.close()
		
