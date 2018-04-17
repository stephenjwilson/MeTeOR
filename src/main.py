'''
Created on Feb 8, 2018

@author: stephen
'''
# Import external packages
import os
# Import internal packages
import crawlPubMedFiles,downloadPubMed,processNetworks,assessNetwork
from NetViz import NetViz
def main(storagedir,datadir,resultdir):
        
    # Ensure the storage dir exists
    try:
        os.mkdir(storagedir)
    except:
        pass
    
    # Download Pubmed Files if needed
    # downloadPubMed.downloadAll(storagedir=storagedir) # Put up on OSF???
    
    # Crawl Downloaded files and make co-occurrence matrix
    #crawlPubMedFiles.main(datadir, storagedir, resultdir)
    
    # Produce networks
    #processNetworks.main(datadir=datadir,resultdir=resultdir)
    
    # Characterize Network
    #assessNetwork.run(resultdir=resultdir) # Generates publications per year plot
    
    # Plot network
    NetViz.main(resultdir) # Produce a dot file that is used by graphviz to make a network visualization that is made as a png
    
    
if __name__ == '__main__':
    storagedir='/home/stephen/ExtraDrive1/MeSH_V3'
    datadir='../data' # Where the package's data folder is stored ('../data' if running from src)
    resultdir='../results' # Where the package's results folder is stored ('../results' if running from src)
    main(storagedir,datadir,resultdir)
