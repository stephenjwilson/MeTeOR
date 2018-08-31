'''
Created on Feb 8, 2018

@author: stephen
'''
# Import external packages
import os
# Import internal packages
import crawlPubMedFiles,downloadPubMed,processNetworks,assessNetwork
from NetViz import NetViz
import click

@click.command()
@click.option('--storagedir', default='../MEDLINE', help='Where to store the downlaoded PubMed XML')
@click.option('--datadir', default='../data',
              help="Where the package's data folder is stored ('../data' if running from src)")
@click.option('--resultdir', default='../results',
              help="Where the package's result folder is stored ('../results' if running from src)")
@click.option('--downloadbulk', default=1,
              help="If 1, this option will download PubMed from an FTP server in whole. If 0, this option will download according to a specific query")
@click.option('--searchterm', default='',
              help="If not empty, this option will specify a specific query to build a network off of. downloadBulk must be set to 0.")
def main(storagedir,datadir,resultdir,downloadbulk,searchterm):
        
    # Ensure the storage dir exists
    try:
        os.mkdir(storagedir)
    except:
        pass
    
    # Download Pubmed Files if needed
    if downloadbulk:
        downloadPubMed.downloadBulk(storagedir=storagedir)
    else:
        downloadPubMed.downloadAll(storagedir=storagedir,searchterm=searchterm)
    
    # Crawl Downloaded files and make co-occurrence matrix
    #crawlPubMedFiles.main(datadir, storagedir, resultdir)
    
    # Produce networks
    processNetworks.main(datadir=datadir,resultdir=resultdir)
    
    # Characterize Network
    assessNetwork.run(resultdir=resultdir) # Generates publications per year plot
    
    # Plot network
    NetViz.main(resultdir) # Produce a dot file that is used by graphviz to make a network visualization that is made as a png
    
    
if __name__ == '__main__':

    main()
