from collections import Counter
import numpy as np
import pylab as pl
from IPython import embed
import powerlaw
import scipy.sparse
fontsize=16

font = {'weight' : 'bold',
        'size'   : fontsize}

import matplotlib
matplotlib.rc('font', **font)


def articlesPerYear(fl='meta.txt',resultdir='../results',title='MeTeOR'):
    myears=[]
    f=open('{}/{}'.format(resultdir,fl))
    for line in f:
        line=line.strip('\n').split('\t')
        year=[ int(x.replace('Year:', '')) for x in line[1].split(';') if 'Year' in x][0] # the year
        myears.append(year)
    cyears=Counter(myears)
    x=cyears.keys()
    y=np.array(list(cyears.values()))/1000.
    pl.scatter(x,y)
    pl.xlabel('Year')
    pl.ylabel('Number of Articles (thousands)')
    pl.tight_layout()
    pl.savefig('{}/ArticleperYear_{}'.format(resultdir,title+'.pdf'))


def degreeDistribution(fl='coOccurrence.npz',mappingfl='UIsMapped.txt',resultdir='../results',title='MeTeOR'):
    
    # Load scipy co-occurrence matrix
    mat=scipy.sparse.load_npz('{}/{}'.format(resultdir,fl)).tocsr()
    # get degree per node by summing row on one dimension
    degrees=np.sum(mat,axis=0)
    # subset to nodes with degree at least 1
    keep=np.where(degrees>0)[1]
    mat=mat[keep,:]
    mat=mat[:,keep]
    degrees=np.sum(mat,axis=0)
    nnz=mat.nnz
    degreelog=open('{}/degreelog.txt'.format(resultdir),'w')
    s='Total edges = {}\n'.format(nnz)
    
    vals=np.array(sorted(np.array(degrees)[0]))
    s+='Total nodes = {}\n'.format(mat.shape[0])
    s+='Mean Degree is {}\n'.format(np.mean(vals))
    s+='Median Degree is {}\n'.format(np.median(vals))
    degreelog.write(s)
    print(s)
    s=''
    f = pl.figure(figsize=(8,11))
    data = vals
    data_inst = 1
    units = 'MeSH Frequency'
    plot_basics(data, data_inst, f, units)
    
    results=powerlaw.Fit(vals)
    s+='a={},xmin={}\n'.format(results.power_law.alpha,results.power_law.xmin)
    
    _, p = results.distribution_compare('power_law', 'lognormal')    
    s+='Prob it is power_law over lognormal {}\n'.format(p)
    _, p = results.distribution_compare('power_law', 'exponential')
    s+='Prob it is power_law over exponential {}\n'.format(p)    
    degreelog.write(s)
    print(s)
    s=''
    
    ### Ignore non-labeled nodes 
    labeled=[]
    for line in open('{}/{}'.format(resultdir,mappingfl)).read().split('\n'):
        line=line.split('\t')
        if len(line)==3:
            labeled.append(line[1]!='Unknown')
    labeled=np.array(labeled)
    labeled=labeled[keep]
    vals=vals[labeled]

    data = vals
    data_inst = 2
    units = 'MeSH Frequency'
    plot_basics(data, data_inst, f, units)
    
    results=powerlaw.Fit(vals)
    s+='a={},xmin={}\n'.format(results.power_law.alpha,results.power_law.xmin)
    
    _, p = results.distribution_compare('power_law', 'lognormal')    
    s+='Prob it is power_law over lognormal {}\n'.format(p)
    _, p = results.distribution_compare('power_law', 'exponential')
    s+='Prob it is power_law over exponential {}\n'.format(p) 
    
    total=float(len(vals))
    s+='Total nodes = {}'.format(total)
    degreelog.write(s)
    print(s)
    degreelog.close()
    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.3, hspace=.2)
    pl.tight_layout()
    f.savefig('{}/DegreeDist_{}'.format(resultdir,title+'.pdf'))



def run(resultdir='../results/'):
    articlesPerYear(resultdir=resultdir)
    degreeDistribution(resultdir=resultdir)
    
def plot_basics(data, data_inst, fig, units):

    ### Setup ###
    from powerlaw import plot_pdf, Fit, pdf
    import pylab
    pylab.rcParams['xtick.major.pad']='8'
    pylab.rcParams['ytick.major.pad']='8'
    #pylab.rcParams['font.sans-serif']='Arial'
    
    from matplotlib.font_manager import FontProperties
    
    panel_label_font = FontProperties().copy()
    panel_label_font.set_weight("bold")
    panel_label_font.set_size(30.0)
    panel_label_font.set_family("sans-serif")
    n_data = 2
    n_graphs = 4
    annotate_coord = (-.4, .95)
    #############
    
    
    ax1 = fig.add_subplot(n_graphs,n_data,data_inst)
    x, y = pdf(data, linear_bins=True)
    ind = y>0
    y = y[ind]
    x = x[:-1]
    x = x[ind]
    ax1.scatter(x, y, color='r', s=.5)
    plot_pdf(data[data>0], ax=ax1, color='b', linewidth=2)
    from pylab import setp
    setp( ax1.get_xticklabels(), visible=False)

    if data_inst==1:
        ax1.annotate("A", annotate_coord, xycoords="axes fraction", fontproperties=panel_label_font)
    
    ax2 = fig.add_subplot(n_graphs,n_data,n_data+data_inst, sharex=ax1)

    plot_pdf(data, ax=ax2, color='b', linewidth=2)
    fit = Fit(data, xmin=1, discrete=True)
    fit.power_law.plot_pdf(ax=ax2, linestyle='--', color='g')
    _ = fit.power_law.pdf()
    ax2.set_xlim((1,max(x)))

    setp( ax2.get_xticklabels(), visible=False)

    if data_inst==1:
        ax2.annotate("B", annotate_coord, xycoords="axes fraction", fontproperties=panel_label_font)        
        ax2.set_ylabel(u"p(X)")# (10^n)")
        
    ax3 = fig.add_subplot(n_graphs,n_data,n_data*2+data_inst)#, sharex=ax1)#, sharey=ax2)
    fit.power_law.plot_pdf(ax=ax3, linestyle='--', color='g')
    fit.exponential.plot_pdf(ax=ax3, linestyle='--', color='r')
    fit.lognormal.plot_pdf(ax=ax3, linestyle=':', color='r')
    fit.plot_pdf(ax=ax3, color='b', linewidth=2)
    
    ax3.set_ylim(ax2.get_ylim())
    ax3.set_xlim(ax1.get_xlim())
    
    if data_inst==1:
        ax3.annotate("C", annotate_coord, xycoords="axes fraction", fontproperties=panel_label_font)

    ax3.set_xlabel(units)

    
if __name__=='__main__':
    run()
