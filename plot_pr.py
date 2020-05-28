#
# Tobias Rubel | rubelato@reed.edu
# CompBio
#
# This program creates precision recall plots given a path followed by a save path,
# followed by a list of directory names to plot from.
#
# It assumes that the names are as follows: [algorithm]_[interactome]_[pathway]
# It assumes that the directory structure for each named directory is as follows:
#
#   algorithm_interactome_pathway/
#   ├── negatives.csv
#   └── pr.csv
import re
import sys
import os
import utils
from itertools import cycle

import matplotlib.pyplot as plt

# routines to verify that it makes sense to compare the precision and recall of
# the different algorithms by ensuring that they were computed on the same
# interactome, that they were predicting the same pathway, and that they used
# the same negatives.

def verify_negatives(lat: list,node_motivation=False,verbose=False) -> bool:
    """
    :lat     list of directory names
    :returns True iff the runs were made with the same negatives
    """
    if node_motivation: ## skip this for now if node_motivation is specified.
        negfilenames = ['negatives-nodes.csv','negatives-nodes-ignoreadj.csv']
    else:
        negfilenames = ['negatives.csv']

    toReturn = {}
    for negfilename in negfilenames:
        #get one of the lists of negatives as a dataframe
        df = utils.read_df(lat[0],negfilename)

        if verbose:
            for j in lat[1:]:
                print(os.path.join(j,negfilename),df.equals(utils.read_df(j,negfilename)))
        #exploiting transitivity of identity, compare the first against all others
        toReturn[negfilename] = all(df.equals(i) for i in [utils.read_df(j,negfilename) for j in lat[1:]])
    return all(list(toReturn.values()))

def verify_pathway(lat: list) -> bool:
    """
    :lat     list of directory names
    :returns True iff the runs were made on the same pathway
    """
    #helper function which gives us the pathway name given a formatted filename
    f = lambda x: x.split('_')[2]
    #get one of the pathways
    p = f(lat[0])
    #exploiting transitivity of identity, compare the first against all others
    return all(p == q for q in [f(l) for l in lat[1:]])

def verify_interactome(lat: list) -> bool:
    """
    :lat     list of directory names
    :returns True iff the runs were made on the same interactome
    """
    #helper function which gives us the pathway name given a formatted filename
    f = lambda x: x.split('_')[1]
    #get one of the pathways
    p = f(lat[0])
    #exploiting transitivity of identity, compare the first against all others
    return all(p == q for q in [f(l) for l in lat[1:]])

def verify_coherence(lat: list,node_motivation=False) -> bool:
    """
    :lat     list of directory names
    :returns True iff all of the coherence checks pass
    """
    return all([verify_negatives(lat,node_motivation),verify_pathway(lat),verify_interactome(lat)])

# routines that handle the plotting of the precision recall plots

def pr(name: str,edges=True,ignore_adj=False) -> (list,list):
    """
    :name    name of directory
    :edges   (trivalent) boolean as to whether we are plotting edges,nodes or both
    :returns recall, precision
    """
    #fetch precision and recall
    if edges == True:
        df = utils.read_df(name,'pr-edges.csv')
    elif edges == False:
        if ignore_adj == True:
            df = utils.read_df(name,'pr-nodes-ignoreadj.csv')
        else:
            df = utils.read_df(name,'pr-nodes.csv')
    df = df.sort_values(by=['recall','precision'],ascending=[True,False])
    #df = df.sort_values('precision',ascending=False)
    recall = list(df['recall'])
    precision = list(df['precision'])
    return recall,precision


def plot(lat: list, spath: str,params=False,edges=True) -> None:
    """
    :lat         list of directory names
    :spath       name of save path
    :edges       (trivalent) boolean as to whether we are plotting edges,nodes or both
    :returns     nothing
    :side-effect saves a plot to spath
    """
    #initialize pyplot figure
    markers = iter(['o','v','^','<','>','1','2','3','4','8','s','p','P','*','h','H','+','x','X','D','d','|','_'][:len(lat)]*50)
    if (edges == True or edges == False):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        #plot each precision recall plot
        for l in lat:
            print(l)
            #get algorithm name for legend
            #lname = '-'.join([l.split('_')[0],l.split('_')[-1]])
            if params:
                lname = l.split('/')[-1]
            else:
                lname = l.split('/')[-1].split('_')[0]
            #plot
            ax.plot(*pr(l,edges),label=lname,marker=next(markers),alpha=0.7)
    elif edges == '#':
        fig = plt.figure(figsize=(12,5))
        ax = fig.add_subplot(121)
        #plot each precision recall plot
        for l in lat:
            #get algorithm name for legend
            lname = l.split('_')[0]
            #plot
            ax.plot(*pr(l,True),label=lname,marker=next(markers),alpha=0.7)
        cax = fig.add_subplot(122)
        #plot each precision recall plot
        for l in lat:
            #get algorithm name for legend
            if params:
                lname = l.split('/')[-1]
            else:
                lname = l.split('/')[-1].split('_')[0]
            #plot
            cax.plot(*pr(l,False),label="_nolegend_",marker=next(markers),alpha=0.7)
    #format figure globally
    #ax.legend()
    if params:
        plt.legend(loc='upper right')
    else:
        plt.legend(loc='center right')
    title = lat[0].split('/')[-1].split('_')[2] + ' Pathway'
    #fig.suptitle(title,fontsize=16)
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_title(title + ' Interactions')
    ax.grid(linestyle='--')
    ax.set_xlim(0,1)
    ax.set_ylim(0.1)
    if edges == '#':
        cax.set_xlabel('Recall')
        cax.set_ylabel('Precision')
        cax.set_title(title + ' Proteins')
        cax.grid(linestyle='--')
        cax.set_xlim(0,1)
        cax.set_ylim(0,1)

    #save the plot
    print(lat)
    lat = [x.replace('HybridLinker','HL') for x in lat]
    lat = [x.replace('PerfectLinker','PeL') for x in lat]

    methods = [x.split('/')[-1].split('_') for x in lat]
    infix = [m[0] for m in methods]+methods[0][1:]
    if edges:
        sname = 'edges-'+'-'.join(infix)+'.png'
    else:
        sname = 'nodes-'+'-'.join(infix)+'.png'

    #in order to incorporate the save path
    #some more work needs to be done.
    plt.tight_layout()
    #plt.subplots_adjust(left=0.21)
    #plt.subplots_adjust(top=0.90)
    print(spath,sname)
    plt.savefig(os.path.join(spath,sname))

def plot_composite(lat: list, spath: str,) -> None:
    """
    :lat         list of directory names
    :spath       name of save path
    :edges       (trivalent) boolean as to whether we are plotting edges,nodes or both
    :returns     nothing
    :side-effect saves a plot to spath
    """
    #initialize pyplot figure
    markers = iter(['o','v','^','<','>','1','2','3','4','8','s','p','P','*','h','H','+','x','X','D','d','|','_'][:len(lat)]*50)
    fig, axs = plt.subplots(2, 2,figsize=(25,15))
    plotloc = axs.flat
    print(axs)
    #turn list of methods into list of tuples of things to plot together:
    combine_key = {'HybridLinker-BTB':'BowtieBuilder','HybridLinker-SP':'ShortestPaths','HybridLinker-PL':'PathLinker','HybridLinker-RN':'ResponseNet','HybridLinker':'PathLinker'}
    partner = lambda x: next(y for y in lat if re.match('^{}_'.format(combine_key[x]),y))
    old_lat = lat
    lat = [(x,partner(x.split('_')[0])) for x in lat if x.split('_')[0] in combine_key]
    print(lat)
    #plot each precision recall plot
    for l in lat:
        loc = next(plotloc)
        print(loc)
        #get algorithm name for legend
        l1name = l[0].split('_')[0]
        #plot
        loc.plot(*pr(l[0],True),label=l1name,marker=next(markers),alpha=0.7)
        #get algorithm name for legend (once more with feeling!)
        l2name = l[1].split('_')[0]
        #plot
        loc.plot(*pr(l[1],True),label=l2name,marker=next(markers),alpha=0.7)
    #cax = fig.add_subplot(122)
    #format figure globally
    #ax.legend()
    fig.legend(loc='center left')
    title = 'Composite Interaction Performance across 29 Pathways'
    fig.suptitle(title,fontsize=16)
    for ax in axs.flat:
        ax.set_xlabel('Recall')
        ax.set_ylabel('Precision')
        ax.grid(linestyle='--')
    #toggle axes
    #plt.setp(axs, xlim=(0,1), ylim=(0,1))
    #save the plot
    old_lat = [x.replace('HybridLinker','HL') for x in old_lat]
    old_lat = [x.replace('PerfectLinker','PeL') for x in old_lat]
    sname = 'composite-blowup.pdf'
    #in order to incorporate the save path
    #some more work needs to be done.
    plt.tight_layout()
    plt.subplots_adjust(left=0.21)
    plt.subplots_adjust(top=0.90)
    plt.savefig(os.path.join('../',spath,sname))

def plot_node_motivation(lat: list, spath: str) -> None:
    """
    :lat         list of directory names
    :spath       name of save path
    :edges       (trivalent) boolean as to whether we are plotting edges,nodes or both
    :returns     nothing
    :side-effect saves a plot to spath
    """
    #initialize pyplot figure
    markers = iter(['o','v','^','<','>','1','2','3','4','8','s','p','P','*','h','H','+','x','X','D','d','|','_'][:len(lat)]*50)
    fig = plt.figure()
    ax = fig.add_subplot(111)

    color_cycle = cycle(['r','g','b', 'c', 'm', 'y', 'k'])
    #plot each precision recall plot
    for l in lat:
        #get algorithm name for legend
        lname = l.split('/')[-1].split('_')[0]
        #plot
        this_color = next(color_cycle)
        recall,precision = pr(l,False,False)
        recall2,precision2 = pr(l,False,True)
        if len(recall) == 1:
            this_marker = next(markers)
            ax.plot(recall,precision,this_marker,color=this_color,label=lname,marker=this_marker,ms=10,alpha=0.3)
            ax.plot(recall2,precision2,this_marker,color=this_color,label=lname +' (nearby)',marker=this_marker,ms=10,alpha=0.9)
            ax.plot([recall[0],recall2[0]],[precision[0],precision2[0]],'--',color=this_color,alpha=0.3)
        else:
            ax.plot(recall,precision,color=this_color,label=lname,lw=4,alpha=0.3)
            ax.plot(recall2,precision2,color=this_color,label=lname +' (nearby)',lw=4,alpha=0.9)
            ax.plot([recall[0],recall2[0]],[precision[0],precision2[0]],'--',color=this_color,alpha=0.3)
            ax.plot([recall[-1],recall2[-1]],[precision[-1],precision2[-1]],'--',color=this_color,alpha=0.3)

    #format figure globally
    #ax.legend()

    title = ' '.join(lat[0].split('_')[1:])

    ax.set_xlabel('Node Recall')
    ax.set_ylabel('Node Precision')
    ax.set_title('Reconstructing Wnt Proteins',fontsize=14)
    ax.set_ylim(0,1.01)
    ax.set_xlim(0,1.01)

    methods = [x.split('/')[-1].split('_') for x in lat]
    infix = [m[0] for m in methods]+methods[0][1:]
    sname = 'node-motivation-' + '-'.join(infix)+'.png'

    #in order to incorporate the save path
    #some more work needs to be done.
    #plt.plot([], [], ' ', label="*No edge-adjacent\nnegatives")
    plt.legend(loc='lower left',fontsize=8)
    plt.tight_layout()
    #plt.subplots_adjust(left=0.21)
    #plt.subplots_adjust(top=0.90)
    plt.savefig(os.path.join(spath,sname))
    print('writing to %s'% (os.path.join(spath,sname)))
    sname = sname.replace('.png','.pdf')
    plt.savefig(os.path.join(spath,sname))
    os.system('pdfcrop %s %s' % (os.path.join(spath,sname),os.path.join(spath,sname)))
    print('writing to %s'% (os.path.join(spath,sname)))

#handle input
def main(args: list) -> None:
    """
    :args        first argument should be path where the directories reside
                 the other arguments should be directories as per the format
                 requirements at the top of this document
    :returns     nothing
    :side-effect plots precision recall if possible
    """
    #get current directory to pass through to plot
    #cdir = os.getcwd()
    #change to working directory
    try:
        path,spath,directories = args[1],args[2],sorted(list(set(args[3:])))
        print('plotting the following: {}'.format(directories))
    except:
        print('arguments are required...')
    directories = [os.path.join(path,d) for d in directories]
    #try:
    #    os.chdir(path)
    #except:
    #    print("path {} either doesn't exist or could not be accessed.".format(path))

    ## these are HARD-CODED - need to make them arguments.
    COMPOSITE=False
    NODE_MOTIVATION=False
    PARAMS=True
    if verify_coherence(directories,NODE_MOTIVATION):
        if COMPOSITE:
            plot_composite(directories,spath)
        elif NODE_MOTIVATION:
            plot_node_motivation(directories,spath)
        else: # plot regular PR
            plot(directories,spath,PARAMS,True)

    else:
        print('Coherence could not be established. Terminating...')
        print('\nVerifying Negatives:')
        print(verify_negatives(directories,NODE_MOTIVATION,True))
        print('\nVerifying Pathway:')
        print(verify_pathway(directories))
        print('\nVerifying Interactome:')
        print(verify_interactome(directories))






if __name__ == '__main__':
    main(sys.argv)
