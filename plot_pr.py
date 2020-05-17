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

import matplotlib.pyplot as plt

# routines to verify that it makes sense to compare the precision and recall of 
# the different algorithms by ensuring that they were computed on the same 
# interactome, that they were predicting the same pathway, and that they used 
# the same negatives.

def verify_negatives(lat: list) -> bool:
    """
    :lat     list of directory names
    :returns True iff the runs were made with the same negatives
    """
    #get one of the lists of negatives as a dataframe
    df = utils.read_df(lat[0],'negatives.csv')
    #exploiting transitivity of identity, compare the first against all others
    return all(df.equals(i) for i in [utils.read_df(j,'negatives.csv') for j in lat[1:]])

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

def verify_coherence(lat: list) -> bool:
    """
    :lat     list of directory names
    :returns True iff all of the coherence checks pass
    """
    return all([verify_negatives(lat),verify_pathway(lat),verify_interactome(lat)])

# routines that handle the plotting of the precision recall plots

def pr(name: str,edges=True) -> (list,list):
    """
    :name    name of directory
    :edges   (trivalent) boolean as to whether we are plotting edges,nodes or both
    :returns recall, precision
    """
    #fetch precision and recall
    if edges == True:
        df = utils.read_df(name,'pr-edges.csv')
    elif edges == False:
        df = utils.read_df(name,'pr-nodes.csv')
    df = df.sort_values(by=['recall','precision'],ascending=[True,False])
    #df = df.sort_values('precision',ascending=False)
    recall = list(df['recall'])
    precision = list(df['precision'])
    return recall,precision


def plot(lat: list, spath: str,edges=True) -> None:
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
            #get algorithm name for legend
            lname = l.split('_')[0]
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
            lname = l.split('_')[0]
            #plot
            cax.plot(*pr(l,False),label="_nolegend_",marker=next(markers),alpha=0.7)
    #format figure globally
    #ax.legend()
    fig.legend(loc='center left')
    title = ' '.join(lat[0].split('_')[1:])
    fig.suptitle(title,fontsize=16)
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_title('Interactions')
    if edges == '#':
        cax.set_xlabel('Recall')
        cax.set_ylabel('Precision')
        cax.set_title('Proteins')
        cax.grid(linestyle='--')
    ax.grid(linestyle='--')
    #save the plot
    lat = [x.replace('HybridLinker','HL') for x in lat]
    lat = [x.replace('PerfectLinker','PeL') for x in lat]
    sname = str(edges)+'-'.join([x.split('_')[0] for x in lat]+lat[0].split('_')[1:])+'.png'
    #in order to incorporate the save path 
    #some more work needs to be done.
    plt.tight_layout()
    plt.subplots_adjust(left=0.21)
    plt.subplots_adjust(top=0.90)
    plt.savefig(os.path.join('../',spath,sname))

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
    cdir = os.getcwd()
    #change to working directory
    try:
        path,spath,directories = args[1],args[2],sorted(list(set(args[3:])))
        print('plotting the following: {}'.format(directories))
    except:
        print('arguments are required...')
    try:
        os.chdir(path)
    except:
        print("path {} either doesn't exist or could not be accessed.".format(path))
    if verify_coherence(directories):
        plot(directories,spath,True)
        #plot_composite(directories,spath)
    else:
        print('coherence could not be established. Terminating...')
        





if __name__ == '__main__':
    main(sys.argv)






