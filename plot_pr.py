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

def pr(name: str) -> (list,list):
    """
    :name    name of directory
    :returns recall, precision
    """
    #fetch precision and recall
    df = utils.read_df(name,'pr.csv')
    recall = list(df['recall'])
    precision = list(df['precision'])
    return recall,precision


def plot(lat: list, spath: str) -> None:
    """
    :lat         list of directory names
    :spath       name of save path
    :returns     nothing
    :side-effect saves a plot to spath
    """
    #initialize pyplot figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #plot each precision recall plot
    for l in lat:
        #get algorithm name for legend
        lname = l.split('_')[0]
        #plot
        ax.plot(*pr(l),label=lname)
    #format figure globally
    ax.legend()
    title = ' '.join(lat[0].split('_')[1:])
    plt.title(title)
    plt.ylabel('Precision')
    #plt.xticks([])
    #plt.yticks([])
    plt.xlabel('Recall')
    plt.grid(linestyle='--')
    #save the plot
    sname = '-'.join(lat)+'.png'
    #in order to incorporate the save path 
    #some more work needs to be done.
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
        plot(directories,spath)
    else:
        print('coherence could not be established. Terminating...')
        





if __name__ == '__main__':
    main(sys.argv)






