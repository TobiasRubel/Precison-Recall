#
# Tobias Rubel | rubelato@reed.edu
# CompBio
#
# This program generates precision and recall data for a given run of an algorithm
# when given a path of the form [algorithm]_[interactome]_[pathway]
# as well as ground truth data about the pathway in question in the form of 
# ...
# It assumes that the directory structure for each named directory is as follows:
#
#   algorithm_interactome_pathway/
#   └── ranked-edges.csv
#
# where the ranked-edges.csv file should look as follows:
# #Tail, Head, k
# change documentation to say it works for many algorithms at once


import utils
import os
import sys
import random
import getopt

import numpy as np
import pandas as pd

def load_df_tab(name:str):
    return pd.read_csv(name,sep='\t')


def get_k(k: int,df: pd.DataFrame) -> pd.DataFrame:
    """
    :k       rank
    :df      dataframe of edges
    :returns sub-dataframe of k-ranked (or less) edges
    """
    try:
        return df[df['k'] <= k].take([0,1],axis=1)
    except:
        return df[df['KSP index'] <= k].take([0,1],axis=1)

def make_edges(df: pd.DataFrame) -> set:
    """
    :df      pandas dataframe
    :returns set of edge tuples 
    """
    return {tuple(x) for x in df.values}

def precision(prediction: set,truth:set,negs: set) -> float:
    prediction = {x for x in prediction if x in truth or x in negs}
    try:
        return len(prediction.intersection(truth))/(len(prediction))
    except:
        #an exception occurs just in case the prediction is neither in 
        #the positive nor negative set. Thus we filter it out.
        return np.NaN

def recall(prediction: set,truth: set,negs: set) -> float:
    prediction = {x for x in prediction if x in truth or x in negs}
    return len(prediction.intersection(truth))/(len(truth))

def pr_edges(predictions: pd.DataFrame,ground: pd.DataFrame,negatives: set,point=False,verbose=False,debug=False) -> dict:
    """
    :prediction dataframe of ranked edges
    :ground     dataframe of "ground truth" edges
    :returns    dictionary of precision keys to recall values
    """
    p = {}
    #turn ground truth into set of edges
    truth = make_edges(ground.take([0,1],axis=1))
    try:
        if point:
            k = max(set(predictions['k']))
            #k = max(set(predictions['KSP index']))
            prediction = make_edges(get_k(k,predictions))
            a = precision(prediction,truth,negatives)
            print(a)
        else:
            for k in set(predictions['k']):
                prediction = make_edges(get_k(k,predictions))
                p[precision(prediction,truth,negatives)] = recall(prediction,truth,negatives)
                b = recall(prediction,truth,negatives)
                p[b] =  a
    except:
        if point:
            k = max(set(predictions['KSP index']))
            prediction = make_edges(get_k(k,predictions))
            a = precision(prediction,truth,negatives)
            b = recall(prediction,truth,negatives)
            p[b] =  a
        else:
            for k in set(predictions['KSP index']):
                if verbose:
                    print('processing k value = {}'.format(k))
                prediction = make_edges(get_k(k,predictions))
                a = precision(prediction,truth,negatives)
                if debug:
                    print('precision: {}'.format(a))
                b = recall(prediction,truth,negatives)
                if debug:
                    print('recall: {}'.format(b))
                if verbose:
                    print('precision: {}\trecall: {}'.format(a,b))
                p[b] =  a
    #default recall = 0 to precision = 1
    p = {k:v for k,v in p.items() if (v != 0 and k != 0)}
    return p

def pr(dname: str,negative:set,edges=True,point=False) -> None:
    """
    :dname       algorithm_interactome_pathway
    :negatives   a set of negatives to use
    :returns     nothing
    :side-effect makes and saves precision-recall data for dname
    """
    #fetch prediction
    try:
        if edges:
            predictions = load_df_tab(os.path.join(dname,'ranked-edges.csv'))
        else:
            predictions = load_df_tab(os.path.join(dname,'ranked-nodes.csv'))
    except:
        print('could not find ranked edges or nodes file for {}'.format(dname))
        return
    #fetch ground truth
    try:
        ground = load_df_tab(os.path.join(dname,'ground.csv'))
    except:
        print('could not find ground truth for {}'.format(dname))
        return
    #load interactome and generate negatives
    try:
        interactome = load_df_tab(os.path.join(dname,'interactome.csv'))
    except:
        print('could not find interactome for {}'.format(dname))
        return
    #negative = negatives(interactome,make_edges(ground))
    if edges:
        p = pr_edges(predictions,ground,negative,point)
    else:
        print('node precision recall not yet implemented. come back soon!')
    #sort dictionary
    p = {k: v for k, v in sorted(p.items(), key=lambda item: item[1])}
    df = pd.DataFrame({'recall':list(p.keys()),'precision':list(p.values())})
    df.to_csv(os.path.join(dname,'pr.csv'),index=False)
    nf = pd.DataFrame({'negatives':list(negative)})
    nf.to_csv(os.path.join(dname,'negatives.csv'),index=False)

def negatives(interactome: pd.DataFrame,positives: set,num:int=0,bivalent=False) -> set:
    """
    :interactome dataframe of interaction data
    :num         number of negatives to randomly generate
    :positives   set of positives
    :returns     k probable negatives
    """
    if num == 0:
        num = len(positives)*50
    edges = make_edges(interactome.take([0,1],axis=1))
    edges = edges - positives
    if bivalent:
        return edges
    return set(random.sample(list(edges),k=num))



def main(argv: str) -> None:
    """
    :args        first argument should be path where the directories reside
                 the other arguments should be directories as per the format
                 requirements at the top of this document
    :returns     nothing
    :side-effect makes precision recall if possible
    """
    try:
        path,directories = argv[1],argv[2:]
        print('generating precision recall data for the following: {}'.format(directories))
    except:
        print('arguments are required...')
        return
    try:
        os.chdir(path)
    except:
        print("path either doesn't exist or could not be accessed.")
    #this is a hack to get the negatives. It should probably be revised.
    #ultimately it would be best if there was one place where interactomes,grounds were stored
    interactome = load_df_tab(os.path.join(directories[0],'interactome.csv'))
    ground = load_df_tab(os.path.join(directories[0],'ground.csv'))
    negative = negatives(interactome,make_edges(ground))
    for d in directories:
        try:
            p = os.path.join(d,'config.conf')
            print(p)
            conf = pd.read_csv(os.path.join(d,'config.conf'),sep=' = ')
        except Exception as e:
            print('no conf file for {}'.format(d))
            print(e)
            return e
        POINT = conf[conf['value'] == 'POINT']['bool'].bool()
        pr(d,negative,True,POINT)



if __name__ == '__main__':
    main(sys.argv)


    
