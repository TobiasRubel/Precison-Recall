#
# Tobias Rubel | rubelato@reed.edu
# Reed CompBio
#
# This program generates composite predictions for all pathways
# for all methods.
#
# It is designed to run without arguments, and can be configured 
# by configuring the global values at the top.
#
# There is a chance that future versions will take arguments. email 
# me if you could make use of such a feature.
import os
import re
import pandas as pd
import shutil

#configure which pathways are being unioned
PATHWAYS = {x.split('-')[0] for x in os.listdir('../../Pathways') if '-nodes' in x}
#the following 3 pathways aren't processable by PathLinker given the 2018 PPI, so we discard them.
PATHWAYS.remove('IL')
PATHWAYS.remove('RAGE')
PATHWAYS.remove('ID')
#whoops
PATHWAYS.remove('RANKL')
print('processing pathways: {}'.format(", ".join(PATHWAYS)))

#configure path to data
DPATH = 'data'

INTERACTOME = '../../../../Interactomes/PathLinker_2018_human-ppi-weighted-cap0_75.txt'
#fetch which algorithms we have data for

ALGORITHMS = {x.split('_')[0] for x in os.listdir(DPATH)}
print('The following algorithms exist: {}'.format(", ".join(ALGORITHMS)))

def sanity_check(A:str) -> bool:
    """
    :A       name of algorithm
    :returns whether data exists for A for all pathways 
    """
    global DPATH, PATHWAYS
    dirs = [x for x in os.listdir(DPATH) if A in x]
    for P in PATHWAYS:
        print(P)
        if not any(P in d for d in dirs):
            return False
    return True


def join(A: str) -> pd.DataFrame:
    global DPATH
    dirs = [x for x in os.listdir(DPATH) if re.match("^{}_".format(A),x) and not 'composit' in x]
    print(A)
    print(dirs)
    print('*'*50)
    edge_lists = [os.path.join(DPATH,os.path.join(d,'ranked-edges.csv')) for d in dirs]
    df = pd.read_csv(edge_lists[0],sep='\t')
    for d in edge_lists[1:]:
        df = pd.concat([df,pd.read_csv(d,sep='\t')])
    #for ranked edgelists...
    ranked = 'KSP index' in df.columns or 'rank' in df.columns
    if ranked:
        try:
            df = df.rename(columns={'KSP index':'rank'})
            #print(df)
        except Exception as e:
            print('*'*50)
            print(e)
            #print(df)
            pass
        df = df.sort_values(by='rank')
        df = df.drop_duplicates().reset_index(drop=True)
    return df

def join_pathways(PLAT: list) -> pd.DataFrame:
    pathways = [os.path.join('../../Pathways','{}-edges.txt'.format(P)) for P in PLAT]
    df = pd.read_csv(pathways[0],sep='\t')
    for d in pathways[1:]:
        df = pd.concat([df,pd.read_csv(d,sep='\t')])
    df = df.drop_duplicates().reset_index(drop=True)
    return df

def join_negatives(dpath,PLAT:list):
    print(PLAT)
    print(dpath)
    neg = os.path.join(next(x for x in os.listdir(dpath) if PLAT[0] in x),'negatives.csv')
    neg = os.path.join(dpath,neg)
    df = pd.read_csv(neg,sep='\t')
    for P in PLAT[1:]:
        df = pd.concat([df,pd.read_csv(os.path.join(dpath,os.path.join(next(x for x in os.listdir(dpath) if P in x),'negatives.csv')))])
    print(df)
    return df

def main():
    global PATHWAYS, DPATH, ALGORITHMS, INTERACTOME
    #perform sanity check to ensure that for all algorithms we have all data
    unsafe_algorithms = []
    for ALG in ALGORITHMS:
        if sanity_check(ALG):
            pass
        else:
            unsafe_algorithms.append(ALG)
    if len(unsafe_algorithms) != 0:
        print('the following algorithms are unsafe: {}.\nPlease ensure that data exists for all algorithms.'.format(', '.join(unsafe_algorithms)))
        return
    #join negatives
    negs = join_negatives(DPATH,list(PATHWAYS))
    #make composit pathway
    composit_pathway = join_pathways(PATHWAYS)
    #make empty directories
    ddict = dict()
    for ALG in ALGORITHMS:
        i = 'PathLinker_2018_human-ppi-weighted-cap0_75'
        name = '_'.join([ALG,i,'composit',str(10000)]) #placeholder k for now, since I need to rewrite kwargs
        DEST = os.path.join(DPATH,name)
        ddict[ALG] = DEST
        try:
            os.mkdir(DEST)
        except:
            print('directory already existed. Not overwriting.')
    #populate directories with union predictions
    for ALG in ALGORITHMS:
        #os.remove(os.path.join(ddict[ALG],'*'))
        join(ALG).to_csv(os.path.join(ddict[ALG],'ranked-edges.csv'),sep='\t',index=False)
        composit_pathway.to_csv(os.path.join(ddict[ALG],'ground.csv'),sep='\t',index=False)
        negs.to_csv(os.path.join(ddict[ALG],'negatives.csv'),sep='\t',index=False)
        try:
            #os.remove(os.path.join(ddict[ALG],'interactome.csv'))
            os.symlink(INTERACTOME,os.path.join(ddict[ALG],'interactome.csv'))
        except Exception as e:
            print(e)
            pass
        shutil.copy('config.conf',os.path.join(ddict[ALG],'config.conf'))







if __name__ == "__main__":
    main()
