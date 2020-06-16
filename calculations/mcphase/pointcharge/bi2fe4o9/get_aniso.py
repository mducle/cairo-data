import numpy as np
import sys
from sklearn.cluster import DBSCAN

def get_aniso(filename):
    dataset = np.loadtxt(filename)
    idx = np.where(dataset[:,7]>=(np.max(dataset[:,7])*0.999))
    clusters = DBSCAN(eps=0.1, min_samples=2).fit(dataset[idx[0],4:7])
    cluster_labels = set(clusters.labels_)
    dirn = []
    for icl in cluster_labels:
        if icl > -1:   # -1 is label for "noise"
            dirn.append( np.mean(dataset[idx[0][np.where(clusters.labels_ == icl)], 4:7], 0) )
            print(dirn[-1])

if __name__ == '__main__':
    get_aniso(sys.argv[1]) 
