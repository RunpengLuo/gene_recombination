from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import pairwise_distances

from kneed import KneeLocator
import numpy as np
import utils.external as external
import os



def get_distance_matrix(read_feature_dict: dict, read_ids: list):
    dist_mat = [[-1 for _ in range(len(read_ids))] for _ in range(len(read_ids))]

    for i, id1 in enumerate(read_ids):
        for j, id2 in enumerate(read_ids):
            if dist_mat[i][j] != -1:
                continue
            if i == j:
                dist_mat[i][j] = 0
            else:
                # reduced bitwise difference
                dist_mat[i][j] = bin(read_feature_dict[id1] ^ read_feature_dict[id2]).count('1')
                dist_mat[j][i] = dist_mat[i][j]
        if i % 1000 == 0:
            print("Up until {0}th turn".format(i))
    return pairwise_distances(dist_mat, metric="precomputed")

def dbscan_clustering(dist_mat: list, min_pts=5):
    print("Pick optimal eps based on nearest k-neighbor distance graph")
    neigh = NearestNeighbors(n_neighbors=min_pts, metric="precomputed")
    nbrs = neigh.fit(dist_mat)
    distances, indices = nbrs.kneighbors(dist_mat)
    distances = np.sort(distances, axis=0)[:,1]
    
    kneedle = KneeLocator(range(1,len(distances)+1), distances, S=1.0, curve="convex", direction="increasing")
    print("Optimal knee: ", round(kneedle.knee, 3))
    print("Optimal epsilon: ", round(kneedle.knee_y, 3))

    print("start DBSCAN clustering")
    clustering = DBSCAN(metric='precomputed', eps=kneedle.knee_y, min_samples=min_pts).fit(dist_mat)

    labels = clustering.labels_

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)

    print("Estimated number of clusters: %d" % n_clusters_)
    print("Estimated number of noise points: %d" % n_noise_)

    return labels, n_clusters_, n_noise_, min_pts, kneedle.knee_y

# legacy
def cd_hit_clustering(read_dict: dict):
    """Use cd-hit to perform clustering

    Args:
        read_dict (dict): _description_
    """
    gfname = external.outdir + "temp_reads.fasta"
    ofname = external.outdir + "cd_hit_res"
    os.system("echo "" > {0}".format(gfname))
    gfd = open(gfname, "w")
    for read_id, read_seq in read_dict.items():
        gfd.write(">{0}\n{1}\n".format(read_id, read_seq))
    gfd.close()

    # run cd-hit
    os.system("{0} -i {1} -o {2} -c 0.90 -n 5".format(external.cd_hit, gfname, ofname))

    return