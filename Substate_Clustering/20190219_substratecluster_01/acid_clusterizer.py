#! /usr/bin/env python

# 20190219
# $ python ./scripts/acid_clusterizer.py
# Will Corcoran, Joseph Ni
"""Clusterizing all alpha keto acids from ZINC"""


import os
from rdkit.Chem import AllChem as Chem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
import pandas as pd
import numpy as np


def read_zinc(dir_input):
    """Read zinc and return smiles, fingerprints"""

    mol_list = []
    zinc_df_output = pd.read_csv(dir_input, sep=' ').set_index('zinc_id')

    for k, v in zinc_df_output.iterrows():
        mol_list.append(Chem.MolFromSmiles(v['#smiles']))

    # generate fingeprints: Morgan fingerprint with radius 2
    return zinc_df_output, [Chem.GetMorganFingerprintAsBitVect(m, 2) for m in mol_list]


def ClusterFps(fps, cutoff=0.2):
    """RDKit"""

    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1-x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)

    return dists, cs


def process_distmtx(dists, cs, ids):
    """Process linear distance into distance matrix with only centroids"""

    dists = np.array(dists)
    centroids_temp = np.array(sorted([c[0] for c in cs]))
    distmtx_full = np.zeros([len(ids), len(ids)])
    distmtx_output = []

    # first compute complete matrix as lower triangle for fast query of centroids
    for i in range(1, len(ids)):
        distmtx_full[i, range(0, i)] = dists[range(0, i)]
        dists = dists[i::]

    # then compute distmtx as list of lists for centroids as lower triangle
    for i, c in enumerate(centroids_temp):
        distmtx_output.append(list(distmtx_full[c, centroids_temp[range(0, i)]]) + [0])

    return distmtx_output


def process_headers(ids, cs):
    """Process centroid headers"""

    centroids_argsort = np.argsort([c[0] for c in cs])

    ids_temp = np.array(ids)[[c[0] for c in cs]][centroids_argsort]
    size_temp = np.array([len(c) for c in cs])[centroids_argsort]

    header_output = ['%s;size=%s;' % ('ZINC|' + ids_temp[i].lstrip('ZINC') + '|', size_temp[i]) for i in range(len(cs))]
    substrate_output = np.array([list(np.array(ids)[list(c)]) for c in cs])[centroids_argsort]

    return header_output, substrate_output


def output_results(distmtx_input, header_input, substrate_input, zinc_df_input, distmtx_path, centroids_path, clusters_dir):
    """Output to distmtx, centroids to file"""

    # output distance matrix
    with open(distmtx_path, 'w') as f:
        for i, dists in enumerate(distmtx_input):
            f.write('%s\t%s\n' % (header_input[i], '\t'.join([str(d) for d in dists])))
        f.write('\t' + '\t'.join(header_input) + '\n')

    # output centroid smiles
    with open(centroids_path, 'w') as f:
        for i, s in enumerate(substrate_input):
            f.write('>%s\n%s\n' % (header_input[i], zinc_df_input.loc[s[0]]['#smiles']))

    # output detailed centroid info in folder
    try:
        os.mkdir(clusters_dir)
    except FileExistsError:
        pass
    for i, substrates in enumerate(substrate_input):
        with open(os.path.sep.join([clusters_dir, str(i)]), 'w') as f:
            for s in substrates:
                f.write('>%s\n%s\n' % ('ZINC|' + s.lstrip('ZINC') + '|', zinc_df_input.loc[s]['#smiles']))


if __name__ == '__main__':

    # user input
    # zinc file directory
    zinc_path = os.path.sep.join(['zinc', 'zinc_alpha_keto_acids.smi'])
    # output distmtx path
    output_distmtx_path = os.path.sep.join(['zinc', 'output', 'substrate_distmtx.txt'])
    # output headers path
    output_headers_path = os.path.sep.join(['zinc', 'output', 'substrate_centroids.txt'])
    # output cluster directory
    output_clusters_dir = os.path.sep.join(['zinc', 'output', 'clusters_substrate'])

    # clustering of all alpha keto acids based on chemical similarity
    zinc_df, zinc_fps = read_zinc(zinc_path)
    dist_linear, cluster_tuple = ClusterFps(zinc_fps, 0.7)
    zinc_ids = list(zinc_df.index)

    # process results for Jon's cluster analysis code
    distmtx_arr = process_distmtx(dist_linear, cluster_tuple, zinc_ids)
    header_list, substrate_list = process_headers(zinc_ids, cluster_tuple)

    # output
    output_results(distmtx_arr, header_list, substrate_list, zinc_df,
                   output_distmtx_path, output_headers_path, output_clusters_dir)
    # cluster tuple info
    print(cluster_tuple)
