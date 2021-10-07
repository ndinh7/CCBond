#! /usr/bin/env python

# 20190301
# $ python ./scripts/find_acid.py
# Will Corcoran, Joseph Ni
"""Find select acids within clustered ZINC acids"""

import os
from Bio import SeqIO
from rdkit.Chem import AllChem as Chem
from rdkit import DataStructs
import numpy as np


def read_clusters(input_dir):
    """Read all cluster information"""

    output_smiles_dict = {}
    output_cluster_dict = {}
    output_all_dict = {}

    cluster_txt = [c for c in os.listdir(input_dir) if c.endswith('')]

    # for each cluster
    for cluster in cluster_txt:
        output_all_dict[cluster] = []
        # for all zincs
        for seq_record in SeqIO.parse(os.path.sep.join([input_dir, cluster]), 'fasta'):
            # trace all zinc smiles
            smiles_temp = Chem.MolToSmiles(Chem.MolFromSmiles(str(seq_record.seq)))
            zinc_temp = ''.join(str(seq_record.id).split('|'))
            output_smiles_dict[smiles_temp] = zinc_temp
            # locate each cluster
            output_cluster_dict[smiles_temp] = cluster
            # find all clusters
            output_all_dict[cluster].append(zinc_temp)

    return output_smiles_dict, output_cluster_dict, output_all_dict


def find_closest(smiles, input_smiles_list, fps_input):
    """
    Find closest smiles of substrate
    :returns smiles: closest smiles of substrate that exists in zinc
    :rtype smiles str
    """

    # if directly in zinc compounds
    if smiles in input_smiles_list:
        return smiles

    # else find closest cpd
    tc_list = DataStructs.BulkTanimotoSimilarity(
        Chem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles), 2), fps_input)
    return input_smiles_list[np.argmax(tc_list)]


def get_cluster(input_smiles, input_smiles_dict, input_cluster_dict, input_all_dict):
    """
    Return detailed info
    :return query_zinc: zinc id of closest substrate
    :rtype zinc id str
    :return querY_cluster: cluster no. that substrate belongs to
    :rtype cluster no. int
    :return query_centroid: zinc id of cluster centroid
    :rtype zinc id str
    """

    return input_smiles_dict[input_smiles], input_cluster_dict[input_smiles], \
           input_all_dict[input_cluster_dict[input_smiles]][0]


if __name__ == '__main__':

    # user input
    # zinc cluster directory
    zinc_dir = os.path.sep.join(['cluster', '20190219_substratecluster_01'])
    # project name
    project_name = 'substrate'
    # smiles entry
    smiles_entry = 'CC(=O)C(=O)O'

    # preprocess info
    smiles_dict, cluster_dict, all_dict = read_clusters(os.path.sep.join([zinc_dir, 'clusters_' + project_name]))

    # sorted smiles and fps list
    smiles_list = sorted(smiles_dict.keys())
    fps_list = [Chem.GetMorganFingerprintAsBitVect(m, 2) for m in [Chem.MolFromSmiles(s) for s in smiles_list]]

    # find closest smiles of substrate
    closest_smiles = find_closest(smiles_entry, smiles_list, fps_list)

    # output detailed info
    query_zinc, query_cluster, query_centroid = get_cluster(closest_smiles, smiles_dict, cluster_dict, all_dict)
