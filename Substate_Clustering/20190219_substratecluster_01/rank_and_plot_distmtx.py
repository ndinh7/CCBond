# Takes a distance matrix (calculated from Calc_distmtx.py) and converts it
# to a numpy array of positions. This is then converted into a ranked list of
#  clusters as a text file with one cluster per line. First argument is
# infile and second is outfile. Third is FASTA file that was used to generate
#  the distance matrix. It must have cluster size annotations for each
# cluster centroid (can get this by using -sizeout option with usearch
# -cluster_fast algorithm).
#
# Example implementation:
# python ./rank_and_plot_distmtx.py distmtx.txt ranked_centroids.txt
#   centroids.fasta -title 'My Title Here' -rank 20 -show 10 -show_pdb True
#   -cluster_file_dir "./cluster_files/" -circle_size_multiplier 10


import numpy as np
from sklearn.manifold import MDS
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Find_PDB_Clusters_SIFTS import *
from get_size import get_size


# Get cluster sizes and generate distance matrix from distance matrix text file
def convert_distance_matrix(filename, num_seqs, out_matrix, in_size_dict):
    line_num = -1
    number = ''
    print(in_size_dict)
    with open(filename, 'r') as inputfile:
        for line in inputfile:
            # Distance matrix text file generated via biopython only includes
            # first 10 characters of the full sequence id. If two sequences
            # have the same first 10 characters for their id (extremely
            # unlikely), then just manually change one of them temporarily.
            
            current_id = '|'.join(line.split('|')[0:2]) + '|'

            print('Line Num: ', line_num, '\tnum_seqs: ', num_seqs)
            if line_num < num_seqs - 1:
                line_num += 1
            else:
                print(labels)
                return out_matrix, size_list, size_sqrt_list, labels
            # Generate cluster labels and sizes for graph
            labels[line_num] = current_id
            size_list[line_num] = in_size_dict[current_id]
            # Arbitrarily multiply by 300 so that circles are larger on graph
            size_sqrt_list[line_num] = 300 * np.sqrt(in_size_dict[current_id])
            # Extract each number from distance matrix text file by finding its
            # decimal point. Then add to our matrix.
            column_num = 0
            for char_num in range(len(line.split('\t')[0]), len(line)):
                if line[char_num] == '.':
                    number += line[char_num - 1]
                    while line[char_num] != '\t' and line[char_num] != '\n':
                        number += line[char_num]
                        char_num += 1
                    number_to_add = float(number)
                    number = ''
                    out_matrix[column_num][line_num] = number_to_add
                    out_matrix[line_num][column_num] = number_to_add
                    column_num += 1


# Get positions (pos) of all clusters based on distance matrix using
# multi-dimensional scaling (mds). Must first convert our list of lists
# (matrix) to a numpy array (np_matrix).
def transform_matrix(in_matrix):
    np_matrix = np.array(in_matrix)
    mds = MDS()
    positions = mds.fit_transform(np_matrix)
    return positions


def get_max_cluster(in_size_dict):
    # Initialize variables
    max_size = 0
    max_cluster = ''
    # Find biggest cluster and its size to calculate weighted size later
    for cluster_id in in_size_dict.keys():
        if in_size_dict[cluster_id] > max_size:
            max_size = in_size_dict[cluster_id]
            max_cluster = cluster_id
    return max_cluster


def get_ranked_list(in_table, number_picked_clusters, number_shown_clusters,
                    starting_cluster, size_weight):
    # Generate ranked list of clusters
    current_cluster = starting_cluster
    ranked_list = []
    i = 0
    # Pairwise maximum distance
    d_pair_max = calc_pairwise_max_distance(in_table)
    # Choose weight for size and distance (must add to 1)
    distance_weight = 1 - size_weight
    while i < number_picked_clusters:
        score_max = 0
        # Maximum sum of distances to all other points from current point
        d_total_max = d_pair_max * (i + 1)
        # Go through each cluster that hasn't yet been picked and calculate its
        # distance to every other cluster that has been picked. Then,
        ranked_list.append(current_cluster)
        in_table.set_value(current_cluster, 'picked', True)
        if i < number_shown_clusters:
            in_table.set_value(current_cluster, 'show', True)
        for cluster in in_table.loc[~in_table['picked']].loc[:, 'id']:
            d = 0
            x1 = in_table.at[cluster, 'x']
            y1 = in_table.at[cluster, 'y']
            # Calculate sum of distances from cluster to all picked clusters
            for picked_cluster in in_table.loc[in_table['picked']].loc[:, 'id']:
                x2 = in_table.at[picked_cluster, 'x']
                y2 = in_table.at[picked_cluster, 'y']
                d += np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2) / d_pair_max
            cluster_score = distance_weight / d_total_max * d + \
                size_weight * size_dict[cluster] / size_dict[starting_cluster]
            if cluster_score > score_max:
                score_max = cluster_score
                current_cluster = cluster
        i += 1
        print(i, 'out of', num_picked_clusters, 'read.')
    return ranked_list


def calc_pairwise_max_distance(in_table):
    max_distance = 0
    x_list = in_table['x']
    y_list = in_table['y']
    # Calculate distance from every point to every other point and record max
    #  distance
    for x1, y1 in zip(x_list, y_list):
        for x2, y2 in zip(x_list, y_list):
            distance = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
            if distance > max_distance:
                max_distance = distance
    return max_distance


def get_pdb_id_dict(ranked_list, directory):
    # Find cluster number for each id
    cluster_file_dict = find_cluster_files(ranked_list,
                                           directory)
    mapping_df = pd.read_csv('/projects/p30511/utils/Cluster_Analysis_Tools/uniprot_pdb.csv')
    mapping_df = mapping_df.set_index('SP_PRIMARY')

    PDB_id_dict = map_all(ranked_list, cluster_file_dict,
                          mapping_df)
    return PDB_id_dict


def fill_in_pdb_info(in_table, PDB_id_dict):
    pdb_list = []
    for centroid_id in in_table['id']:
        if centroid_id in PDB_id_dict:
            pdb_list.append(PDB_id_dict[centroid_id])
        else:
            pdb_list.append(None)
    in_table['PDB'] = pdb_list


if __name__ == '__main__':
    # Parse command line arguments from stdin
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs=1, type=str,
                        help='Must be a distance matrix text file (i.e. '
                             'calculated with calc_distmtx.py).')
    parser.add_argument('outfile', nargs=1, type=str,
                        help='Must be a text file.')
    parser.add_argument('fig_filepath', nargs=1, type=str,
                        help='Filepath to output figure (.svg)')
    parser.add_argument('size_fasta', nargs=1, type=str,
                        help='Centroids FASTA file used to generate '
                             'alignment file is needed to get size of each '
                             'cluster. Must use -sizeout option with usearch '
                             'to generate such a FASTA file.')
    parser.add_argument('-csv', nargs=1, type=str)
    parser.add_argument('-title', nargs=1, type=str,
                        help='Title for plot - default is an empty string.')
    parser.add_argument('-rank', nargs=1, type=int,
                        help='Number of clusters to rank - default is all.')
    parser.add_argument('-size_weight', nargs=1, type=float,
                        help='Must be a number between 0 and 1. Indicates how '
                             'much size will be weighted in ranking. Values '
                             'that seem to \"look\" the best tend to fall '
                             'between 0.2 and 0.5. If zero, only distance '
                             'will be used to determine ranking.')
    parser.add_argument('-show', nargs=1, type=int,
                        help='Number of clusters to label - default is all.')
    parser.add_argument('-show_pdb', nargs=1, type=bool, default=False,
                        help='If true, highlight clusters for which PDB '
                             'structural information known for at least one '
                             'sequence within that cluster (May take a '
                             'while). Make sure to combine this with the '
                             '-cluster_file_dir option.')
    parser.add_argument('-cluster_file_dir', nargs=1, type=str,
                        help='If you want to show PDB information on plot, '
                             'make sure show_pdb is checked and provide the '
                             'directory of all the cluster fasta files ('
                             'generated in usearch via the -cluster option.)')
    parser.add_argument('-circle_size_multiplier', nargs=1, type=float,
                        default=1,
                        help='Can increase or decrease circle size by putting '
                             'in a multiplier - less than 1 decreases circle '
                             'size and greater than 1 increases circle size '
                             'proportionally.')

    args = parser.parse_args()

    # Initialize all lists, matrices, and dicts to be used in the script
    with open(vars(args)['infile'][0], 'r') as infile:
        # Get number of clusters by counting lines in script except last line
        # which starts with a tab
        num_clusters = 0
        print('reading infile')
        for line in infile:
            if line[0] != '\t':
                num_clusters += 1
        empty_matrix = [[0 for j in range(num_clusters)] for i in
                        range(num_clusters)]
        size_list = [0 for k in range(num_clusters)]
        size_sqrt_list = [0 for l in range(num_clusters)]
        labels = ['' for m in range(num_clusters)]

    # Get dict of cluster ids (keys) and their respective sizes (values)
    size_dict = get_size(vars(args)['size_fasta'][0])

    # Fill our empty matrix with values from the distance matrix text file.
    # Also generates 3 ordered lists: cluster size (the absolute area of the
    # cluster); the square root of that size (to account for the fact that r
    # is proportional to the square root of size for a circle, and r is used
    # to create the circles on the graph, not area); and labels, which are
    # just the cluster ids to be used in the graph labels.
    matrix, size_list, size_sqrt_list, labels = convert_distance_matrix(
        vars(args)['infile'][0], num_clusters, empty_matrix, size_dict)

    # Get a set of x, y coordinates from our set of distances in our distance
    #  matrix using mds via scikit learn.
    pos = transform_matrix(matrix)

    # Create pandas dataframe to easily view and manipulate the data as a table
    data_dict = {'id': labels, 'size': size_list, 'show': [False for n in
                                                           range(len(labels))],
                 'x': pos[:, 0], 'y': pos[:, 1], 'PDB': [None for n in range(
                                                         len(labels))],
                 'picked': [False for n in range(len(labels))]}
    table = pd.DataFrame(data_dict, index=labels)

    # Start our ranking at the largest cluster.
    largest_cluster = get_max_cluster(size_dict)
    table.set_value(largest_cluster, 'picked', True)
    pd.set_option('max_rows', 100)

    # Set how many clusters we want to choose (specified by '-rank' option).
    # Default to all clusters.
    if vars(args)['rank']:
        num_picked_clusters = vars(args)['rank'][0]
    else:
        num_picked_clusters = num_clusters

    if vars(args)['size_weight']:
        cluster_size_weight = vars(args)['size_weight'][0]
    else:
        cluster_size_weight = 0

    if vars(args)['show']:
        num_shown_clusters = vars(args)['show'][0]
    else:
        num_shown_clusters = 10

    if vars(args)['show_pdb']:
        show_pdb_info = True
        try:
            cluster_file_dir = vars(args)['cluster_file_dir'][0]
        except TypeError:
            raise Exception('If you want to show PDB information, you must '
                            'also use -cluster_file_dir option.')
    else:
        show_pdb_info = False
        cluster_file_dir = None

    # Generate ranked list from our table of cluster positions, the amount of
    #  clusters we would like to rank, and the starting cluster.
    ranked_list_centroids = get_ranked_list(table, num_picked_clusters,
                                            num_shown_clusters, largest_cluster,
                                            cluster_size_weight)

    # If desired, get PDB information (if any) and add to table
    if show_pdb_info:
        print('Extracting PDB information...')
        pdb_id_dict = get_pdb_id_dict(table['id'], cluster_file_dir)
        fill_in_pdb_info(table, pdb_id_dict)

    # Generate a list of colors based on which clusters were picked.
    color_list = ['' for i in range(num_clusters)]
    for i in range(num_clusters):
        if table.iloc[i].loc['PDB'] and table.iloc[i].loc['show']:
            color_list[i] = 'magenta'
        elif table.iloc[i].loc['PDB']:
            color_list[i] = 'green'
        elif table.iloc[i].loc['show']:
            color_list[i] = 'red'
        else:
            color_list[i] = 'blue'

    # Functionality to increase/decrease circle size
    if vars(args)['circle_size_multiplier']:
        circle_size_multiplier = float(vars(args)['circle_size_multiplier'][0])
        for index, size in enumerate(size_list):
            size_list[index] = circle_size_multiplier * size

    # Plot the clusters using x and y positions (pos) with circle sizes
    # taken from size_list and colors taken from color_list.
    plt.figure(figsize=(15, 15))
    plt.scatter(pos[:, 0], pos[:, 1], size_list, alpha=0.7, c=color_list)
    try:
        plt.title(vars(args)['title'][0], fontsize=60)
    except TypeError:
        print('No plot title chosen. Add a title with the -title option.\n')
    plt.grid(False)

    # Add ranking to list of cluster labels (form is 'cluster_id ranking_num')
    for label_num in range(len(labels)):
        for seq_id in ranked_list_centroids:
            if labels[label_num] == seq_id:
                labels[label_num] += (
                    '  ' + str(ranked_list_centroids.index(seq_id) + 1))

    # Add labels to plot
    for label, x, y in zip(labels, pos[:, 0], pos[:, 1]):
        # Only want labels for picked clusters
        if table.iloc[labels.index(label)].loc['show']:
            if table.iloc[labels.index(label)].loc['PDB']:
                plt.annotate(label, xy=(x, y), xytext=(10, 10), fontsize=22,
                             textcoords='offset points', ha='left', va='top',
                             bbox=dict(boxstyle='round', pad=0.5, fc='magenta',
                                       alpha=0.3))
            else:
                plt.annotate(label, xy=(x, y), xytext=(-10, 10), fontsize=22,
                             textcoords='offset points', ha='right',
                             va='bottom', bbox=dict(boxstyle='round,pad=0.5',
                                                    fc='red', alpha=0.3))
        elif table.iloc[labels.index(label)].loc['PDB']:
            plt.annotate(label, xy=(x, y),
                         xytext=(25, 25), textcoords='offset points',
                         ha='left', va='top', fontsize=22,
                         bbox=dict(boxstyle='round', pad=0.5, fc='green',
                                   alpha=0.3))

    # Add axis labels
    plt.xlabel('Relative Distance (arbitrary units)', size=50)
    plt.ylabel('Relative Distance (arbitrary units)', size=50)

    # Add legend to plot
    blue_patch = mpatches.Patch(color='blue', label='No PDB, not ranked in '
                                                    'top ' +
                                                    str(num_shown_clusters))
    red_patch = mpatches.Patch(color='red', label='No PDB, ranked in top ' +
                                                  str(num_shown_clusters))
    green_patch = mpatches.Patch(color='green', label='PDB structure exists, '
                                                      'not ranked in top ' +
                                                      str(num_shown_clusters))
    magenta_patch = mpatches.Patch(color='magenta',
                                   label='PDB structure exists, and ranked '
                                         'in top ' + str(num_shown_clusters))

    # plt.legend(handles=[blue_patch, red_patch, green_patch, magenta_patch],
    #           fontsize=30)

    # Generate ranked list text output file (to be used by find_reviewed.py)
    with open(vars(args)['outfile'][0], 'w') as outfile:
        for item in ranked_list_centroids:
            outfile.write(item + '\n')

    print(table.head())
    if vars(args)['csv']:
        table.to_csv(vars(args)['csv'][0])

    plt.savefig(vars(args)['fig_filepath'][0])
