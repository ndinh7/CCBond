# modified from Jon's script to run cluster analysis

# User Input

# path to cluster analysis tool folder
CA_TOOLPATH=../utils/Cluster_Analysis_Tools

# path to project folder
CA_PROJECT=./cluster/20190219_substratecluster_01

# name of project
GROUP='substrate'
OUT_PREFIX=${GROUP}
PLOT_TITLE="ZINC Alpha-Keto Acids Clustering"

# How many clusters you want to rank (the more, the longer it will take)
RANK=50

# How heavily you want to weight cluster size when ranking clusters. The higher the value, the more
# likely big clusters will appear near the top of the list, even if they are close together.
# Values around 0.1 - 0.3 seem to work well.
SIZE_WEIGHT=0.3

# How many clusters you want to show on the plot (e.g. SHOW=5 will color in the top 5 clusters)
# I wouldn't show a lot because the plot can get really crowded if you do.
SHOW=0

# We then use the distance matrix to plot the clusters and rank them in the order we should test
# them in.
python $CA_TOOLPATH/rank_and_plot_distmtx.py \
	-csv ${CA_PROJECT}/${OUT_PREFIX}_table.csv \
	-title "$PLOT_TITLE" \
	-rank $RANK \
	-size_weight $SIZE_WEIGHT \
	-show $SHOW \
	-cluster_file_dir ${CA_PROJECT}/clusters_${OUT_PREFIX} \
	-circle_size_multiplier 10 \
	${CA_PROJECT}/${OUT_PREFIX}_distmtx.txt \
	${CA_PROJECT}/${OUT_PREFIX}_ranked.txt \
	${CA_PROJECT}/${OUT_PREFIX}_cluster_plot \
	${CA_PROJECT}/${OUT_PREFIX}_centroids.txt
