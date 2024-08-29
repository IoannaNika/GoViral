import argparse
import sys
import pandas as pd
import igraph as ig
import random
from sklearn.metrics import homogeneity_score, completeness_score
import os

def get_sequence_ids(results):
    sequence_ids_1= results["Sequence_1_id"].unique()
    sequence_ids_2= results["Sequence_2_id"].unique()
    sequence_ids = set(sequence_ids_1).union(set(sequence_ids_2))
    return list(sequence_ids)

def make_output_file(output):
    # if output file exists, remove it
    if os.path.exists(output):
        os.remove(output)
    # open file to write the predictions
    file = open(output, "w")

    file.write("Community\tGenomic_regions\tSequence_id\tSequence\n")

    return file

def find_number_of_vertices(results):
    sequence_ids = get_sequence_ids(results)
    n = len(sequence_ids)
    print("Number of unique sequence ids: ", n)
    return n

def create_mappings_from_sequence_id_to_vertex_and_back(n, results):
    sequence_ids = get_sequence_ids(results)
    sequence_id_to_vertex = {}
    vertex_to_sequence_id = {}

    for i in range(n):
        sequence_id_to_vertex[sequence_ids[i]] = i
        vertex_to_sequence_id[i] = sequence_ids[i]

    return sequence_id_to_vertex, vertex_to_sequence_id

def from_sequence_id_to_vertex(sequence_id,  sequence_id_to_vertex):
        return sequence_id_to_vertex[sequence_id]

def from_vertex_to_sequence_id(vertex, vertex_to_sequence_id):
    return vertex_to_sequence_id[vertex]

def create_graph(results, sequence_id_to_vertex, n):
    edges = []
    weights = []
    for i in range(len(results)):
        if results.iloc[i]["Predicted_label"] == 1:
            edges.append((from_sequence_id_to_vertex(results.iloc[i]["Sequence_1_id"], sequence_id_to_vertex), from_sequence_id_to_vertex(results.iloc[i]["Sequence_2_id"], sequence_id_to_vertex)))
            weight = results.iloc[i]["Predicted_probability"]
            weights.append(weight)
        
    g = ig.Graph()
    g.add_vertices(n)
    g.add_edges(edges)
    g.es["weight"] = weights
    print("Number of predicted edges: ", len(edges))
    return g

def colors_for_communities(communities):

    colors_dict = {}
    for i in range(len(communities)):
        # get a  random color for each community
        color = (random.random(), random.random(), random.random())
        for j in range(len(communities[i])):
            colors_dict[communities[i][j]] =  color

    return colors_dict

def get_true_labels_and_predicted_labels_NCBI(predicted_communities, vertex_to_sequence_id):
    true_labels = []
    predicted_labels = []
    for i in range(len(predicted_communities)):
        community = predicted_communities[i]

        true_community_labels = []

        for vertex in community:
            sequence_id = vertex_to_sequence_id[vertex]
            true_community_labels.append(sequence_id.split("_")[0])
        
        # make list with same length as the true labels
        predicted_community_labels = [i for j in range(len(community))]

        true_labels.extend(true_community_labels)
        predicted_labels.extend(predicted_community_labels)

    return true_labels, predicted_labels


def write_predicted_communities_to_file(predicted_communities, vertex_to_sequence_id, genomic_region, output_file, results):
    for i in range(len(predicted_communities)):
        community = predicted_communities[i]
        for vertex in community:
            sequence_id = vertex_to_sequence_id[vertex]
            # find the sequence that corresponds to the sequence id
            sequence_ids_1= results["Sequence_1_id"].unique()
            if sequence_id in sequence_ids_1:
                sequence = results[results["Sequence_1_id"] == sequence_id]["Sequence_1"].values[0]
            else:
                sequence = results[results["Sequence_2_id"] == sequence_id]["Sequence_2"].values[0]

            output_file.write("{}\t{}\t{}\t{}\n".format(i, genomic_region, sequence_id, sequence))


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--results', dest = 'results', required=True, type=str, help="tsv file with predictions")
    parser.add_argument('--output', dest = 'output', required=True, type=str, help="output file with community predictions")
    args = parser.parse_args()
    
    results = pd.read_csv(args.results, sep='\t', header=0)
    output = args.output
    output_file = make_output_file(output)

    # take all possible genomic regions
    genomic_regions = results["Genomic_region"].unique()

    for genomic_region in genomic_regions:
        # do it for one genomic region
        results_per_gr = results[results["Genomic_region"] == genomic_region]
        # number of vertices is the number unique sequene ids, thus number of sequencing reads
        n = find_number_of_vertices(results_per_gr)

        # create a mapping from sequence ids to vertex indices
        sequence_id_to_vertex, vertex_to_sequence_id = create_mappings_from_sequence_id_to_vertex_and_back(n, results_per_gr)

        # create graph
        g = create_graph(results_per_gr, sequence_id_to_vertex, n)

        output_dir = "/".join(output.split("/")[:-1])
        # make folder for the graphs
        folder_path= output_dir + "/graphs/{}/".format(genomic_region)
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)

        # apply community detection
        # predicted_communities = g.community_multilevel(weights="weight")
        predicted_communities = g.community_leiden(weights = "weight", objective_function = "modularity", resolution = 0.75)
        print("Number of predicted communities: ", len(predicted_communities))

        # define colors for each community in the graph
        colors_dict = colors_for_communities(predicted_communities) 
        # assign colors to vertices
        g.vs["color"] = [colors_dict[i] for i in range(n)]
        # assign the labels to the vertices
        g.vs["label"] = [vertex_to_sequence_id[i].split("_")[0] for i in range(n)]
        # visualize the graph with the communities with labels, and position the labels next to the vertices, not in the middle of the vertices 
        ig.plot(g, folder_path + "graph_predicted.pdf", vertex_size=10, vertex_color=g.vs["color"], edge_width=0.3, edge_color="black", vertex_label_dist=2, vertex_label_size=0, margin=50, rotation=30)
        
        # write the predicted communities to a file
        write_predicted_communities_to_file(predicted_communities, vertex_to_sequence_id, genomic_region, output_file, results_per_gr)

if __name__ == "__main__":
    sys.exit(main())