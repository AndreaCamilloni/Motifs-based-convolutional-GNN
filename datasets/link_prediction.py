from __future__ import division
from math import floor
import math
import os
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import scipy.sparse as sp
import torch
from torch.utils.data import Dataset
from torchvision import models
import torch.nn as nn

from PIL import Image


import networkx as nx
from itertools import combinations
 

import utils
# try:
#     from neoConnector import all_cells_with_n_hops_in_area, get_all_edges
# except ImportError:
#     from .neoConnector import all_cells_with_n_hops_in_area, get_all_edges

np.random.seed(0)

class_map = {'inflammatory': 0, 'lymphocyte' : 1, 'fibroblast and endothelial': 2,
               'epithelial': 3, 'apoptosis / civiatte body': 4}


class KIGraphDataset2(Dataset):

    def __init__(self, path, mode='train',
                 num_layers=2,
                 data_split=[0.8, 0.2], add_self_edges=False):
        """
        Parameters
        ----------
        path : list
            List with filename, coordinates and path to annotation. For example, ['P7_HE_Default_Extended_1_1', (0, 2000, 0, 2000), 'datasets/annotations/P7_annotated/P7_HE_Default_Extended_1_1.txt']
        mode : str
            One of train, val or test. Default: train.
        num_layers : int
            Number of layers in the computation graph. Default: 2.
        data_split: list
            Fraction of edges to use for graph construction / train / val / test. Default: [0.85, 0.08, 0.02, 0.03].
        """
        super().__init__()

        self.path = path
        self.mode = mode
        self.num_layers = num_layers
        self.data_split = data_split

        print('--------------------------------')
        print('Reading edge dataset from {}'.format(self.path[0]))

        ########## MINE ###########
        # Cells, distance_close_to_edges
        edge_path = path[1]
        node_path = path[2]

        # with glob
        edges = pd.read_csv(edge_path)
        nodes = pd.read_csv(node_path)

        if add_self_edges:
            for i in range(len(nodes)):
                new_row = {'source': i, 'target': i, 'type': 0, 'distance': 0}
                # append row to the dataframe
                edges = edges.append(new_row, ignore_index=True)

        edges_crossing = edges.copy()
        edges_crossing = edges_crossing[edges_crossing["type"] == 1]

        edges['type'] = edges['type'].replace(1, 0)

        col_row_len = len(nodes['id'])
        distances_close_to_edges = pd.DataFrame(0, index=np.arange(col_row_len), columns=np.arange(col_row_len))

        for i, row in edges.iterrows():
            source = row['source']
            target = row['target']
            distance = row['distance']
            distances_close_to_edges[source][target] = distance
            distances_close_to_edges[target][source] = distance

        distances_close_to_edges = np.array(distances_close_to_edges)

        # coords
        coords = nodes[["x", "y"]].to_numpy()

        # all_labels_cell_types
        nodes["gt"].replace({'inflammatory': 0, 'lymphocyte': 1, 'fibroblast and endothelial': 2, 'epithelial': 3}, inplace=True) # hover-net
                     

        # nuclei features
        #nuclei_feat = nodes[["area", "perim"]].to_numpy()

        all_labels_cell_types = nodes["gt"].to_numpy()
        nodes_with_types_zero_one = nodes.copy()
        nodes_with_types_prob = nodes.copy()
        for i, row in nodes.iterrows():
            ##The iloc positions are based on the csv cell positions
            # if row['gt'] == 'null':
            nodes_with_types_zero_one.iloc[i, 3] = 1 if row['gt'] == 2 else 0
            nodes_with_types_zero_one.iloc[i, 4] = 1 if row['gt'] == 0 else 0
            nodes_with_types_zero_one.iloc[i, 1] = 1 if row['gt'] == 1 else 0
            nodes_with_types_zero_one.iloc[i, 2] = 1 if row['gt'] == 3 else 0
    
        print('<<<<<<<<<<<<<', nodes_with_types_zero_one)

        
        # cell_types_scores
        cell_types_scores = nodes_with_types_zero_one[['inf', 'lym', 'fib', 'epi']] #One-hot encoding of GT data

        cell_types_scores = cell_types_scores.to_numpy()
        print(cell_types_scores.shape)

        # adjacency_matrix_close_to_edges
        adjacency_matrix_close_to_edges = np.copy(distances_close_to_edges)
        adjacency_matrix_close_to_edges[adjacency_matrix_close_to_edges != 0] = 1

        # edge_list_close_to_edge
        edge_list_close_to_edge = edges[["source", "target"]]
        edge_list_close_to_edge = edge_list_close_to_edge.to_numpy()

        # edge_list_crossing_edges TODO
        edge_list_crossing_edges = edges_crossing.to_numpy()

        self.am_close_to_edges_including_distances = distances_close_to_edges
        self.classes = all_labels_cell_types
        self.class_scores = cell_types_scores
        self.coords = coords

        print('Finished reading data.')

        print('Setting up graph.')
        vertex_id = {j: i for (i, j) in enumerate(range(len(coords)))}

        idxs = [floor(v * edge_list_crossing_edges.shape[0]) for v in np.cumsum(data_split)]

        edges_t, pos_examples_crossing_edges = edge_list_close_to_edge, edge_list_crossing_edges

        edges_t[:, :2] = np.array([vertex_id[u] for u in edges_t[:, :2].flatten()]).reshape(edges_t[:, :2].shape)
        edges_t_no_duplicates = np.unique(edges_t[:, :2], axis=0)  # Filter duplicate edges

        self.nodes_count = len(vertex_id)  # Count vertices
        self.edges_count = edges_t_no_duplicates.shape[0]  # Count edges

        adjacency_matrix_close_to_edges = sp.coo_matrix(
            (np.ones(self.edges_count), (edges_t_no_duplicates[:, 0], edges_t_no_duplicates[:, 1])),
            shape=(self.nodes_count, self.nodes_count),
            dtype=np.float32)

        self.adjacency_matrix_close_to_edges_as_coo_to_lil = adjacency_matrix_close_to_edges.tolil()

        self.node_neighbors = self.adjacency_matrix_close_to_edges_as_coo_to_lil.rows  # Neighbors

        cell_density = nodes['Cell_density'].to_numpy() 
        cell_density = np.array(cell_density)
        cell_density = cell_density.astype(np.float64)

        cell_entropy = nodes['Node_Entropy'].to_numpy()
        cell_entropy = np.array(cell_entropy)
        cell_entropy = cell_entropy.astype(np.float64)

        #mean_neigh_entropy = nodes['Mean_Neighborhood_Entropy'].to_numpy() 
        #mean_neigh_entropy = np.array(mean_neigh_entropy)
        #mean_neigh_entropy = mean_neigh_entropy.astype(np.float64)

        #graph_node_features = np.concatenate((cell_types_scores, cell_density[:,None], cell_entropy[:,None]), axis=1 )  ### Concatenate all features 
        #self.features = torch.from_numpy(graph_node_features).float()  # Cell features 

        self.features = torch.from_numpy(cell_types_scores).float()  # Cell features   


        print('self.features.shape:', self.features.shape)
        # [2] end

        print('Finished setting up graph.')

        print('Setting up examples.')

        if len(pos_examples_crossing_edges) > 0:
            pos_examples_crossing_edges = pos_examples_crossing_edges[:, :2]
            pos_examples_crossing_edges = np.unique(pos_examples_crossing_edges, axis=0)

        # Generate negative examples not in cell edges crossing path
        neg_examples_close_to_edges = []
        cur = 0
        n_count, _choice = self.nodes_count, np.random.choice
        neg_seen = set(tuple(e[:2]) for e in edge_list_crossing_edges)  # Dont sample positive edges
        adj_tuple = set(tuple(e[:2]) for e in edge_list_close_to_edge)  # List all edges

        if self.mode != 'train':  # Add all edges except positive edges if validation/test
            print("self.mode != 'train'")
            for example in edge_list_close_to_edge:
                if (example[0], example[1]) in neg_seen:
                    continue
                neg_examples_close_to_edges.append(example)
            neg_examples_close_to_edges = np.array(neg_examples_close_to_edges, dtype=np.int64)
        else:  # Undersample negative samples from adjacency edges not in positive
            num_neg_examples = pos_examples_crossing_edges.shape[0]
            while cur < num_neg_examples:
                u, v = _choice(n_count, 2, replace=False)
                if (u, v) in neg_seen or (u, v) not in adj_tuple:
                    continue
                cur += 1
                neg_examples_close_to_edges.append([u, v])
            neg_examples_close_to_edges = np.array(neg_examples_close_to_edges, dtype=np.int64)

        x = np.vstack((pos_examples_crossing_edges, neg_examples_close_to_edges))
        y = np.concatenate((np.ones(pos_examples_crossing_edges.shape[0]),
                            np.zeros(neg_examples_close_to_edges.shape[0])))
        perm = np.random.permutation(x.shape[0])
        x, y = x[perm, :], y[perm]  # ERROR HERE -> IndexError: too many indices for array: array is 1-dimensional,
        # but 2 were indexed
        x, y = torch.from_numpy(x).long(), torch.from_numpy(y).long()
        self.x, self.y = x, y

        print('Finished setting up examples.')

        print('Dataset properties:')
        print('Mode: {}'.format(self.mode))
        print('Number of vertices: {}'.format(self.nodes_count))
        print('Number of edges: {}'.format(self.edges_count))
        print('Number of positive/negative datapoints: {}/{}'.format(pos_examples_crossing_edges.shape[0],
                                                                     neg_examples_close_to_edges.shape[0]))
        print('Number of examples/datapoints: {}'.format(self.x.shape[0]))
        print('--------------------------------')

    def __len__(self):
        return len(self.x)

    def __getitem__(self, idx):
        return self.x[idx], self.y[idx]

    def get_coords_and_class(self):
        return self.coords, self.classes

    def _form_computation_graph(self, idx):
        """
        Parameters
        ----------
        idx : int or list
            Indices of the node for which the forward pass needs to be computed.
        Returns
        -------
        node_layers : list of numpy array
            node_layers[i] is an array of the nodes in the ith layer of the
            computation graph.
        mappings : list of dictionary
            mappings[i] is a dictionary mapping node v (labelled 0 to |V|-1)
            in node_layers[i] to its position in node_layers[i]. For example,
            if node_layers[i] = [2,5], then mappings[i][2] = 0 and
            mappings[i][5] = 1.
        """
        _list, _set = list, set
        if type(idx) is int:
            node_layers = [np.array([idx], dtype=np.int64)]
        elif type(idx) is list:
            node_layers = [np.array(idx, dtype=np.int64)]

        for _ in range(self.num_layers):
            prev = node_layers[-1]
            arr = [node for node in prev]
            arr.extend([e for node in arr for e in self.node_neighbors[node]])  # add neighbors to graph
            arr = np.array(_list(_set(arr)), dtype=np.int64)
            node_layers.append(arr)
        node_layers.reverse()

        mappings = [{j: i for (i, j) in enumerate(arr)} for arr in node_layers]

        return node_layers, mappings

    def collate_wrapper(self, batch):
        """
        Parameters
        ----------
        batch : list
            A list of examples from this dataset. An example is (edge, label).
        Returns
        -------
        edges : numpy array
            The edges in the batch.
        features : torch.FloatTensor
            An (n' x input_dim) tensor of input node features.
        node_layers : list of numpy array
            node_layers[i] is an array of the nodes in the ith layer of the
            computation graph.
        mappings : list of dictionary
            mappings[i] is a dictionary mapping node v (labelled 0 to |V|-1)
            in node_layers[i] to its position in node_layers[i]. For example,
            if node_layers[i] = [2,5], then mappings[i][2] = 0 and
            mappings[i][5] = 1.
        rows : numpy array
            Each row is the list of neighbors of nodes in node_layers[0].
        labels : torch.LongTensor
            Labels (1 or 0) for the edges in the batch.
        """
        idx = list(set([v.item() for sample in batch for v in sample[0][:2]]))

        node_layers, mappings = self._form_computation_graph(idx)

        rows = self.node_neighbors[node_layers[0]]
        features = self.features[node_layers[0], :]

        dist = torch.from_numpy(self.am_close_to_edges_including_distances)
        labels = torch.FloatTensor([sample[1] for sample in batch])
        edges = np.array([sample[0].numpy() for sample in batch])
        edges = np.array([mappings[-1][v] for v in edges.flatten()]).reshape(edges.shape)

        # TODO: Pin memory. Change type of node_layers, mappings and rows to
        # tensor?

        return edges, features, node_layers, mappings, rows, labels, dist

    def get_dims(self):
        print("self.features.shape: {}".format(self.features.shape))
        print("input_dims (input dimension) -> self.features.shape[1] = {}".format(self.features.shape[1]))
        return self.features.shape[1], 1

    def parse_points(self, fname):
        with open(fname, 'r') as f:
            lines = f.readlines()
        lines = [line[:-1].split(',') for line in lines]  # Remove \n from line
        return lines

def adj_to_edge(adj):
    edges = []
    for i in range(len(adj)):
        edges += ([[i,index] for index, element in enumerate(adj[i]) if element == 1])

    return edges

def get_intersections(points, coords, adj):
    # Loop through cells
    intersections = []
    count = 0
    for i in range(len(coords)):
        # Get ids of all neighbors
        nbrs = [index for index, element in enumerate(adj[i]) if element == 1]
        for j in range(len(nbrs)):
            passed = False
            for k in range(len(points)-2):
                if len(points[k]) == 2 and len(points[k+1]) == 2:
                    L1 = line(coords[i], coords[nbrs[j]]) # Line between node and neighbor
                    L2 = line([int(float(point)) for point in points[k]], [int(float(point)) for point in points[k+1]]) # Line between two points of path
                    inter = intersection(L1, L2) # Get x-coordinate for intersection or False if none
                    if inter != False:
                        if ( (inter > max( min(coords[i][0],coords[nbrs[j]][0]), min(int(float(points[k][0])),int(float(points[k+1][0]))) )) and
                            (inter < min( max(coords[i][0],coords[nbrs[j]][0]), max(int(float(points[k][0])),int(float(points[k+1][0]))) )) ): # If intersection is inside line segments
                            intersections.append([i, nbrs[j]])
                            passed = True
                            break
            #if passed == False: # If no intersections between cell and neighbor
            #    intersections.append([i, nbrs[j], 0])
            #
    return intersections

def line(p1, p2):
    A = (p1[1] - p2[1])
    B = (p2[0] - p1[0])
    C = (p1[0]*p2[1] - p2[0]*p1[1])
    return A, B, -C

def intersection(L1, L2):
    D  = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    if D != 0:
        x = Dx / D
        return x
    else:
        return False


# dataset for Graph Convolutional Neural Networks
class KIGraphDatasetGCN(Dataset):

    def __init__(self, path, mode='train',
                 num_layers=2,
                 data_split=[0.8, 0.2], add_self_edges=False):
        """
        Parameters
        ----------
        path : list
            List with filename, coordinates and path to annotation. For example, ['P7_HE_Default_Extended_1_1', (0, 2000, 0, 2000), 'datasets/annotations/P7_annotated/P7_HE_Default_Extended_1_1.txt']
        mode : str
            One of train, val or test. Default: train.
        num_layers : int
            Number of layers in the computation graph. Default: 2.
        data_split: list
            Fraction of edges to use for graph construction / train / val / test. Default: [0.85, 0.08, 0.02, 0.03].
        """
        super().__init__()

        self.path = path
        self.mode = mode
        self.num_layers = num_layers
        self.data_split = data_split

        print('--------------------------------')
        print('Reading edge dataset from {}'.format(self.path[0]))

        ########## MINE ###########
        # Cells, distance_close_to_edges
        edge_path = path[1]
        node_path = path[2]

        # with glob
        edges = pd.read_csv(edge_path)
        nodes = pd.read_csv(node_path)

        if add_self_edges:
            for i in range(len(nodes)):
                new_row = {'source': i, 'target': i, 'type': 0, 'distance': 0}
                # append row to the dataframe
                edges = edges.append(new_row, ignore_index=True)

        edges_crossing = edges.copy()
        edges_crossing = edges_crossing[edges_crossing["type"] == 1]

        edges['type'] = edges['type'].replace(1, 0)

        col_row_len = len(nodes['id'])
        distances_close_to_edges = pd.DataFrame(0, index=np.arange(col_row_len), columns=np.arange(col_row_len))
        delta_entropy_edges = pd.DataFrame(0, index=np.arange(col_row_len), columns=np.arange(col_row_len))
        neighborhood_similarity_edges = pd.DataFrame(0, index=np.arange(col_row_len), columns=np.arange(col_row_len))

        for i, row in edges.iterrows():
            source = row['source']
            target = row['target']
            
            distance = float(row['distance'])
           
            delta_entropy = float(row['Delta_Entropy'])
            sorenson_neigh_similarity = float(row['Sorenson_Similarity'])
        
            distances_close_to_edges[source][target] = distance
            distances_close_to_edges[target][source] = distance

            delta_entropy_edges[source][target] = delta_entropy
            delta_entropy_edges[target][source] = delta_entropy

            neighborhood_similarity_edges[source][target] = sorenson_neigh_similarity
            neighborhood_similarity_edges[target][source] = sorenson_neigh_similarity


        distances_close_to_edges = np.array(distances_close_to_edges)
        delta_entropy_edges = np.array(delta_entropy_edges)
        neighborhood_similarity_edges = np.array(neighborhood_similarity_edges)

        # coords
        coords = nodes[["x", "y"]].to_numpy()

        # process neighborhood densities
        density_types = ["Cell_density"]
        #entropy_types = ["Node_Entropy"]

        densities = nodes[density_types].to_numpy()
        edge_density = np.zeros((col_row_len, col_row_len))
        edge_densities = np.empty((0, col_row_len, col_row_len))

        
        for i in range(len(density_types)):
            for _, row in edges.iterrows():
                source = int(row['source'])
                target = int(row['target'])

                edge_density[source][target] = float(densities[:, i][target]) - float(densities[:, i][source])
                edge_density[target][source] = float(densities[:, i][source]) - float(densities[:, i][target])

            edge_densities = np.append(edge_densities, edge_density.reshape(-1, col_row_len, col_row_len), axis=0)

        #print('*************')
        #print('Edge_density Shape : ' + str(edge_densities.shape))
        
        distances_close_to_edges = distances_close_to_edges.reshape(-1, col_row_len, col_row_len)
        delta_entropy_edges = delta_entropy_edges.reshape(-1, col_row_len, col_row_len)
        neighborhood_similarity_edges = neighborhood_similarity_edges.reshape(-1, col_row_len, col_row_len)

        #print('Edge_entropy Shape : ' + str(delta_entropy_edges.shape))
        #print('Edge_distance Shape : ' + str(distances_close_to_edges.shape))
        #print('Neighborhood Similarity Shape : ' + str(neighborhood_similarity_edges.shape))


        #edge_features = np.concatenate((edge_densities, delta_entropy_edges, neighborhood_similarity_edges, distances_close_to_edges), axis=0) #ablation1234
        edge_features = np.concatenate((edge_densities, delta_entropy_edges, distances_close_to_edges), axis=0)                                 #ablation 123
        #edge_features = np.concatenate((edge_densities, distances_close_to_edges), axis=0)                                                      #ablation 12
        #edge_features = np.concatenate((distances_close_to_edges), axis=0)                                                                      #ablation 1
        #edge_features = delta_entropy_edges
        #print(edge_features)
        # self.edge_features = utils.normalize_edge_feature_doubly_stochastic(edge_features) ### not to be used
        self.edge_features = utils.normalize_edge_features_rows(edge_features) ### Use it to normalise the edge features
        #self.edge_features = edge_features  ### Use only if not using the normalization feature above


        ## To DO 

        # Change utils.normalize_edge_features_rows to log function to the base e

        #####

        self.channel = edge_features.shape[0]

        self.dist = utils.normalize_edge_features_rows(distances_close_to_edges.reshape(-1, col_row_len, col_row_len))


        # all_labels_cell_types
        nodes["gt"].replace({'inflammatory': 0, 'lymphocyte': 1, 'fibroblast and endothelial': 2, 'epithelial': 3}, inplace=True) # hover-net
                     

        # nuclei features
        #nuclei_feat = nodes[["area", "perim"]].to_numpy()

        all_labels_cell_types = nodes["gt"].to_numpy()

        nodes_with_types_zero_one = nodes.copy()
        nodes_with_types_prob = nodes.copy()
        for i, row in nodes.iterrows():
            ##The iloc positions are based on the csv cell positions
            # if row['gt'] == 'null':
            nodes_with_types_zero_one.iloc[i, 3] = 1 if row['gt'] == 2 else 0
            nodes_with_types_zero_one.iloc[i, 4] = 1 if row['gt'] == 0 else 0
            nodes_with_types_zero_one.iloc[i, 1] = 1 if row['gt'] == 1 else 0
            nodes_with_types_zero_one.iloc[i, 2] = 1 if row['gt'] == 3 else 0
    
        #print('<<<<<<<<<<<<<', nodes_with_types_zero_one)

        
        # cell_types_scores
        cell_types_scores = nodes_with_types_zero_one[['inf', 'lym', 'fib', 'epi']] #One-hot encoding of GT data

        cell_types_scores = cell_types_scores.to_numpy()
        #print(cell_types_scores.shape)

        # adjacency_matrix_close_to_edges
        adjacency_matrix_close_to_edges = np.copy(distances_close_to_edges)
        adjacency_matrix_close_to_edges[adjacency_matrix_close_to_edges != 0] = 1
        self.adj = adjacency_matrix_close_to_edges

        # edge_list_close_to_edge
        edge_list_close_to_edge = edges[["source", "target"]]
        edge_list_close_to_edge = edge_list_close_to_edge.to_numpy()

        # edge_list_crossing_edges
        edge_list_crossing_edges = edges_crossing.to_numpy()

        self.am_close_to_edges_including_distances = distances_close_to_edges
        self.classes = all_labels_cell_types
        self.class_scores = cell_types_scores
        self.coords = coords

        print('Finished reading data.')

        print('Setting up graph.')
        vertex_id = {j: i for (i, j) in enumerate(range(len(coords)))}

        idxs = [floor(v * edge_list_crossing_edges.shape[0]) for v in np.cumsum(data_split)]

        edges_t, pos_examples_crossing_edges = edge_list_close_to_edge, edge_list_crossing_edges

        edges_t[:, :2] = np.array([vertex_id[u] for u in edges_t[:, :2].flatten()]).reshape(edges_t[:, :2].shape)
        edges_t_no_duplicates = np.unique(edges_t[:, :2], axis=0)  # Filter duplicate edges

        self.nodes_count = len(vertex_id)  # Count vertices
        self.edges_count = edges_t_no_duplicates.shape[0]  # Count edges

        adjacency_matrix_close_to_edges = sp.coo_matrix(
            (np.ones(self.edges_count), (edges_t_no_duplicates[:, 0], edges_t_no_duplicates[:, 1])),
            shape=(self.nodes_count, self.nodes_count),
            dtype=np.float32)

        self.adjacency_matrix_close_to_edges_as_coo_to_lil = adjacency_matrix_close_to_edges.tolil()

        self.node_neighbors = self.adjacency_matrix_close_to_edges_as_coo_to_lil.rows  # Neighbors

        #### Code to add node features ### 

        cell_density = nodes['Cell_density'].to_numpy() 
        cell_density = np.array(cell_density)
        cell_density = cell_density.astype(np.float64)

        cell_entropy = nodes['Node_Entropy'].to_numpy()
        cell_entropy = np.array(cell_entropy)
        cell_entropy = cell_entropy.astype(np.float64)

        #mean_neigh_entropy = nodes['Mean_Neighborhood_Entropy'].to_numpy() 
        #mean_neigh_entropy = np.array(mean_neigh_entropy)
        #mean_neigh_entropy = mean_neigh_entropy.astype(np.float64)

        #graph_node_features = np.concatenate((cell_types_scores, cell_density[:,None], cell_entropy[:,None]), axis=1 )  ### Concatenate all features 
        #self.features = torch.from_numpy(graph_node_features).float()  # Cell features 

        ###### Code to add node features ends here ##### 

        ### Use this self feature if only one-hot embedding is required as node feature set
        self.features = torch.from_numpy(cell_types_scores).float()  # Cell features with just one-hot encoding 

        print('self.features.shape:', self.features.shape)
        # [2] end

        print('Finished setting up graph.')

        print('Setting up examples.')

        if len(pos_examples_crossing_edges) > 0:
            pos_examples_crossing_edges = pos_examples_crossing_edges[:, :2]
            pos_examples_crossing_edges = np.unique(pos_examples_crossing_edges, axis=0)

        # Generate negative examples not in cell edges crossing path
        neg_examples_close_to_edges = []
        cur = 0
        n_count, _choice = self.nodes_count, np.random.choice
        neg_seen = set(tuple(e[:2]) for e in edge_list_crossing_edges)  # Dont sample positive edges
        adj_tuple = set(tuple(e[:2]) for e in edge_list_close_to_edge)  # List all edges

        if self.mode != 'train':  # Add all edges except positive edges if validation/test
            print("self.mode != 'train'")
            for example in edge_list_close_to_edge:
                if (example[0], example[1]) in neg_seen:
                    continue
                neg_examples_close_to_edges.append(example)
            neg_examples_close_to_edges = np.array(neg_examples_close_to_edges, dtype=np.int64)
        else:  # Undersample negative samples from adjacency edges not in positive

            num_neg_examples = pos_examples_crossing_edges.shape[0] # for BCE loss function 

            #If using Focal Loss function use the next line else comment it and uncomment the line above for BCE loss 
            #num_neg_examples = int(pos_examples_crossing_edges.shape[0]) # Increasing the size of neg_samples for focal loss function as it can handle class imbalance 
            
            while cur < num_neg_examples:
                u, v = _choice(n_count, 2, replace=False)
                if (u, v) in neg_seen or (u, v) not in adj_tuple:
                    continue
                cur += 1
                neg_examples_close_to_edges.append([u, v])
            neg_examples_close_to_edges = np.array(neg_examples_close_to_edges, dtype=np.int64)

        x = np.vstack((pos_examples_crossing_edges, neg_examples_close_to_edges))
        y = np.concatenate((np.ones(pos_examples_crossing_edges.shape[0]),
                            np.zeros(neg_examples_close_to_edges.shape[0])))
        perm = np.random.permutation(x.shape[0])
        x, y = x[perm, :], y[perm]  # ERROR HERE -> IndexError: too many indices for array: array is 1-dimensional,
        # but 2 were indexed
        x, y = torch.from_numpy(x).long(), torch.from_numpy(y).long()
        self.x, self.y = x, y

        print('Finished setting up examples.')

        print('Dataset properties:')
        print('Mode: {}'.format(self.mode))
        print('Number of vertices: {}'.format(self.nodes_count))
        print('Number of edges: {}'.format(self.edges_count))
        print('Number of positive/negative datapoints: {}/{}'.format(pos_examples_crossing_edges.shape[0],
                                                                     neg_examples_close_to_edges.shape[0]))
        print('Number of examples/datapoints: {}'.format(self.x.shape[0]))
        print('--------------------------------')

    def __len__(self):
        return len(self.x)

    def __getitem__(self, idx):
        return self.x[idx], self.y[idx]

    def get_coords_and_class(self):
        return self.coords, self.classes

    def _form_computation_graph(self, idx):
        """
        Parameters
        ----------
        idx : int or list
            Indices of the node for which the forward pass needs to be computed.
        Returns
        -------
        node_layers : list of numpy array
            node_layers[i] is an array of the nodes in the ith layer of the
            computation graph.
        mappings : list of dictionary
            mappings[i] is a dictionary mapping node v (labelled 0 to |V|-1)
            in node_layers[i] to its position in node_layers[i]. For example,
            if node_layers[i] = [2,5], then mappings[i][2] = 0 and
            mappings[i][5] = 1.
        """
        _list, _set = list, set
        if type(idx) is int:
            node_layers = [np.array([idx], dtype=np.int64)]
        elif type(idx) is list:
            node_layers = [np.array(idx, dtype=np.int64)]

        for _ in range(self.num_layers):
            prev = node_layers[-1]
            arr = [node for node in prev]
            arr.extend([e for node in arr for e in self.node_neighbors[node]])  # add neighbors to graph
            arr = np.array(_list(_set(arr)), dtype=np.int64)
            node_layers.append(arr)
        node_layers.reverse()

        mappings = [{j: i for (i, j) in enumerate(arr)} for arr in node_layers]

        return node_layers, mappings

    def collate_wrapper(self, batch):
        """
        Parameters
        ----------
        batch : list
            A list of examples from this dataset. An example is (edge, label).
        Returns
        -------
        adj : torch.Tensor
            adjacency matrix of entire graph
        features : torch.FloatTensor
            A (n' x input_dim) tensor of input node features.
        edge_features : torch.FloatTensor
            A 3d tensor of edge features.
        edges : numpy array
            The edges in the batch.
        labels : torch.LongTensor
            Labels (1 or 0) for the edges in the batch.
        dist : torch.Tensor
            A distance matrix
        """
        adj = torch.from_numpy(self.adj).float()

        features = self.features
        edge_features = torch.from_numpy(self.edge_features).float()
        edges = np.array([sample[0].numpy() for sample in batch])
        labels = torch.FloatTensor([sample[1] for sample in batch])
        dist = torch.from_numpy(self.dist)

        return adj, features, edge_features, edges, labels, dist

    def get_dims(self):
        print("self.features.shape: {}".format(self.features.shape))
        print("input_dims (input dimension) -> self.features.shape[1] = {}".format(self.features.shape[1]))
        return self.features.shape[1], 1

    def get_channel(self):
        return self.channel

    def parse_points(self, fname):
        with open(fname, 'r') as f:
            lines = f.readlines()
        lines = [line[:-1].split(',') for line in lines]  # Remove \n from line
        return lines

def adj_to_edge(adj):
    edges = []
    for i in range(len(adj)):
        edges += ([[i,index] for index, element in enumerate(adj[i]) if element == 1])

    return edges

def get_intersections(points, coords, adj):
    # Loop through cells
    intersections = []
    count = 0
    for i in range(len(coords)):
        # Get ids of all neighbors
        nbrs = [index for index, element in enumerate(adj[i]) if element == 1]
        for j in range(len(nbrs)):
            passed = False
            for k in range(len(points)-2):
                if len(points[k]) == 2 and len(points[k+1]) == 2:
                    L1 = line(coords[i], coords[nbrs[j]]) # Line between node and neighbor
                    L2 = line([int(float(point)) for point in points[k]], [int(float(point)) for point in points[k+1]]) # Line between two points of path
                    inter = intersection(L1, L2) # Get x-coordinate for intersection or False if none
                    if inter != False:
                        if ( (inter > max( min(coords[i][0],coords[nbrs[j]][0]), min(int(float(points[k][0])),int(float(points[k+1][0]))) )) and
                            (inter < min( max(coords[i][0],coords[nbrs[j]][0]), max(int(float(points[k][0])),int(float(points[k+1][0]))) )) ): # If intersection is inside line segments
                            intersections.append([i, nbrs[j]])
                            passed = True
                            break
            #if passed == False: # If no intersections between cell and neighbor
            #    intersections.append([i, nbrs[j], 0])
            #
    return intersections

def line(p1, p2):
    A = (p1[1] - p2[1])
    B = (p2[0] - p1[0])
    C = (p1[0]*p2[1] - p2[0]*p1[1])
    return A, B, -C

def intersection(L1, L2):
    D  = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    if D != 0:
        x = Dx / D
        return x
    else:
        return False


 

# Function to find triangles in a graph using nodes and edges DataFrames
def find_triangles_old(nodes_df, edges_df):
    # Create a NetworkX graph from the nodes and edges DataFrames
    G = nx.Graph()
    G.add_nodes_from(nodes_df['id'])
    G.add_edges_from(edges_df[['source', 'target']].values)

    # Find triangles
    triangles = []
    for node in G.nodes:
        neighbors = list(G.neighbors(node))
        for u, v in combinations(neighbors, 2):
            if G.has_edge(u, v):
                triangles.append(sorted([node, u, v]))
    
    # return a set of unique triangles
    return set(map(tuple, triangles))
    #return [list(x) for x in set(tuple(x) for x in triangles)]

def find_triangles_by_edge(nodes_df, edges_df, DROP_EDGES=False):
    # Create a NetworkX graph from the nodes and edges DataFrames
    G = nx.Graph()
    G.add_nodes_from(nodes_df['id'])
    G.add_edges_from(edges_df[['source', 'target']].values)

    # Find triangles
    triangles = set()
    triangles_dict = dict()
    for (u,v) in G.edges:
        tmp = set()
        neighbors_u = set(G.neighbors(u))
        neighbors_v = set(G.neighbors(v))

        vertices = sorted(neighbors_u.intersection(neighbors_v))

        if len(vertices) == 1: 
            if DROP_EDGES:
                print("Dropping edge ({},{})".format(u,v))
                # drop edge from df
                edges_df = edges_df.drop(edges_df[(edges_df['source'] == u) & (edges_df['target'] == v)].index)
                edges_df = edges_df.drop(edges_df[(edges_df['source'] == v) & (edges_df['target'] == u)].index)
                continue
            else:
                #triangles_dict[frozenset((u,v))] = 
                t1 = frozenset((u,v,vertices[0]))
                tmp.add(t1)
                triangles.add(t1)
                # MAYBE WE CAN DEAL WITH PADDING LATER IN THE CODE
                #PADDING: add random vertex = u,v or z    
                rand = np.random.rand()
                if rand < .5:
                    tmp.add((u,v,u))
                    z=u
                else:
                    tmp.add((u,v,v)) # we cannot use frozenset with 2 equal nodes
                    z=v
                #else: tmp.add(frozenset(u,v,vertices[0]))

                triangles_dict[frozenset((u,v))] = [vertices[0], z]

        elif len(vertices) == 2:
            t1, t2 = frozenset((u,v,vertices[0])), frozenset((u,v,vertices[1]))
            tmp.add(t1)
            tmp.add(t2)
            triangles.add(t1)
            triangles.add(t2)

            triangles_dict[frozenset((u,v))] = vertices



        elif len(vertices) > 2:
            #print("More than 2 neighbors for edge ({},{})".format(u,v))
            #print("Neighbors: ", vertices)
            # Get the edges information
            filtered_edges = edges_df[(edges_df['source'].isin([u, v, *vertices])) | (edges_df['target'].isin([u, v, *vertices]))]
            # convert distance to float
            filtered_edges['distance'] = filtered_edges['distance'].astype(float)
            distances = filtered_edges.groupby(['source', 'target'])['distance'].min().reset_index()
            #print("ok1")
            neighbor_distances = distances[(distances['source'].isin(vertices)) & (distances['target'].isin([u,v])) |
                                            (distances['target'].isin(vertices)) & (distances['source'].isin([u,v]))]
            #print("Neighbor distances: ", neighbor_distances)
            #print("ok2")
            sorted_distances = neighbor_distances.sort_values('distance')
            #print("Sorted distances: ", sorted_distances)
            #print("ok3")
            
            closest_nodes = []

            # Iterate through the sorted distances
            for _, row in sorted_distances.iterrows():
                source, target = row['source'], row['target']

                # append the closest nodes to the list if in vertices and common neighbors of u and v
                if (source in vertices) and (target == u or target == v) and (source not in closest_nodes):
                    closest_nodes.append(source)
                elif (target in vertices) and (source == u or source == v) and (target not in closest_nodes):
                    closest_nodes.append(target)


               # Break the loop if we have found the required number of closest nodes
                if len(closest_nodes) >= 2:
                    break
            
            #print("ok4")
            #print("Closest nodes: ", closest_nodes)
            t1, t2 = frozenset((u,v,closest_nodes[0])), frozenset((u,v,closest_nodes[1]))
            triangles_dict[frozenset((u,v))] = closest_nodes
            #print("Triangles: ", t1, t2)

            triangles.add(t1)
            triangles.add(t2)


        """elif len(vertices) == 0:
            
            print("-------------ERROR--------------------")
            # This shpuldn't happen - SKIPP this part for now, i just kept this code for later 
            t1, t2, t3 = frozenset((u,v,vertices[0])), frozenset((u,v,vertices[1])), frozenset((u,v,vertices[2]))
            tmp.add(t1)
            tmp.add(t2)
            tmp.add(t3)

            for n in vertices:
                if G.degree(n) == 3:
                    print("------ Removing outer triangles -------")
                    try: 
                        tmp.remove(frozenset(G.neighbors(n)))
                        vertices.remove(n)
                    except:
                        print("------No Outer Triangle for: ", list(G.neighbors(n)), " ----------")

            triangles.add(list(tmp)[0])
            triangles.add(list(tmp)[1])

            triangles_dict[frozenset((u,v))] = vertices    """        

    
    return triangles, triangles_dict, edges_df


def extract_features(xc, yc, patch_size, image, model):
    if xc-patch_size//2 < 0:
        xc = patch_size//2
    if yc-patch_size//2 < 0:
        yc = patch_size//2
    if xc+patch_size//2 > image.shape[0]:
        xc = image.shape[0] - patch_size//2
    if yc+patch_size//2 > image.shape[1]:
        yc = image.shape[1] - patch_size//2

    # Extract the patch from the image
    patch = image[xc-patch_size//2:xc+patch_size//2, yc-patch_size//2:yc+patch_size//2, :]
    # Convert the patch to a tensor
    patch_tensor = torch.from_numpy(patch).permute(2, 0, 1).float()
    #patch_tensor = torch.FloatTensor(patch.transpose(2, 0, 1))[None, :]
    # Normalize the patch
    patch_tensor = patch_tensor / 255.0
    # Add a dimension to the tensor to represent the batch
    patch_tensor = patch_tensor.unsqueeze(0)
    # Extract the features of the patch using the pretrained model
    features = model(patch_tensor)
    # Remove the batch dimension
    features = features.squeeze(0)
    # Convert the features to a numpy array
    features = features.detach().numpy()
    return features

def extract_edge_features(x1, y1, x2, y2, image, model):
    """
    if x2 < x1:
        x1, x2 = x2, x1
    if y2 < y1:
        y1, y2 = y2, y1
    # if the patch is smaller than 32x32 then extend it where it is possible
    if x2-x1 < 32:
        x1 = x1 - (32 - (x2-x1))//2
        x2 = x2 + (32 - (x2-x1))//2
    if y2-y1 < 32:
        y1 = y1 - (32 - (y2-y1))//2
        y2 = y2 + (32 - (y2-y1))//2

    # if the patch is outside the image, then move it inside
    if x1 < 0:
        x1 = 0
        x2 = 32
    if y1 < 0:
        y1 = 0
        y2 = 32
    if x2 > image.shape[0]:
        x2 = image.shape[0]
        x1 = x2 - 32
    if y2 > image.shape[1]:
        y2 = image.shape[1]
        y1 = y2 - 32
    """
    xc, yc = (x1+x2)//2, (y1+y2)//2
    patch_size = 128
    # check if the patch is outside the image
    if xc-patch_size//2 < 0:
        xc = patch_size//2
    if yc-patch_size//2 < 0:
        yc = patch_size//2
    if xc+patch_size//2 > image.shape[0]:
        xc = image.shape[0] - patch_size//2
    if yc+patch_size//2 > image.shape[1]:
        yc = image.shape[1] - patch_size//2

    x1, y1 = xc-patch_size//2, yc-patch_size//2
    x2, y2 = xc+patch_size//2, yc+patch_size//2

        
    # Extract the patch from the image
    patch = image[x1:x2, y1:y2, :]
    # Convert the patch to a tensor
    patch_tensor = torch.from_numpy(patch).permute(2, 0, 1).float()
    #patch_tensor = torch.FloatTensor(patch.transpose(2, 0, 1))[None, :]
    # Normalize the patch
    patch_tensor = patch_tensor / 255.0
    # Add a dimension to the tensor to represent the batch
    patch_tensor = patch_tensor.unsqueeze(0)
    # Extract the features of the patch using the pretrained model
    features = model(patch_tensor)
    # Remove the batch dimension
    features = features.squeeze(0)
    # Convert the features to a numpy array
    features = features.detach().numpy()

    
    return features



ADD_NODE_FEATURES = False
ADD_EDGE_FEATURES = False #Result in out of memory error
ADD_MOTIF_FEATURES = False
TRIANGLES_ext = True

class KIGraphDatasetSUBGCN(Dataset):

    def __init__(self, path, mode='train',
                 num_layers=2,
                 data_split=[0.8, 0.2], add_self_edges=False):
        """
        Parameters
        ----------
        path : list
            List with filename, coordinates and path to annotation. For example, ['P7_HE_Default_Extended_1_1', (0, 2000, 0, 2000), 'datasets/annotations/P7_annotated/P7_HE_Default_Extended_1_1.txt']
        mode : str
            One of train, val or test. Default: train.
        num_layers : int
            Number of layers in the computation graph. Default: 2.
        data_split: list
            Fraction of edges to use for graph construction / train / val / test. Default: [0.85, 0.08, 0.02, 0.03].
        """
        super().__init__()

        self.path = path
        self.mode = mode
        self.num_layers = num_layers
        self.data_split = data_split

        print('--------------------------------')
        print('Reading edge dataset from {}'.format(self.path[0]))

        ########## MINE ###########
        # Cells, distance_close_to_edges
        edge_path = path[1]
        node_path = path[2]
    

        # with glob
        edges = pd.read_csv(edge_path)
        nodes = pd.read_csv(node_path)


        if TRIANGLES_ext:
            self.triangles, self.triangles_dict, edges = find_triangles_by_edge(nodes, edges)
        else:
            self.triangles = set()
            self.triangles_dict = dict()

        ## TODO: ADD a condintion to check if the node features are required or not
        if ADD_NODE_FEATURES or ADD_EDGE_FEATURES or ADD_MOTIF_FEATURES:
            #image_path = '..\\..\\intelligraph\\slides\\' + path[0].split('\\')[1] + '.tif'
            image_path = 'datasets\\images\\' + path[0].split('\\')[1] + '.tif'
            image = Image.open(image_path)
            image = image.convert('RGB')
            image = np.array(image)
            self.image = image[:, :, 0:3]
            
            self.model = models.resnet18(weights='ResNet18_Weights.DEFAULT')
            self.model.fc = nn.Identity()
            self.model.eval()

        if add_self_edges:
            for i in range(len(nodes)):
                new_row = {'source': i, 'target': i, 'type': 0, 'distance': 0}
                # append row to the dataframe
                edges = edges.append(new_row, ignore_index=True)

        edges_crossing = edges.copy()
        edges_crossing = edges_crossing[edges_crossing["type"] == 1]

        edges['type'] = edges['type'].replace(1, 0)

        col_row_len = len(nodes['id'])
        
         

       

        distances_close_to_edges = pd.DataFrame(0, index=np.arange(col_row_len), columns=np.arange(col_row_len))
        delta_entropy_edges = pd.DataFrame(0, index=np.arange(col_row_len), columns=np.arange(col_row_len))
        neighborhood_similarity_edges = pd.DataFrame(0, index=np.arange(col_row_len), columns=np.arange(col_row_len))

        for i, row in edges.iterrows():
            source = row['source']
            target = row['target']
            
            distance = float(row['distance'])
           
            delta_entropy = float(row['Delta_Entropy'])
            sorenson_neigh_similarity = float(row['Sorenson_Similarity'])
        
            distances_close_to_edges[source][target] = distance
            distances_close_to_edges[target][source] = distance

            delta_entropy_edges[source][target] = delta_entropy
            delta_entropy_edges[target][source] = delta_entropy

            neighborhood_similarity_edges[source][target] = sorenson_neigh_similarity
            neighborhood_similarity_edges[target][source] = sorenson_neigh_similarity


        distances_close_to_edges = np.array(distances_close_to_edges)
        delta_entropy_edges = np.array(delta_entropy_edges)
        neighborhood_similarity_edges = np.array(neighborhood_similarity_edges)

        # coords
        coords = nodes[["x", "y"]].to_numpy()

        # process neighborhood densities
        density_types = ["Cell_density"]
        #entropy_types = ["Node_Entropy"]

        densities = nodes[density_types].to_numpy()
        edge_density = np.zeros((col_row_len, col_row_len))
        edge_densities = np.empty((0, col_row_len, col_row_len))

        
        for i in range(len(density_types)):
            for _, row in edges.iterrows():
                source = int(row['source'])
                target = int(row['target'])

                edge_density[source][target] = float(densities[:, i][target]) - float(densities[:, i][source])
                edge_density[target][source] = float(densities[:, i][source]) - float(densities[:, i][target])

            edge_densities = np.append(edge_densities, edge_density.reshape(-1, col_row_len, col_row_len), axis=0)

        #print('*************')
        #print('Edge_density Shape : ' + str(edge_densities.shape))
        
        distances_close_to_edges = distances_close_to_edges.reshape(-1, col_row_len, col_row_len)
        delta_entropy_edges = delta_entropy_edges.reshape(-1, col_row_len, col_row_len)
        neighborhood_similarity_edges = neighborhood_similarity_edges.reshape(-1, col_row_len, col_row_len)

        #print('Edge_entropy Shape : ' + str(delta_entropy_edges.shape))
        #print('Edge_distance Shape : ' + str(distances_close_to_edges.shape))
        #print('Neighborhood Similarity Shape : ' + str(neighborhood_similarity_edges.shape))

        if ADD_EDGE_FEATURES:
            edges['morph_features'] = None

            # create sparse tensor of size col_row_len x col_row_len x 512
            morph_features = torch.zeros((col_row_len, col_row_len, 512))

            for _, row in edges.iterrows():
                x1, y1 = nodes.loc[nodes['id'] == row['source']]['x'], nodes.loc[nodes['id'] == row['source']]['y']
                x2, y2 = nodes.loc[nodes['id'] == row['target']]['x'], nodes.loc[nodes['id'] == row['target']]['y']
                features = extract_edge_features(int(x1), int(y1), int(x2), int(y2), self.image, self.model)
                edges.at[_, 'morph_features'] = features

                features = torch.from_numpy(features)
                source = row['source']
                target = row['target']
                
                morph_features[source][target] = features
                morph_features[target][source] = features

            #edge_morph_features = edges['morph_features'].to_numpy()
            #edge_morph_features = np.array(edge_morph_features)


            
            edge_features = np.concatenate((edge_densities, delta_entropy_edges, neighborhood_similarity_edges, distances_close_to_edges, morph_features.permute(2, 0, 1)), axis=0)
        
        else:
            #edge_features = distances_close_to_edges #e1
            #edge_features = np.concatenate((edge_densities, distances_close_to_edges), axis=0) #E12
            #edge_features = np.concatenate((edge_densities, delta_entropy_edges, distances_close_to_edges), axis=0) #E123
            edge_features = np.concatenate((edge_densities, delta_entropy_edges, neighborhood_similarity_edges, distances_close_to_edges), axis=0)
        
        #edge_features = delta_entropy_edges
        #print(edge_features)
        # self.edge_features = utils.normalize_edge_feature_doubly_stochastic(edge_features) ### not to be used
        self.edge_features = utils.normalize_edge_features_rows(edge_features) ### Use it to normalise the edge features
        #self.edge_features = edge_features  ### Use only if not using the normalization feature above


        ## To DO 

        # Change utils.normalize_edge_features_rows to log function to the base e

        #####

        self.channel = edge_features.shape[0]

        self.dist = utils.normalize_edge_features_rows(distances_close_to_edges.reshape(-1, col_row_len, col_row_len))


        # all_labels_cell_types
        nodes["gt"].replace({'inflammatory': 0, 'lymphocyte': 1, 'fibroblast and endothelial': 2, 'epithelial': 3}, inplace=True) # hover-net
                     

        # nuclei features
        #nuclei_feat = nodes[["area", "perim"]].to_numpy()

        all_labels_cell_types = nodes["gt"].to_numpy()

        nodes_with_types_zero_one = nodes.copy()
        nodes_with_types_prob = nodes.copy()
        for i, row in nodes.iterrows():
            ##The iloc positions are based on the csv cell positions
            # if row['gt'] == 'null':
            nodes_with_types_zero_one.iloc[i, 3] = 1 if row['gt'] == 2 else 0
            nodes_with_types_zero_one.iloc[i, 4] = 1 if row['gt'] == 0 else 0
            nodes_with_types_zero_one.iloc[i, 1] = 1 if row['gt'] == 1 else 0
            nodes_with_types_zero_one.iloc[i, 2] = 1 if row['gt'] == 3 else 0
    
        #print('<<<<<<<<<<<<<', nodes_with_types_zero_one)

        
        # cell_types_scores
        cell_types_scores = nodes_with_types_zero_one[['inf', 'lym', 'fib', 'epi']] #One-hot encoding of GT data

        cell_types_scores = cell_types_scores.to_numpy()
        #print(cell_types_scores.shape)

        # adjacency_matrix_close_to_edges
        adjacency_matrix_close_to_edges = np.copy(distances_close_to_edges)
        adjacency_matrix_close_to_edges[adjacency_matrix_close_to_edges != 0] = 1
        self.adj = adjacency_matrix_close_to_edges

        # edge_list_close_to_edge
        edge_list_close_to_edge = edges[["source", "target"]]
        edge_list_close_to_edge = edge_list_close_to_edge.to_numpy()

        # edge_list_crossing_edges
        edge_list_crossing_edges = edges_crossing.to_numpy()

        self.am_close_to_edges_including_distances = distances_close_to_edges
        self.classes = all_labels_cell_types
        self.class_scores = cell_types_scores
        self.coords = coords

        print('Finished reading data.')

        print('Setting up graph.')
        vertex_id = {j: i for (i, j) in enumerate(range(len(coords)))}

        idxs = [floor(v * edge_list_crossing_edges.shape[0]) for v in np.cumsum(data_split)]

        edges_t, pos_examples_crossing_edges = edge_list_close_to_edge, edge_list_crossing_edges

        edges_t[:, :2] = np.array([vertex_id[u] for u in edges_t[:, :2].flatten()]).reshape(edges_t[:, :2].shape)
        edges_t_no_duplicates = np.unique(edges_t[:, :2], axis=0)  # Filter duplicate edges

        self.nodes_count = len(vertex_id)  # Count vertices
        self.edges_count = edges_t_no_duplicates.shape[0]  # Count edges

        adjacency_matrix_close_to_edges = sp.coo_matrix(
            (np.ones(self.edges_count), (edges_t_no_duplicates[:, 0], edges_t_no_duplicates[:, 1])),
            shape=(self.nodes_count, self.nodes_count),
            dtype=np.float32)

        self.adjacency_matrix_close_to_edges_as_coo_to_lil = adjacency_matrix_close_to_edges.tolil()

        self.node_neighbors = self.adjacency_matrix_close_to_edges_as_coo_to_lil.rows  # Neighbors

        #### Code to add node features ### 

        cell_density = nodes['Cell_density'].to_numpy() 
        cell_density = np.array(cell_density)
        cell_density = cell_density.astype(np.float64)

        cell_entropy = nodes['Node_Entropy'].to_numpy()
        cell_entropy = np.array(cell_entropy)
        cell_entropy = cell_entropy.astype(np.float64)

        #mean_neigh_entropy = nodes['Mean_Neighborhood_Entropy'].to_numpy() 
        #mean_neigh_entropy = np.array(mean_neigh_entropy)
        #mean_neigh_entropy = mean_neigh_entropy.astype(np.float64)

        graph_node_features = np.concatenate((cell_types_scores, cell_density[:,None]), axis=1 )# , cell_entropy[:,None]), axis=1 )  ### Concatenate all features 
        



        
        if ADD_NODE_FEATURES:
            nodes['morph_features'] = None
            # For each node, extract the features from the image
            for _, row in nodes.iterrows():
                xc, yc = row['x'], row['y']
                features_resnet = extract_features(int(xc), int(yc), 64, self.image, self.model)
                nodes.at[_, 'morph_features'] = features_resnet

            cell_morph_features = nodes['morph_features'].to_numpy()
            cell_morph_features = np.array(cell_morph_features)
    

            graph_node_features = np.concatenate((cell_types_scores,  np.stack(cell_morph_features, axis=0).astype(np.float64)), axis=1 )  ### Concatenate all features 
            self.features = torch.from_numpy(graph_node_features).float()  # Cell features 
        ###### Code to add node features ends here ##### 
        else:
            ### Use this self feature if only one-hot embedding is required as node feature set
            self.features = torch.from_numpy(cell_types_scores).float()  # Cell features with just one-hot encoding 
            #self.features = torch.from_numpy(graph_node_features).float()  # Cell features 

        self.triangle_morph_features = dict()
        if ADD_MOTIF_FEATURES:
            edges['morph_features'] = None

            # For each edge, extract the features from the image
            for _, row in edges.iterrows():
                x1, y1 = nodes.loc[nodes['id'] == row['source']]['x'], nodes.loc[nodes['id'] == row['source']]['y']
                x2, y2 = nodes.loc[nodes['id'] == row['target']]['x'], nodes.loc[nodes['id'] == row['target']]['y']
                features = extract_edge_features(int(x1), int(y1), int(x2), int(y2), self.image, self.model)
                edges.at[_, 'morph_features'] = features
                u,v = int(row['source']), int(row['target'])
                
                self.triangle_morph_features[frozenset((u,v))] =torch.from_numpy( features).float()
                 

            #triangle_morph_features = edges['morph_features'].to_numpy()
            #self.triangle_morph_features = np.stack( np.array(triangle_morph_features), axis=0).astype(np.float64)
            #self.triangle_morph_features = torch.from_numpy(self.triangle_morph_features).float()  # Cell features

        print('self.features.shape:', self.features.shape)
        # [2] end

        print('Finished setting up graph.')

        print('Setting up examples.')

        if len(pos_examples_crossing_edges) > 0:
            pos_examples_crossing_edges = pos_examples_crossing_edges[:, :2]
            pos_examples_crossing_edges = np.unique(pos_examples_crossing_edges, axis=0)

        # Generate negative examples not in cell edges crossing path
        neg_examples_close_to_edges = []
        cur = 0
        n_count, _choice = self.nodes_count, np.random.choice
        neg_seen = set(tuple(e[:2]) for e in edge_list_crossing_edges)  # Dont sample positive edges
        adj_tuple = set(tuple(e[:2]) for e in edge_list_close_to_edge)  # List all edges

        if self.mode != 'train':  # Add all edges except positive edges if validation/test
            print("self.mode != 'train'")
            for example in edge_list_close_to_edge:
                if (example[0], example[1]) in neg_seen:
                    continue
                neg_examples_close_to_edges.append(example)
            neg_examples_close_to_edges = np.array(neg_examples_close_to_edges, dtype=np.int64)
        else:  # Undersample negative samples from adjacency edges not in positive

            num_neg_examples = pos_examples_crossing_edges.shape[0] # for BCE loss function 

            #If using Focal Loss function use the next line else comment it and uncomment the line above for BCE loss 
            #num_neg_examples = int(pos_examples_crossing_edges.shape[0]) # Increasing the size of neg_samples for focal loss function as it can handle class imbalance 
            
            while cur < num_neg_examples:
                u, v = _choice(n_count, 2, replace=False)
                if (u, v) in neg_seen or (u, v) not in adj_tuple:
                    continue
                cur += 1
                neg_examples_close_to_edges.append([u, v])
            neg_examples_close_to_edges = np.array(neg_examples_close_to_edges, dtype=np.int64)

        x = np.vstack((pos_examples_crossing_edges, neg_examples_close_to_edges))
        y = np.concatenate((np.ones(pos_examples_crossing_edges.shape[0]),
                            np.zeros(neg_examples_close_to_edges.shape[0])))
        perm = np.random.permutation(x.shape[0])
        x, y = x[perm, :], y[perm]  # ERROR HERE -> IndexError: too many indices for array: array is 1-dimensional,
        # but 2 were indexed
        x, y = torch.from_numpy(x).long(), torch.from_numpy(y).long()
        self.x, self.y = x, y

        print('Finished setting up examples.')

        print('Dataset properties:')
        print('Mode: {}'.format(self.mode))
        print('Number of vertices: {}'.format(self.nodes_count))
        print('Number of edges: {}'.format(self.edges_count))
        print('Number of triangles: {}'.format(len(self.triangles)))
        print('Number of positive/negative datapoints: {}/{}'.format(pos_examples_crossing_edges.shape[0],
                                                                     neg_examples_close_to_edges.shape[0]))
        print('Number of examples/datapoints: {}'.format(self.x.shape[0]))

        print('--------------------------------')

    def __len__(self):
        return len(self.x)

    def __getitem__(self, idx):
        return self.x[idx], self.y[idx]

    def get_coords_and_class(self):
        return self.coords, self.classes

    def _form_computation_graph(self, idx):
        """
        Parameters
        ----------
        idx : int or list
            Indices of the node for which the forward pass needs to be computed.
        Returns
        -------
        node_layers : list of numpy array
            node_layers[i] is an array of the nodes in the ith layer of the
            computation graph.
        mappings : list of dictionary
            mappings[i] is a dictionary mapping node v (labelled 0 to |V|-1)
            in node_layers[i] to its position in node_layers[i]. For example,
            if node_layers[i] = [2,5], then mappings[i][2] = 0 and
            mappings[i][5] = 1.
        """
        _list, _set = list, set
        if type(idx) is int:
            node_layers = [np.array([idx], dtype=np.int64)]
        elif type(idx) is list:
            node_layers = [np.array(idx, dtype=np.int64)]

        for _ in range(self.num_layers):
            prev = node_layers[-1]
            arr = [node for node in prev]
            arr.extend([e for node in arr for e in self.node_neighbors[node]])  # add neighbors to graph
            arr = np.array(_list(_set(arr)), dtype=np.int64)
            node_layers.append(arr)
        node_layers.reverse()

        mappings = [{j: i for (i, j) in enumerate(arr)} for arr in node_layers]

        return node_layers, mappings

    def collate_wrapper(self, batch):
        """
        Parameters
        ----------
        batch : list
            A list of examples from this dataset. An example is (edge, label).
        Returns
        -------
        adj : torch.Tensor
            adjacency matrix of entire graph
        features : torch.FloatTensor
            A (n' x input_dim) tensor of input node features.
        edge_features : torch.FloatTensor
            A 3d tensor of edge features.
        edges : numpy array
            The edges in the batch.
        labels : torch.LongTensor
            Labels (1 or 0) for the edges in the batch.
        dist : torch.Tensor
            A distance matrix
        """
        adj = torch.from_numpy(self.adj).float()

        features = self.features
        edge_features = torch.from_numpy(self.edge_features).float()
        edges = np.array([sample[0].numpy() for sample in batch])
        labels = torch.FloatTensor([sample[1] for sample in batch])
        dist = torch.from_numpy(self.dist)

        return adj, features, edge_features, edges, labels, dist, self.triangles_dict, self.triangle_morph_features

    def get_dims(self):
        print("self.features.shape: {}".format(self.features.shape))
        print("input_dims (input dimension) -> self.features.shape[1] = {}".format(self.features.shape[1]))
        return self.features.shape[1], 1

    def get_channel(self):
        return self.channel

    def parse_points(self, fname):
        with open(fname, 'r') as f:
            lines = f.readlines()
        lines = [line[:-1].split(',') for line in lines]  # Remove \n from line
        return lines
##############################################################################################################



def line(p1, p2):
    A = (p1[1] - p2[1])
    B = (p2[0] - p1[0])
    C = (p1[0]*p2[1] - p2[0]*p1[1])
    return A, B, -C

def intersection(L1, L2):
    D  = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    if D != 0:
        x = Dx / D
        return x
    else:
        return False

