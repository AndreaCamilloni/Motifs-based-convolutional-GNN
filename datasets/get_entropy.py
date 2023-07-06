
import csv
import os
from tempfile import tempdir
from tkinter import E
import pandas as pd
import math

 
def get_entropy(file_names, csv_path, outPath):
     
    for file in file_names:
        edge_csv_file = file.strip()+'_edges.csv'
        node_csv_file = file.strip() + '_nodes.csv'
        adj_list = {}
        with open (csv_path + edge_csv_file) as edges:
            print(file)
            edge_reader = csv.reader(edges, delimiter =',')

            node_cell_type_count = {}
            for eachEdge in edge_reader:
                source , target = eachEdge[0], eachEdge[1]

                node_reader = csv.reader(open(csv_path+node_csv_file, 'r', encoding='iso-8859-1'))
                if source != 'source':    
                    for eachRow in node_reader:
                        if source == eachRow[0]: 
                            source_type = eachRow[7] #get cell-type of the source node
                        if target == eachRow[0]:
                            target_type = eachRow[7] #get cell-type of the target node   

                    ### Add neighbours of source              
                    if source in node_cell_type_count.keys():
                        temp_source = node_cell_type_count[source]
                    else: 
                        temp_source = {'lymphocyte': 0, 'epithelial':0, 'fibroblast and endothelial': 0, 'inflammatory':0,'apoptosis / civiatte body':0, 'total': 0}

                    for eachKey in temp_source.keys():
                        if target_type == eachKey:
                            temp_source[target_type] = int(temp_source[target_type]) + 1 # add the count of the node_type
                            temp_source['total'] = int(temp_source['total']) + 1 # add the total count of all neighbours
                    

                    node_cell_type_count[source] = temp_source

                    ### Add neighbours of target
                    if target in node_cell_type_count.keys():
                        temp = node_cell_type_count[target]
                    else: 
                        temp = {'lymphocyte': '0', 'epithelial':'0', 'fibroblast and endothelial': '0', 'inflammatory':'0','apoptosis / civiatte body':0,'total': '0'}

                    for eachKey in temp.keys(): 
                        if source_type == eachKey: 
                            temp[source_type] = int(temp[source_type]) + 1
                            temp['total'] = int(temp['total']) + 1
                

                    node_cell_type_count[target] = temp
                        
            
            node_reader = csv.reader(open(csv_path+node_csv_file, 'r', encoding='iso-8859-1'))
            writer = csv.writer(open(outPath+node_csv_file, 'w', newline=''))   
            node_entropy_dict = {}
            for eachRow in node_reader: 
                eachRow = eachRow
                node_entropy = 'Node_entropy'
                id = eachRow[0]
             
                if id in node_cell_type_count.keys(): 
                    
                    id_type = eachRow[7]
                    entropy_list = node_cell_type_count[id]

                    ###Add self node type for entropy calculation
                    entropy_list[id_type] = int(entropy_list[id_type]) +1
                    entropy_list['total'] = int(entropy_list['total']) +1

                    ###Calculate node type probabilities
                    prob_epi = float(entropy_list['epithelial'])/ float(entropy_list['total'])
                    prob_lym = float(entropy_list['lymphocyte'])/ float(entropy_list['total'])
                    prob_fib = float(entropy_list['fibroblast and endothelial'])/ float(entropy_list['total'])
                    prob_inf = float(entropy_list['inflammatory'])/ float(entropy_list['total'])
                    prob_apop = float(entropy_list['apoptosis / civiatte body'])/float(entropy_list['total'])

                    if prob_epi == 0.0: 
                        entropy_epi = 0
                    else: 
                        entropy_epi = abs(prob_epi * math.log(prob_epi)) 

                    if prob_lym == 0.0: 
                        entropy_lym = 0
                    else: 
                        entropy_lym = abs(prob_lym * math.log(prob_lym)) 

                    if prob_fib == 0.0: 
                        entropy_fib = 0
                    else: 
                        entropy_fib = abs(prob_fib * math.log(prob_fib)) 

                    if prob_inf == 0.0: 
                        entropy_inf = 0
                    else: 
                        entropy_inf = abs(prob_inf * math.log(prob_inf)) 

                    if prob_apop == 0.0:
                        entropy_apop = 0.0
                    else: 
                        entropy_apop = abs(prob_apop * math.log(prob_apop)) 


                    node_entropy = entropy_epi + entropy_fib + entropy_inf + entropy_lym + entropy_apop
                    node_entropy_dict[id] = node_entropy
                
                else:
                   node_entropy = 'Node_Entropy' 
                
                eachRow.append(node_entropy)
                writer.writerow(eachRow)
        edges.close()

        with open (csv_path + edge_csv_file) as edges:
            edge_reader = csv.reader(edges, delimiter =',')
            edge_writer = csv.writer(open(outPath+edge_csv_file, 'w',newline=''))
            delta_entropy = 'Delta_Entropy'
            for eachEdge in edge_reader:
                source , target = eachEdge[0], eachEdge[1]
                if source != 'source' and target != 'target' and node_entropy_dict[source] != 'Node_Entropy' and node_entropy_dict[target] != 'Node_Entropy':
                    delta_entropy = abs(float(node_entropy_dict[source]) - float(node_entropy_dict[target]))
                else: 
                    delta_entropy = 'Delta_Entropy'
                eachEdge.append(delta_entropy)
                edge_writer.writerow(eachEdge)
        edges.close()

 