from Swing import Swing
import numpy as np
import pandas as pd
import time
from datetime import datetime
import os
import random
import sys

start_instance = int(sys.argv[1])
end_instance = int(sys.argv[2])
number_instances = end_instance - start_instance + 1
genes = int(sys.argv[3])
degree = int{sys.argv[4]}

width = 10
method = 'RandomForest'
n_trees = 1000   
k_min = 1
k_max = 3
gene_start_column = 1
time_label = "Time"
separator = "\t"
gene_end = None
header = ["regulator-target","mean_importance"]

path_results = "Results_{genes}genes_100insilico_networks_{degree}_degree"
try:
    os.mkdir(path_results)
except OSError as error:
    print(f"Folder {path_results} already created.")

path_results_rand = "Results_{genes}genes_100insilico_networks_rand_{degree}_degree"
try:
    os.mkdir(path_results_rand)
except OSError as error:
    print(f"Folder {path_results_rand} already created.")
    
path_results_crand = "Results_{genes}genes_100insilico_networks_crand_{degree}_degree"
try:
    os.mkdir(path_results_crand)
except OSError as error:
    print(f"Folder {path_results_crand} already created.")

for instance in range(number_instances):

    value_iter = start_instance + instance

    data_SWING = pd.read_csv(f"Networks_in_out_degree_{degree}_num_nodes_{genes}/in_and_out_degree_{degree}_num_nodes_{genes}_case_{value_iter}/goldstandard_signed_{value_iter}_dream4_timeseries.tsv", sep = '\t')
    data_SWING_rand = pd.DataFrame(columns = data_SWING.iloc[:,1:].columns.values.tolist())
    data_SWING_crand = pd.DataFrame(columns = data_SWING.iloc[:,1:].columns.values.tolist())
    
    number_time_points = 21
    iterator = 0
    
    for i in range(10):
    
        data2 = []
        for i in range(len(data_SWING.columns.values.tolist()[1:])):
            data2.append(data_SWING[f"{data_SWING.columns.values.tolist()[i+1]}"].iloc[iterator:number_time_points+iterator].sample(frac = 1).values.tolist())
        data3 = pd.DataFrame(np.array(data2).T, columns = data_SWING.columns.values.tolist()[1:])
        
        data_SWING_rand = data_SWING_rand.append(data_SWING.iloc[:,1:].iloc[iterator:iterator+number_time_points].sample(frac = 1, axis = 0))
        data_SWING_crand = data_SWING_crand.append(data3)
        iterator = iterator + number_time_points 
        
    data_SWING_rand.insert(0, "Time", data_SWING["Time"].values.tolist())
    data_SWING_rand.to_csv(f"data_rand_SWING_{value_iter}.tsv",sep ="\t", index = False)
    
    data_SWING_crand.insert(0, "Time", data_SWING["Time"].values.tolist())
    data_SWING_crand.to_csv(f"data_crand_SWING_{value_iter}.tsv",sep ="\t", index = False)
    
    file_path_data = f"Networks_in_out_degree_{degree}_num_nodes_{genes}/in_and_out_degree_{degree}_num_nodes_{genes}_case_{value_iter}/goldstandard_signed_{value_iter}_dream4_timeseries.tsv"
    file_path_data_rand = f"data_rand_SWING_{value_iter}.tsv"
    file_path_data_crand = f"data_crand_SWING_{value_iter}.tsv"
    
    #Creating the SWING object
    sg = Swing(file_path_data, gene_start_column, gene_end, time_label, separator, min_lag=k_min, max_lag=k_max, window_width=width, window_type=method)
    sg_rand = Swing(file_path_data_rand, gene_start_column, gene_end, time_label, separator, min_lag=k_min, max_lag=k_max, window_width=width, window_type=method)
    sg_crand = Swing(file_path_data_crand, gene_start_column, gene_end, time_label, separator, min_lag=k_min, max_lag=k_max, window_width=width, window_type=method)
    
    #Partition of the data into windows
    sg.create_windows()
    sg_rand.create_windows()
    sg_crand.create_windows()
    
    sg.window_list[0]
    sg_rand.window_list[0]
    sg_crand.window_list[0]
    
    #Initialize/Optimize parameters for SWING
    sg.optimize_params()
    sg_rand.optimize_params()
    sg_crand.optimize_params()
    
    #Estimating edge importance
    sg.fit_windows(n_trees=n_trees, show_progress=False, n_jobs=-1)
    sg_rand.fit_windows(n_trees=n_trees, show_progress=False, n_jobs=-1)
    sg_crand.fit_windows(n_trees=n_trees, show_progress=False, n_jobs=-1)
    
    sg.window_list[0].edge_importance
    sg_rand.window_list[0].edge_importance
    sg_crand.window_list[0].edge_importance
    
    #Compile all edges into one table
    sg.compile_edges()
    sg_rand.compile_edges()
    sg_crand.compile_edges()
    
    sg.full_edge_list
    sg_rand.full_edge_list
    sg_crand.full_edge_list
    
    # Get agregated model
    sg.make_static_edge_dict(self_edges=False, lag_method='mean_mean')
    sg_rand.make_static_edge_dict(self_edges=False, lag_method='mean_mean')
    sg_crand.make_static_edge_dict(self_edges=False, lag_method='mean_mean')
    
    ranked_edges= sg.make_sort_df(sg.edge_dict)
    ranked_edges_rand= sg_rand.make_sort_df(sg_rand.edge_dict)
    ranked_edges_crand= sg_crand.make_sort_df(sg_crand.edge_dict)
    
    ranked_edges.sort_values(by=["regulator-target"]).to_csv(f"{path_results}/mean_importance_SWING_{value_iter}.tsv", sep = "\t", index= False, header = header)
    ranked_edges_rand.sort_values(by=["regulator-target"]).to_csv(f"{path_results_rand}/mean_importance_SWING_{value_iter}_rand.tsv", sep = "\t", index= False, header = header)
    ranked_edges_crand.sort_values(by=["regulator-target"]).to_csv(f"{path_results_crand}/mean_importance_SWING_{value_iter}_crand.tsv", sep = "\t", index= False, header = header)

    os.remove(file_path_data)
    os.remove(file_path_data_rand)
    os.remove(file_path_data_crand)