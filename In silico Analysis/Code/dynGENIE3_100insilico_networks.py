from dynGENIE3 import *
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
degree = int(sys.argv[4])

path_results = f"Results_{genes}genes_100insilico_networks_{degree}_degree"
try:
    os.mkdir(path_results)
except OSError as error:
    print(f"Folder {path_results} already created.")

path_results_rand = f"Results_{genes}genes_100insilico_networks_rand_{degree}_degree"
try:
    os.mkdir(path_results_rand)
except OSError as error:
    print(f"Folder {path_results_rand} already created.")
    
path_results_crand = f"Results_{genes}genes_100insilico_networks_crand_{degree}_degree"
try:
    os.mkdir(path_results_crand)
except OSError as error:
    print(f"Folder {path_results_crand} already created.")
    
header_list = ["Cause Gene", "Effect Gene", "mean_importance"]
tree_method = 'RF'    #Random Forest
ntrees = 1000   

number_time_points = 21
number_repetitions = 10

for instance in range(number_instances):

    value_iter = start_instance + instance
    data = pd.read_csv(f"Networks_in_out_degree_{degree}_num_nodes_{genes}/in_and_out_degree_{degree}_num_nodes_{genes}_case_{value_iter}/goldstandard_signed_{value_iter}_dream4_timeseries.tsv", sep = '\t')
    data_dynGENIE3 = []
    data_dynGENIE3_rand = []
    time_list = []
    iterator = 0
    gene_names = data.columns.values[1:].tolist()
    
    data_dynGENIE3_crand = []
    
    for j in range(number_repetitions):
    
        data_dynGENIE3.append(data.iloc[:,1:].iloc[iterator:number_time_points+iterator].values)
        data_dynGENIE3_rand.append(data.iloc[:,1:].iloc[iterator:number_time_points+iterator].sample(frac = 1, axis = 0).values)
        time_list.append(data.iloc[:,0].iloc[iterator:number_time_points+iterator].values)
        
        data2 = []
        
        for i in range(len(data.columns.values.tolist()[1:])):
            data2.append(data[f"{data.columns.values.tolist()[i+1]}"].iloc[iterator:number_time_points+iterator].sample(frac = 1).values.tolist())
        
        data3 = pd.DataFrame(np.array(data2).T, columns = data.columns.values.tolist()[1:])
        
        data_dynGENIE3_crand.append(data3.values)
        
        iterator = iterator + number_time_points
        
    (VIM, alphas, prediction_score, stability_score, treeEstimators) = dynGENIE3(data_dynGENIE3, time_list, tree_method=tree_method, ntrees=ntrees, compute_quality_scores=True)
    (VIM_rand, alphas_rand, prediction_score_rand, stability_score_rand, treeEstimators_rand) = dynGENIE3(data_dynGENIE3_rand, time_list, tree_method=tree_method, ntrees=ntrees, compute_quality_scores=True)
    (VIM_crand, alphas_crand, prediction_score_crand, stability_score_crand, treeEstimators_crand) = dynGENIE3(data_dynGENIE3_crand, time_list, tree_method=tree_method, ntrees=ntrees, compute_quality_scores=True)
    
    get_link_list(VIM,gene_names=gene_names,file_name=f'{path_results}/mean_importance_dynGENIE3_{value_iter}.tsv')
    get_link_list(VIM_rand,gene_names=gene_names,file_name=f'{path_results_rand}/mean_importance_dynGENIE3_{value_iter}_rand.tsv')
    get_link_list(VIM_crand,gene_names=gene_names,file_name=f'{path_results_crand}/mean_importance_dynGENIE3_{value_iter}_crand.tsv')
    


