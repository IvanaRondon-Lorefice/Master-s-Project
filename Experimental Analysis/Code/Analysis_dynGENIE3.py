from dynGENIE3 import *
import numpy as np
import pandas as pd
from scipy import stats
import re
import os
import time
import sys
import collections
from datetime import datetime
from random import shuffle
            
start_time = datetime.now()

start_instance = int(sys.argv[1])    
end_instance = int(sys.argv[2])
experiment_ID = str(sys.argv[3])
number_genes = int(sys.argv[4])
number_instances = int(sys.argv[5])

path_Data = "Data" 
algorithm = "dynGENIE3"

instances = [i+1 for i in range(end_instance - start_instance +1)]


path_results = f"Results_Analysis_Experimental_Data_instances-{number_instances}"
try:
    os.mkdir(path_results)
except OSError as error:
    print(f"Folder {path_results} already created.")
            
folder_results = f"{path_results}/Results_{algorithm}_{experiment_ID}_genes-{number_genes}"
try:
    os.mkdir(folder_results)
except OSError as error:
    print(f"Folder {folder_results} already created.")
        
folder_mean_importance = f"{folder_results}/Mean_Importance_{experiment_ID}"
try:
    os.mkdir(folder_mean_importance)
except OSError as error:
    print(f"Folder {folder_mean_importance} already created")
            
folder_mean_importance_rand = f"{folder_results}/Mean_Importance_{experiment_ID}_rand"
try:
    os.mkdir(folder_mean_importance_rand)
except OSError as error:
    print(f"Folder {folder_mean_importance_rand} already created")

folder_mean_importance_crand = f"{folder_results}/Mean_Importance_{experiment_ID}_crand"
try:
    os.mkdir(folder_mean_importance_crand)
except OSError as error:
    print(f"Folder {folder_mean_importance_crand} already created")
    
folder_columns_data = f"{folder_results}/Columns_data_{experiment_ID}"
try:
    os.mkdir(folder_columns_data)
except OSError as error:
    print(f"Folder {folder_columns_data} already created")   

header_list = ["Cause Gene", "Effect Gene", "mean_importance"]
tree_method = 'RF'    #Random Forest
ntrees = 1000       

data_experiment = pd.read_csv(f"{path_Data}/{experiment_ID}.tsv", "\t")
data_experiment = data_experiment.iloc[0:number_genes,:]  
            
data_time_series = {}
data_time_series_rand = {}

for instance in instances:
                
    columns_data = data_experiment.columns[1:]
    string_time = list(collections.Counter([re.search('_T(.+?)_', columns_data[i]).group(1) for i in range(len(columns_data))]))
    time_points = len(string_time)
    step_time = int(string_time[1]) - int(string_time[0])
    columns_list = [list(filter(re.compile(".*_T{}_*".format(i)).match, columns_data)) for i in string_time]
    columns_list = [columns_list[i][np.random.randint(0, len(columns_list[i]))] for i in range(len(string_time))]
                
    step = instance +  start_instance - 1
                
    columns_list_file = open(f"{folder_columns_data}/colums_data_instance-{step}.txt", "w")
                
    for element in columns_list:
        columns_list_file.write(str(element) + "\n")
    columns_list_file.close()
                
    gene_ID = data_experiment["gene_id"].values.tolist()
   
    data_time_series[f"data_time_series_{instance}"] = pd.DataFrame(data_experiment[columns_list].values,index= gene_ID, columns = columns_list).T
    data_time_series[f"data_time_series_{instance}"].insert(0, "'Time'", [int(i) for i in string_time])
    

data_dynGENIE3 = []
data_dynGENIE3_rand = []
data_dynGENIE3_crand = []
time_list = []
iterator = 0
number_time_points = len(string_time)


for instance in instances:
    
    data_dynGENIE3.append(data_time_series[f"data_time_series_{instance}"].iloc[:,1:].values)
    data_dynGENIE3_rand.append(data_time_series[f"data_time_series_{instance}"].iloc[:,1:].sample(frac=1, axis = 0).values)
    time_list.append(data_time_series[f"data_time_series_{instance}"].iloc[:,0].values) 
    
    data2 = []
        
    for i in range(len(data_time_series[f"data_time_series_{instance}"].iloc[:,1:].columns.values.tolist())):
        data2.append(data_time_series[f"data_time_series_{instance}"].iloc[:,i+1].sample(frac = 1).values.tolist())
    
    data3 = pd.DataFrame(np.array(data2).T, columns = data_time_series[f"data_time_series_{instance}"].iloc[:,1:].columns.values.tolist())

    data_dynGENIE3_crand.append(data3.values)
 
          

(VIM_crand, alphas_crand, prediction_score_crand, stability_score_crand, treeEstimators_crand) = dynGENIE3(data_dynGENIE3_crand, time_list, tree_method=tree_method,compute_quality_scores=True, ntrees=ntrees)
get_link_list(VIM_crand, gene_names = gene_ID, file_name=f'{folder_mean_importance_crand}/ranking_{algorithm}_crand_{experiment_ID}_range{start_instance}-{end_instance}.txt')
           
(VIM, alphas, prediction_score, stability_score, treeEstimators) = dynGENIE3(data_dynGENIE3, time_list, tree_method=tree_method,compute_quality_scores=True, ntrees=ntrees)
get_link_list(VIM, gene_names = gene_ID, file_name=f'{folder_mean_importance}/ranking_{algorithm}_{experiment_ID}_range{start_instance}-{end_instance}.txt')

(VIM_rand, alphas_rand, prediction_score_rand, stability_score_rand, treeEstimators_rand) = dynGENIE3(data_dynGENIE3_rand, time_list, tree_method=tree_method,compute_quality_scores=True, ntrees=ntrees)
get_link_list(VIM_rand, gene_names = gene_ID, file_name=f'{folder_mean_importance_rand}/ranking_{algorithm}_rand_{experiment_ID}_range{start_instance}-{end_instance}.txt')



end_time = datetime.now()
time = end_time - start_time
print(f'Duration : {time}')
final_message = print("DONE")
        
textfile_time = open(f"{folder_results}/execution_time_{algorithm}_{experiment_ID}_range{start_instance}-{end_instance}.txt", "w")
textfile_time.write(str(time))
textfile_time.close()