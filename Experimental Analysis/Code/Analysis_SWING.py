from Swing import Swing
import numpy as np
import pandas as pd
from scipy import stats
import re
import os
import sys
import time
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
algorithm = "SWING"

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
    
data_experiment = pd.read_csv(f"{path_Data}/{experiment_ID}.tsv", "\t")
data_experiment = data_experiment.iloc[0:number_genes,:]  
            
data_time_series = {}
data_time_series_rand = {}
data_time_series_crand = {}

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
    
    data_time_series_rand[f"data_time_series_rand_{instance}"] = data_time_series[f"data_time_series_{instance}"].sample(frac=1, axis=0)
    
    data_time_series[f"data_time_series_{instance}"].insert(0, "Time", [int(i) for i in string_time])
    data_time_series_rand[f"data_time_series_rand_{instance}"].insert(0, "Time", [int(i) for i in string_time])
    
    data2 = []
    for i in range(len(gene_ID)):
        data2.append(data_experiment[columns_list].iloc[i,:].sample(frac = 1).values.tolist())
    
    
    data_time_series_crand[f"data_time_series_crand_{instance}"] = pd.DataFrame(data2, columns = columns_list, index= gene_ID).T
    data_time_series_crand[f"data_time_series_crand_{instance}"].insert(0, "Time", [int(i) for i in string_time])
    print(data_time_series[f"data_time_series_{instance}"])
    print(data_time_series_rand[f"data_time_series_rand_{instance}"])
    print(data_time_series_crand[f"data_time_series_crand_{instance}"])

data_SWING = pd.DataFrame(columns = data_time_series["data_time_series_1"].columns.values.tolist())
data_SWING_rand = pd.DataFrame(columns = data_time_series_rand["data_time_series_rand_1"].columns.values.tolist())
data_SWING_crand = pd.DataFrame(columns = data_time_series_crand["data_time_series_crand_1"].columns.values.tolist())

        
for instance in instances:

    data_SWING = data_SWING.append(data_time_series[f"data_time_series_{instance}"])
    data_SWING_rand = data_SWING_rand.append(data_time_series_rand[f"data_time_series_rand_{instance}"])
    data_SWING_crand = data_SWING_rand.append(data_time_series_crand[f"data_time_series_crand_{instance}"])

data_SWING.to_csv(f"{folder_results}/data_{experiment_ID}_range{start_instance}-{end_instance}.tsv",sep ="\t", index = False)
data_SWING_rand.to_csv(f"{folder_results}/data_rand_{experiment_ID}_range{start_instance}-{end_instance}.tsv",sep ="\t", index = False)
data_SWING_crand.to_csv(f"{folder_results}/data_crand_{experiment_ID}_range{start_instance}-{end_instance}.tsv",sep ="\t", index = False)

file_path_data = f"{folder_results}/data_{experiment_ID}_range{start_instance}-{end_instance}.tsv"
file_path_data_rand = f"{folder_results}/data_rand_{experiment_ID}_range{start_instance}-{end_instance}.tsv"
file_path_data_crand = f"{folder_results}/data_crand_{experiment_ID}_range{start_instance}-{end_instance}.tsv"

width = round(time_points/2)
method = 'RandomForest'
n_trees = 1000  
k_min = 1
k_max = 3
gene_start_column = 1
time_label = "Time"
separator = "\t"
gene_end = None
header = ["regulator-target","mean_importance"]

#Creating the SWING object
sg_crand = Swing(file_path_data_crand, gene_start_column, gene_end, time_label, separator, min_lag=k_min, max_lag=k_max, window_width=width, window_type=method)
sg = Swing(file_path_data, gene_start_column, gene_end, time_label, separator, min_lag=k_min, max_lag=k_max, window_width=width, window_type=method)
sg_rand = Swing(file_path_data_rand, gene_start_column, gene_end, time_label, separator, min_lag=k_min, max_lag=k_max, window_width=width, window_type=method)


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

ranked_edges = sg.make_sort_df(sg.edge_dict)
ranked_edges_rand = sg_rand.make_sort_df(sg_rand.edge_dict)
ranked_edges_crand = sg_crand.make_sort_df(sg_crand.edge_dict)

ranked_edges.sort_values(by=["regulator-target"]).to_csv(f"{folder_mean_importance}/ranking_{algorithm}_{experiment_ID}_range{start_instance}-{end_instance}.txt", sep = "\t", index = False)
ranked_edges_rand.sort_values(by=["regulator-target"]).to_csv(f"{folder_mean_importance_rand}/ranking_{algorithm}_rand_{experiment_ID}_range{start_instance}-{end_instance}.txt", sep = "\t", index = False)
ranked_edges_crand.sort_values(by=["regulator-target"]).to_csv(f"{folder_mean_importance_crand}/ranking_{algorithm}_crand_{experiment_ID}_range{start_instance}-{end_instance}.txt", sep = "\t", index = False)

os.remove(file_path_data)
os.remove(file_path_data_rand)
os.remove(file_path_data_crand)

end_time = datetime.now()
time = end_time - start_time
print(f'Duration : {time}')
final_message = print("DONE")
        
textfile_time = open(f"{folder_results}/execution_time_{algorithm}_{experiment_ID}_range{start_instance}-{end_instance}.txt", "w")
textfile_time.write(str(time))
textfile_time.close()