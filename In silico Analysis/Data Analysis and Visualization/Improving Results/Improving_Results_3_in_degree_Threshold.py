import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import style
from datetime import datetime
import os
style.use('fast')
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import average_precision_score  
from statistics import mean

def ROC_PR_AUC_Curves(Values_real,Value_infered, number_genes, file_path_data, algorithm, color_ROC="indigo", color_PR= "springgreen", show_plot=False):
    """
        This functions returns the area under the ROC and PR curve 
        and plots the tpr vs fpr and the precision vs recall.
    """
    
    #ROC
    fpr, tpr, thresholds = roc_curve(Values_real, Value_infered)
    AUC_ROC = auc(fpr, tpr)
    
    
    #PR
    precision, recall, thresholds = precision_recall_curve(Values_real, Value_infered)
    AUC_PR = auc(recall, precision)
    
    if show_plot==True:
        
        print("No plot")
    
    return AUC_ROC, AUC_PR

# Variating over different thresholds
header_list = ["Cause Gene", "Effect Gene", "mean_importance"]
Threshold = np.arange(0,0.15,0.001)
Nodes = [10,20,50,100]
In_degree = [3]
         
Algorithm = ["SWING", "dynGENIE3"]
AUROC = {}
AUPRC = {}
ranking = pd.DataFrame()

for nodes in Nodes:
    
    print(nodes)
    
    for degree in In_degree:
        for algorithm in Algorithm:
            for threshold in Threshold:
                
                AUROC[f"{algorithm}_{nodes}Nodes_{degree}Degree_{threshold}"] = []
                AUPRC[f"{algorithm}_{nodes}Nodes_{degree}Degree_{threshold}"] = []
                
                AUROC[f"{algorithm}_{nodes}Nodes_{degree}Degree_{threshold}_rand"] = []
                AUPRC[f"{algorithm}_{nodes}Nodes_{degree}Degree_{threshold}_rand"] = []
                
                for instance in range(100):
            
                    instance = instance + 1
        
                    
                    if nodes == 10:                                                                                                                                   
                        gold_standard = pd.read_csv(f"networks_in_out_degree_3_num_nodes_{nodes}/Networks_in_out_degree_3_num_nodes_{nodes}/in_and_out_degree_2_num_nodes_{nodes}_case_{instance}/goldstandard_signed_{instance}_goldstandard.tsv", sep="\t", names= header_list)
                    else:
                        gold_standard = pd.read_csv(f"networks_in_out_degree_3_num_nodes_{nodes}/Networks_in_out_degree_3_num_nodes_{nodes}/in_and_out_degree_3_num_nodes_{nodes}_case_{instance}/goldstandard_signed_{instance}_goldstandard.tsv", sep="\t", names= header_list)
                    
                    gold_standard = gold_standard.sort_values(by=["Cause Gene", "Effect Gene"])

                    if algorithm == "SWING":    
                                                                
                        mean_importance_file = pd.read_csv(f"networks_in_out_degree_{degree}_num_nodes_{nodes}/Results_{nodes}genes_100insilico_networks/mean_importance_SWING_{instance}.tsv", sep ="\t").sort_values(by=["regulator-target"])
                        mean_importance_file_rand = pd.read_csv(f"networks_in_out_degree_{degree}_num_nodes_{nodes}/Results_{nodes}genes_100insilico_networks_rand/mean_importance_SWING_{instance}_rand.tsv", sep ="\t").sort_values(by=["regulator-target"])
                        
                        rank = mean_importance_file["mean_importance"]
                        rank_rand = mean_importance_file_rand["mean_importance"]
                    
                        new_ranking = pd.DataFrame(columns = ["regulator-target","mean_importance"])
                        
                        for i in range(len(rank)):
                        
                            if rank.iloc[i] > threshold:
                                if (rank.iloc[i]) > rank_rand.iloc[i]:
                                    new_ranking.loc[i] = [mean_importance_file["regulator-target"].iloc[i], 1]
                                else:
                                    new_ranking.loc[i] = [mean_importance_file["regulator-target"].iloc[i], rank.iloc[i]]
                            else:
                                new_ranking.loc[i] = [mean_importance_file["regulator-target"].iloc[i], rank.iloc[i]] 
                            
                        new_ranking = new_ranking.sort_values(by=["regulator-target"])
                        
                    else:
                        mean_importance_file = pd.read_csv(f"networks_in_out_degree_{degree}_num_nodes_{nodes}/Results_{nodes}genes_100insilico_networks/mean_importance_dynGENIE3_{instance}", sep ="\t").sort_values(by=["Cause Gene", "Effect Gene"])
                        mean_importance_file_rand = pd.read_csv(f"networks_in_out_degree_{degree}_num_nodes_{nodes}/Results_{nodes}genes_100insilico_networks_rand/mean_importance_dynGENIE3_{instance}_rand.tsv", sep ="\t", names = header_list).sort_values(by=["Cause Gene", "Effect Gene"])
                        
                        rank = mean_importance_file["mean_importance"]
                        rank_rand = mean_importance_file_rand["mean_importance"]
                    
                        new_ranking = pd.DataFrame(columns = ["Cause Gene","Effect Gene","mean_importance"])
                        
                        for i in range(len(rank)):
                        
                            if rank.iloc[i] > threshold:
                                if (rank.iloc[i]) > rank_rand.iloc[i]:
                                    new_ranking.loc[i] = [mean_importance_file["Cause Gene"].iloc[i],mean_importance_file["Effect Gene"].iloc[i], 1]
                                else:
                                    new_ranking.loc[i] = [mean_importance_file["Cause Gene"].iloc[i],mean_importance_file["Effect Gene"].iloc[i], rank.iloc[i]]
                            else:
                                new_ranking.loc[i] = [mean_importance_file["Cause Gene"].iloc[i],mean_importance_file["Effect Gene"].iloc[i], rank.iloc[i]]
                            
                        new_ranking = new_ranking.sort_values(by=["Cause Gene","Effect Gene"])
                        
                    AUROC_value, AUPRC_value = ROC_PR_AUC_Curves(gold_standard["mean_importance"], new_ranking["mean_importance"] , nodes, "Nothing to plot", f"{algorithm} Run {instance}", show_plot=False)
                    AUROC[f"{algorithm}_{nodes}Nodes_{degree}Degree_{threshold}"].append(AUROC_value)
                    AUPRC[f"{algorithm}_{nodes}Nodes_{degree}Degree_{threshold}"].append(AUPRC_value)

                    AUROC_value_rand, AUPRC_value_rand = ROC_PR_AUC_Curves(gold_standard["mean_importance"], rank_rand, nodes, "Nothing to plot", f"{algorithm} Run {instance}", show_plot=False)
                    AUROC[f"{algorithm}_{nodes}Nodes_{degree}Degree_{threshold}_rand"].append(AUROC_value_rand)
                    AUPRC[f"{algorithm}_{nodes}Nodes_{degree}Degree_{threshold}_rand"].append(AUPRC_value_rand)

            
            try:
                os.mkdir(f"Values_Threshold_{degree}_in_degree")
            except OSError as error:
                print(f"folder Values_Threshold_{degree}_in_degree already created")
    
            values_threshold = [(np.power(np.prod(AUROC[f"{algorithm}_{nodes}Nodes_{degree}Degree_{threshold}"]),1/100) + np.power(np.prod(AUPRC[f"{algorithm}_{nodes}Nodes_{degree}Degree_{threshold}"]),1/100))/2 for threshold in Threshold]
            values_threshold_files = open(f"Values_Threshold_{degree}_in_degree/values_threshold_{degree}_{nodes}_{algorithm}.txt", "w")
            for element in values_threshold:
                values_threshold_files.write(str(element)+"\n")
            values_threshold_files.close()
