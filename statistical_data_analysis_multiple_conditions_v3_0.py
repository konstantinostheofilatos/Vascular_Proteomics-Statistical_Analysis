

#python statistical_data_analysis_timeseries.py data_example.txt 4 tags.txt
#python statistical_data_analysis_timeseries.py miRNA_protein_data.txt 12 miRNA_protein_tags.txt
#python statistical_data_analysis_timeseries.py tissue_crhistian_data.txt 6 christian_labels.txt
import os
import statistics
import numpy as np
import scipy.stats as st
import time
import math
import sys
import csv
from numpy import matrix
import matplotlib.pyplot as plt
import plotly.plotly as py
import statsmodels.stats.multitest as smm
import plotly.tools as tls
import plotly
from scipy.stats import mstats
import copy
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler
from itertools import compress
from sklearn.preprocessing import Imputer
from knnimpute import (
    knn_impute_few_observed,
    knn_impute_with_argpartition,
    knn_impute_optimistic,
    knn_impute_reference,
)
os.environ['R_HOME'] = 'C:/Program Files/R/R-3.4.2'
os.environ['R_USER'] = 'C:/Users/Konstantinos/AppData/Local/Programs/Python/Python36/Lib/site-packages/rpy2'

from rpy2 import *
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
ro.numpy2ri.activate()
from rpy2.robjects import r
import rpy2.robjects.lib.ggplot2 as ggplot2



gplots = importr("gplots")
import rpy2.robjects as robjects 
from rpy2.robjects.packages import importr 
limma = importr('limma') 
statmod=importr('statmod')
ggrepel=importr('ggrepel')
lattice=importr('lattice')

def parse_data(data_filename):
	num_of_lines=0
	proteins=list()
	data=list()
	with open(data_filename) as data_fname:
			for line in csv.reader(data_fname, delimiter="\t"):
				if num_of_lines==0:
					for j in range(len(line)):
						proteins.append(line[j].strip().replace(" ", "").replace('/',''))
				else:
					
					data.append([])
					for j in range(len(line)):
						if line[j]!='' and line[j]!=-1000:
							data[num_of_lines-1].append(float(line[j]))
						elif line[j]!='':
							data[num_of_lines-1].append(-1000)
						else:
							data[num_of_lines-1].append('')
				num_of_lines+=1
	print('Data were successfully parsed!')
	return [proteins,data]

def parse_samples(samples_filename):
	samples=list()
	num_of_lines=0
	with open(samples_filename) as samples_fname:
			for line in csv.reader(samples_fname, delimiter="\t"):
				if num_of_lines>0:
					samples.append(line[1])
				num_of_lines+=1
	return samples
	
def filter_data(data, samples):
	filtered_data=list()
	for sample in samples:
		filtered_data.append(data[sample])
	return filtered_data

def fold_data(data, timepoints):
	data_new=copy.deepcopy(data)
	for j in range(len(data_new[0])):
		for i in range(len(data_new)):
			if i%timepoints==0:
				for k in range(timepoints):
					if data_new[i][j]==0:
						data_new[i+timepoints-i][j]=0
					else:
						data_new[i+timepoints-i][j]/=data_new[i][j]
	return data_new
def data_per_timepoints(data, timepoints):
	splitted_data=list()
	for i in range(timepoints):
		temp_table=[]
		for j in range(len(data)):
			if j%timepoints==i:
				temp_table.append(data[j])
		splitted_data.append(temp_table)		
	return splitted_data
	
def print_boxplots(proteins, data, annotation, folder, folded, tags):
	for j in range(len(proteins)):
	
		splitted_data=list()
		for i in range(len(tags)):
			temp_table=[]
			for k in range(len(data)):
				if annotation[k]==tags[i]:
					temp_table.append(data[k][j])
			splitted_data.append(temp_table)	
		mpl_fig = plt.figure(figsize=((len(tags)+1), 5))
		mpl_fig.subplots_adjust(bottom=0.25,left=0.2)
		ax = mpl_fig.add_subplot(111)
		ax.boxplot(splitted_data, widths=0.02, positions= (np.arange(len(tags))+0.5)/(len(tags)+1) )
		#ax.set_aspect(1.2)
		ax.set_title(proteins[j])
		if folded==1:
			ax.set_ylabel('Folded Quantities')
		else:
			ax.set_ylabel('Logged Relative Quantities')
		ax.set_xlabel('Phenotypes')
		labels = [item.get_text() for item in ax.get_xticklabels()]
		for i in range(len(tags)):
			labels[i]=tags[i]
		ax.set_xticklabels(labels, rotation=45)
		ax.set_xticks((np.arange(len(labels))+0.5)/(len(tags)+1))
		ax.set_xlim(right=0.9, left=0)
		if not os.path.exists(folder):
			os.makedirs(folder)
		mpl_fig.savefig(folder+proteins[j]+'_boxplot.png')
		mpl_fig.clf()
		plt.close('all')
		

def print_significant_boxplots(proteins, data, annotation, folder, folded,pvals,tags,p_value_threshold):
	for j in range(len(proteins)):
		if(pvals[j]<float(p_value_threshold)):
			splitted_data=list()
			for i in range(len(tags)):
				temp_table=[]
				for k in range(len(data)):
					if annotation[k]==tags[i]:
						temp_table.append(data[k][j])
				splitted_data.append(temp_table)	
			mpl_fig = plt.figure(figsize=((len(tags)+1), 5))
			mpl_fig.subplots_adjust(bottom=0.25,left=0.2)
			ax = mpl_fig.add_subplot(111)
			ax.boxplot(splitted_data, widths=0.02, positions= (np.arange(len(tags))+0.5)/(len(tags)+1) )
			#ax.set_aspect(1.2)
			ax.set_title(proteins[j])
			if folded==1:
				ax.set_ylabel('Folded Quantities')
			else:
				ax.set_ylabel('Logged Relative Quantities')
			ax.set_xlabel('Phenotypes')
			labels = [item.get_text() for item in ax.get_xticklabels()]
			for i in range(len(tags)):
				labels[i]=tags[i]
			ax.set_xticklabels(labels, rotation=45)
			ax.set_xticks((np.arange(len(labels))+0.5)/(len(tags)+1))
			ax.set_xlim(right=0.9, left=0)
			if not os.path.exists(folder):
				os.makedirs(folder)
			mpl_fig.savefig(folder+proteins[j]+'_boxplot.png')
			mpl_fig.clf()
			plt.close('all')		

def kruskal_wallis_test(data,labels,markers):
	Hs=list()
	pvals=list()
	#tags=list(set(labels))
	#tags=['Fibroatheroma','Complex','Calcified','Fibrotic']
	tags=['A','C','D','v1H','v2H']
	for j in range(len(data[0])):
		splitted_data=dict()
		for i in range(len(tags)):
			splitted_data[i]=list()
			for k in range(len(data)):
				if labels[k]==tags[i]:
					splitted_data[i].append(data[k][j])
		
		try:
			H, pval =mstats.kruskalwallis(*splitted_data.values())
			
		except:
			H=0
			pval=1
		Hs.append(H)
		pvals.append(pval)
	return [Hs,pvals]
	
def anova_test(data,labels,markers):
	Hs=list()
	pvals=list()
	#tags=list(set(labels))
	#tags=['Fibroatheroma','Complex','Calcified','Fibrotic']
	tags=['A','C','D','v1H','v2H']
	for j in range(len(data[0])):
		splitted_data=dict()
		
		
		for i in range(len(tags)):
			splitted_data[i]=list()
			for k in range(len(data)):
				if labels[k]==tags[i]:
					splitted_data[i].append(data[k][j])
		
		try:
			H, pval =mstats.f_oneway(*splitted_data.values())
		
		except:
			H=0
			pval=1
		Hs.append(H)
		pvals.append(pval)
	return [Hs,pvals]
	


def create_clustering_dataset(proteins, data, filename):
	data_for_clustering=list()
	for j in range(len(proteins)):
		data_for_clustering.append([])
		data_for_clustering[j].append(proteins[j])
		splitted_data=list()
		for i in range(timepoints):
			temp_table=[]
			for k in range(len(data)):
				if k%timepoints==i:
					temp_table.append(data[k][j])
			data_for_clustering[j].append(sum(temp_table)/float(len(temp_table)))
		splitted_data.append(temp_table)
	with open(filename,'wb') as resultFile:
		wr = csv.writer(resultFile, dialect='excel')
		wr.writerows(data_for_clustering)
		
def cluster_proteins(proteins,data, filename,name):
	data_for_clustering=list()
	for j in range(len(proteins)):
		data_for_clustering.append([])
		#data_for_clustering[j].append(proteins[j])
		splitted_data=list()
		for i in range(timepoints):
			temp_table=[]
			for k in range(len(data)):
				if k%timepoints==i:
					temp_table.append(data[k][j])
			data_for_clustering[j].append(sum(temp_table)/float(len(temp_table)))
		splitted_data.append(temp_table)
	db = DBSCAN(eps=1, min_samples=2).fit(data_for_clustering)
	core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
	core_samples_mask[db.core_sample_indices_] = True
	labels = db.labels_
	print(labels)
	print(len(labels))
	# Number of clusters in labels, ignoring noise if present.
	n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

	print('Estimated number of clusters: %d' % n_clusters_)
	
	print("Silhouette Coefficient: %0.3f"
		  % metrics.silhouette_score(data_for_clustering, labels))


	

	# Black removed and is used for noise instead.
	unique_labels = set(labels)
	print(unique_labels)
	colors = [plt.cm.Spectral(each)
			  for each in np.linspace(0, 1, len(unique_labels))]
	resultFile= open(filename,'w') 
	message=''
	for k, col in zip(unique_labels, colors):
		

		class_member_mask = (labels == k)
		data_for_clustering_filtered=list()
		mpl_fig = plt.figure()
		ax = mpl_fig.add_subplot(111)
		protein_list=''
		for iter in range(len(data_for_clustering)):
			#print(labels)
			#print(iter)
			if labels[iter]==k:
				#data_for_clustering_filtered.append(copy.deepcopy(data_for_clustering_filtered))
				protein_list=protein_list+'\t'+proteins[iter]
				ax.plot(data_for_clustering[iter])
		
		
		
		
		ax.set_title('Cluster_'+str(k))
		ax.set_ylabel('Folded Concetrations')
		ax.set_xlabel('Time-points')
		#labels_n=list()
		labels_n = [item.get_text() for item in ax.get_xticklabels()]
		
		labels_n [1] = '0'
		labels_n [3] ='60min'
		labels_n [5] = '8h'
		labels_n [7] ='24h'
		ax.set_xticklabels(labels_n)
		#mpl_fig.Title(proteins[j])
		mpl_fig.savefig(name+'_Cluster_'+str(k)+'_figure.png')
		mpl_fig.clf()
		
		message+='Cluster_'+str(k)+':'+protein_list+'\n'
		#
	resultFile.write(message)
	resultFile.close()	

def wilcoxon_rank_sum_test(proteins,control, condition, paired_flag ):
	pvals=list()
	Zs=list()
	upregulated_proteins=list()
	downregulated_proteins=list()
	num_of_ups=0
	num_of_downs=0
	folds=list()
	stdevs=list()
	for i in range(len(proteins)):
		control_data_per_protein=list()
		condition_data_per_protein=list()
		for j in range(len(control)):
			control_data_per_protein.append(control[j][i])
		for j in range(len(condition)):
			condition_data_per_protein.append(condition[j][i])
		
		#[z,pval]=st.wilcoxon(control_data_per_protein, condition_data_per_protein)
		if statistics.stdev(control_data_per_protein)==0 and statistics.stdev(condition_data_per_protein)==0 and statistics.mean(control_data_per_protein)==statistics.mean(condition_data_per_protein):
			pval=1
		else:
			if paired_flag==0:
				[z,pval]=st.ranksums(control_data_per_protein, condition_data_per_protein)
			else:
				[z,pval]=st.wilcoxon(control_data_per_protein, condition_data_per_protein)
			
		pvals.append(pval)
		
		if paired_flag==1:
			fold_per_sample=list()
				
			for k in range(len(condition_data_per_protein)):
				fold_per_sample.append(condition_data_per_protein[k]-control_data_per_protein[k])
		if statistics.mean(control_data_per_protein)==0:
			folds.append(0)
			stdevs.append(0)
		else:
			folds.append(statistics.mean(condition_data_per_protein)-statistics.mean(control_data_per_protein))
			if paired_flag==1:
				stdevs.append(statistics.stdev(fold_per_sample))
			else:
				stdevs.append(0)
	if paired_flag==0:
		stdevs=[0]*len(proteins)
	return [pvals,folds,stdevs]
	
def t_test(proteins,control, condition, paired_flag ):
	pvals=list()
	Zs=list()
	upregulated_proteins=list()
	downregulated_proteins=list()
	num_of_ups=0
	num_of_downs=0
	folds=list()
	stdevs=list()
	for i in range(len(proteins)):
		control_data_per_protein=list()
		condition_data_per_protein=list()
		for j in range(len(control)):
			control_data_per_protein.append(control[j][i])
		for j in range(len(condition)):
			condition_data_per_protein.append(condition[j][i])
		
		
		if statistics.stdev(control_data_per_protein)==0 and statistics.stdev(condition_data_per_protein)==0 and statistics.mean(control_data_per_protein)==statistics.mean(condition_data_per_protein):
			pval=1
		else:
			if paired_flag==1:
				[z,pval]=st.ttest_rel(control_data_per_protein, condition_data_per_protein)
			else:
				[z,pval]=st.ttest_ind(control_data_per_protein, condition_data_per_protein)
			
		pvals.append(pval)
		if paired_flag==1:
			fold_per_sample=list()
				
			for k in range(len(condition_data_per_protein)):
				fold_per_sample.append(condition_data_per_protein[k]-control_data_per_protein[k])
		if statistics.mean(control_data_per_protein)==0:
			folds.append(0)
			stdevs.append(0)
		else:
			folds.append(statistics.mean(condition_data_per_protein)-statistics.mean(control_data_per_protein))
			if paired_flag==1:
				stdevs.append(statistics.stdev(fold_per_sample))
			else:
				stdevs.append(0)
	if paired_flag==0:
		stdevs=[0]*len(proteins)
	return [pvals,folds,stdevs]

def perform_pairwise_analysis(proteins, data, name_prefix,annotation,tags,parametric_flag,p_value_threshold,output_message,paired_flag):
	
	splitted_data=list()
	for i in range(len(tags)):
		temp_table=list()
		num_of_samples=0
		for k in range(len(data)):
			if annotation[k]==tags[i]:
				temp_table.append([])
				for j in range(len(data[0])):
					temp_table[num_of_samples].append(data[k][j])
				num_of_samples+=1
		splitted_data.append(temp_table)
	
	
	for i in range(len(tags)):
		for j in range(len(tags)):
			if j>i:
				if parametric_flag=='0' or parametric_flag=='2':
					output_message+='Wilcoxon test was used for pairwise comparisons\n'
					
					[pvals, folds,stdevs]=wilcoxon_rank_sum_test(proteins,splitted_data[i], splitted_data[j],paired_flag )
				else:
					[pvals, folds,stdevs]=t_test(proteins,splitted_data[i], splitted_data[j], paired_flag )
					output_message+='Paired t-test was used for pairwise comparisons\n'
				for k in range(len(pvals)):
					#print(pvals[k])
					if 'nan' in str(pvals[k]):
						
						pvals[k]=1
				pvalues=smm.multipletests(pvals,method='fdr_bh')
				pvals2=pvalues[1]
				message=''
				message+='IDs\tInitial Pvalue\tAdjusted Pvalue\tFold Change\tStandard Deviation of Fold Changes\n'				
				for k in range(len(pvals)):
					message+=proteins[k]+'\t'+str(pvals[k])+'\t'+str(pvals2[k])+'\t'+str(folds[k])+'\t'+str(stdevs[k])+'\n'
				resultFile= open(name_prefix+'all_pvals_'+str(tags[i])+'_vs_'+str(tags[j])+'.tsv','w') 
				resultFile.write(message)
				resultFile.close()
			
	return output_message



		
def cluster_proteins_combined(proteins,data, filename,name,prots,data_proteins):
	data_for_clustering=list()

	for j in range(len(proteins)):
		data_for_clustering.append([])

		splitted_data=list()
		for i in range(timepoints):
			if i>0:
				temp_table=[]
				for k in range(len(data)):
					if k%timepoints==i:
						temp_table.append(data[k][j])
				data_for_clustering[j].append(sum(temp_table)/float(len(temp_table)))
		splitted_data.append(temp_table)
	for j in range(len(prots)):
		data_for_clustering.append([])
	
		splitted_data=list()
		for i in range(timepoints):
			if i>0:
				temp_table=[]
				for k in range(len(data_proteins)):
					if k%timepoints==i:
						temp_table.append(data_proteins[k][j])
				data_for_clustering[j+len(proteins)].append(sum(temp_table)/float(len(temp_table)))
		splitted_data.append(temp_table)

	db = DBSCAN(eps=0.5, min_samples=2).fit(data_for_clustering)
	core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
	core_samples_mask[db.core_sample_indices_] = True
	labels = db.labels_

	# Number of clusters in labels, ignoring noise if present.
	n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

	print('Estimated number of clusters: %d' % n_clusters_)
	
	print("Silhouette Coefficient: %0.3f"
		  % metrics.silhouette_score(data_for_clustering, labels))


	

	# Black removed and is used for noise instead.
	unique_labels = set(labels)
	print(unique_labels)
	colors = [plt.cm.Spectral(each)
			  for each in np.linspace(0, 1, len(unique_labels))]
	resultFile= open(filename,'w') 
	message=''
	for k, col in zip(unique_labels, colors):
		

		class_member_mask = (labels == k)
		data_for_clustering_filtered=list()
		mpl_fig = plt.figure()
		ax = mpl_fig.add_subplot(111)
		protein_list=''
		for iter in range(len(data_for_clustering)):
			#print(labels)
			#print(iter)
			if labels[iter]==k:
				#data_for_clustering_filtered.append(copy.deepcopy(data_for_clustering_filtered))
				if (iter>=len(proteins)):
					protein_list=protein_list+'\t'+prots[iter-len(proteins)]
				else:
					protein_list=protein_list+'\t'+proteins[iter]
				ax.plot(data_for_clustering[iter])
		
		
		
		
		ax.set_title('Cluster_'+str(k))
		ax.set_ylabel('Average Folded Concetrations')
		ax.set_xlabel('Time-points')
		#labels_n=list()
		labels_n = [item.get_text() for item in ax.get_xticklabels()]
		
		labels_n [1] = '0'
		labels_n [3] ='60min'
		labels_n [5] = '8h'
		labels_n [7] ='24h'
		ax.set_xticklabels(labels_n)
		#mpl_fig.Title(proteins[j])
		mpl_fig.savefig(name+'_Combined_Cluster_'+str(k)+'_figure.png')
		mpl_fig.clf()
		
		message+='Cluster_'+str(k)+':'+protein_list+'\n'
		#
	resultFile.write(message)
	resultFile.close()	
	
	
		
def create_combined_data(proteins,data, filename,name,prots,data_proteins):
	data_for_clustering=list()#
	variables=list()
	#print(data)
	for j in range(len(proteins)):
		variables.append(proteins[j])
		data_for_clustering.append([])
		#data_for_clustering[j].append(proteins[j])
		splitted_data=list()
		for i in range(timepoints):
			if i>0:
				temp_table=[]
				for k in range(len(data)):
					if k%timepoints==i:
						data_for_clustering[j].append(data[k][j])
				
	for j in range(len(prots)):
		variables.append(prots[j])
		data_for_clustering.append([])
		#data_for_clustering[j].append(proteins[j])
		splitted_data=list()
		for i in range(timepoints):
			if i>0:
				temp_table=[]
				for k in range(len(data_proteins)):
					if k%timepoints==i:
						data_for_clustering[j+len(proteins)].append(data_proteins[k][j])
		
	return [variables,data_for_clustering]
	
	

		
def create_combined_data_plasma(proteins,data, filename,name,prots,data_proteins):
	data_for_clustering=list()#
	variables=list()
	#print(data)
	
				
	for j in range(len(prots)):
		variables.append(prots[j])
		data_for_clustering.append([])
		#data_for_clustering[j].append(proteins[j])
		splitted_data=list()
		for i in range(timepoints):
			if i>0:
				temp_table=[]
				for k in range(len(data_proteins)):
					if k%timepoints==i:
						data_for_clustering[j].append(data_proteins[k][j])
		
	return [variables,data_for_clustering]
	
def spearman_correlation(markers,data,filesuffix):
	
	spearman_table= open(filesuffix+'spearman_table.txt','w')
	spearman_network=open(filesuffix+'spearman_network.txt','w')
	num_of_pairs=0
	for i in range(len(markers)):
		spearman_table.write('\t'+markers[i])
	spearman_table.write('\n')
	for i in range(len(markers)):
		spearman_table.write(markers[i])
		for j in range(len(markers)):
			if j>i:
				
				[r,p]=st.spearmanr(data[i],data[j])
				spearman_table.write('\t'+str(r))
				if p<0.05 and r>0.5:
					spearman_network.write(markers[i]+'\t'+markers[j]+'\t'+str(r)+'\n')
			elif j==i:
				spearman_table.write('\t1')
			else:
				spearman_table.write('\t')
		spearman_table.write('\n')
					
def filter_dataset(dataset_initial,proteins,percentage,output_message):
	new_data=list()
	selected=0
	new_proteins=list()
	missing_proteins=0
	proteins_missing_values_percentage=0
	
	
	for i in range(len(dataset_initial)):
		missing=0
		
		for j in range (len(dataset_initial[0])):
			if (dataset_initial[i][j]=='') or dataset_initial[i][j]==-1000:
				missing+=1
				proteins_missing_values_percentage+=1
					
		
		if ((missing/float(len(dataset_initial[0])))<=(percentage)):
			#print(i)
			new_data.append([])
			for k in range(len(dataset_initial[i])):
				new_data[selected].append(dataset_initial[i][k])
			selected+=1
			new_proteins.append(proteins[i])
		else:
			
			missing_proteins+=1
	print('Data were successfully filtered!')
	print('Total Number of Molecules='+str(len(dataset_initial)))
	print('Total Number of Molecules with missing values less than less than allowed threshold='+str(selected))
	print('Percentage of Missing Values in all molecules='+str(proteins_missing_values_percentage/float(len(dataset_initial)*len(dataset_initial[0]))))
	output_message+='Total Number of Molecules='+str(len(dataset_initial))+'\n'
	output_message+='Total Number of Molecules with missing values less than less than allowed threshold='+str(selected)+'\n'
	output_message+='Percentage of Missing Values in all molecules='+str(proteins_missing_values_percentage/float(len(dataset_initial)*len(dataset_initial[0])))+'\n'
	return [new_data,new_proteins,output_message]

	
def perform_missing_values_imputation(dataset_initial,missing_imputation_method, folder_name,output_message):
	#missing values imputation
	averages=[0]*(len(dataset_initial))
	if missing_imputation_method=="1":
		#average imputation here
		num_of_non_missing_values=[0]*(len(dataset_initial))
		for j in range(len(dataset_initial)):
			for i in range((len(dataset_initial[0]))):
				if dataset_initial[j][i]!=-1000 and dataset_initial[j][i]!='':
					#print(dataset_initial[j][i+1])
					averages[j]+=float(dataset_initial[j][i])
					num_of_non_missing_values[j]+=1
			
			averages[j]=averages[j]/float(num_of_non_missing_values[j])
		write_one_dimensional_list_to_tab_delimited_file(averages, folder_name+'averages_for_missing_values_imputation.csv')
		output_message+='Average imputation method was used!\n'
		for i in range(len(dataset_initial)):
			for j in range((len(dataset_initial[0]))):
				if dataset_initial[i][j]==-1000:
					dataset_initial[i][j]=averages[i]
		return dataset_initial
	else:
		#KNN-impute 
		dataset_initial=list(map(list, zip(*dataset_initial)))
		for i in range(len(dataset_initial)):
			for j in range(len(dataset_initial[0])):
				if dataset_initial[i][j]=='' or dataset_initial[i][j]==-1000:
					dataset_initial[i][j]= np.NaN
		dataset = knn_impute_optimistic(np.asarray(dataset_initial), np.isnan(np.asarray(dataset_initial)), k=3)
		dataset=list(map(list, zip(*dataset)))
		output_message+='KNN imputation method was used!\n'
		return [dataset,output_message]
	

def normalize_dataset(dataset_initial, folder_name,output_message,normalization_method):	
	
	
	if normalization_method=='1':
		#arithmetic sample-wise normalization
		maximums=[-1000.0]*(len(dataset_initial[0]))
		minimums=[1000.0]*(len(dataset_initial[0]))
		for i in range(len(dataset_initial)):
			for j in range((len(dataset_initial[0]))):
				if(dataset_initial[i][j]!="" or dataset_initial[i][j]!=-1000):
					if float(dataset_initial[i][j])>maximums[j]:
						maximums[j]=float(dataset_initial[i][j])
					if float(dataset_initial[i][j])<minimums[j]:
						minimums[j]=float(dataset_initial[i][j])
		max1=max(maximums)
		min1=min(minimums)
		print('Maximum Quantity Value:'+str(max1))
		print('Minimum Quantity Value:'+str(min1))
		for i in range(len(dataset_initial)):
			for j in range((len(dataset_initial[0]))):
				if (dataset_initial[i][j]!="" or dataset_initial[i][j]!=-1000):
					dataset_initial[i][j]=0+(1/(max1-min1))*(float(dataset_initial[i][j])-min1)
		output_message+='Arithmetic normalization was used!\n'
		return [dataset_initial,output_message]
	else:
		logged_data=list()
	
		for i in range(len(dataset_initial)):
			#print('i='+str(i))
			logged_data.append([])
			for j in range(len(dataset_initial[0])):
				#print('j='+str(j))
				if dataset_initial[i][j]=='' or dataset_initial[i][j]==-1000:
					logged_data[i].append('')
				else:
					if(dataset_initial[i][j]==0):
						logged_data[i].append(0)
					else:
						logged_data[i].append(math.log2(dataset_initial[i][j]))
		output_message+='Logarithmic normalization was used!\n'
		return[logged_data,output_message]

def average_duplicate_measurements(dataset_initial,markers):
	dataset={}
	dict_of_occurences={}
	num_of_elements=0
	for i in range(len(dataset_initial)):
		if dataset_initial[i][0] not in dataset:
			dict_of_occurences[markers[i]]=1
			dataset[markers[i]]=list()
			for j in range(len(dataset_initial[0])):
				dataset[markers[i]].append(float(dataset_initial[i][j]))
		else:
			
			dict_of_occurences[markers[i]]+=1
			for j in range(len(dataset_initial[0])):
				dataset[markers[i]][j]=dataset[markers[i]][j]+float(dataset_initial[i][j])
		num_of_elements+=1
	element=0
	for key in dataset:
		for j in range(len(dataset[key])):
			dataset[key][j]=dataset[key][j]/dict_of_occurences[key]
	data=list()
	markers=list()
	num_of_markers=0
	for key,vals in dataset.items():
		data.append([])
		markers.append(key)
		for i in range(len(vals)):
			data[num_of_markers].append(vals[i])
		num_of_markers+=1
	return [data,markers]
	
def print_data(data,markers,tags,folder_name, filename,number_of_phenotypes):
	file=open(folder_name+filename,'w')
	message=''
	for i in range(len(data[0])):
		message=message+'\t'+tags[i%number_of_phenotypes]
	message+='\n'
	for i in range(len(data)):
		message+=markers[i]
		for j in range(len(data[0])):
			message+='\t'+str(data[i][j])
		message+='\n'
	file.write(message)
	file.close()

def print_heatmap(proteins,data,labels,unique_labels, samples,scaling, filename):
	samples_colors=list()
	col_palet=list()
	col_palet.append("#0000FF")
	col_palet.append("#FF0000")
	col_palet.append("#00FF00")
	col_palet.append("#FFFF00")
	col_palet.append("#00FFFF")
	
	if len(unique_labels)<=5:
		max_num_of_colors=len(unique_labels)
	else:
		max_num_of_colors=5
	for i in range(len(data[0])):
		for k in range(max_num_of_colors):
			if labels[i]==unique_labels[k]:
				samples_colors.append(col_palet[k])
		
		
	
	print('Heatmaps as PNG')
	r.png(filename, width = 10, height = 10,units = 'in',res = 1200)
	data_array=np.asarray(data)
	
	#r.heatmap(data_array)
	print(len(data_array))
	print(len(np.asarray(samples_colors)))
	r.heatmap(data_array, cexRow=scaling, labCol=np.asarray(samples),labRow=np.asarray(proteins),ColSideColors = np.asarray(samples_colors), col = ro.r.redgreen(75), show_heatmap_legend = True, show_annotation_legend = True)
	r("dev.off()")	
	
if __name__ == "__main__":
	#python statistical_data_analysis_multiple_conditions_v3_0.py example_dataset.txt example_labels.txt example_samples.txt 0 2 2 0.3 0.05 0 example_results/

	

	data_filename=sys.argv[1]
	tags_filename=sys.argv[2]
	samples_filename=sys.argv[3]
	parametric_flag=sys.argv[4]
	normalization_method=sys.argv[5]
	missing_imputation_method=sys.argv[6]
	missing_threshold=float(sys.argv[7])
	p_value_threshold=float(sys.argv[8])
	paired_flag=int(sys.argv[9])
	folder_name=sys.argv[10]
	
	output_message=''
	if not os.path.exists(folder_name):
		os.makedirs(folder_name)
	labels=list()
	with open(tags_filename) as tags_fname:
		for line in tags_fname:
			words=line.split('\t')
			for label in words:
				labels.append(label.strip())
		
	print(labels)
	samples=list()
	with open(samples_filename) as samples_fname:
		for line in samples_fname:
			words=line.split('\t')
			for sample in words:
				samples.append(sample.strip())
		
	
	[markers,data]=parse_data(data_filename)
	#data has number_of_samples * timepoints rows and number of markrs columns
	
	data_transp=list(map(list, zip(*data)))
	
	[data_transp,markers,output_message]=filter_dataset(data_transp,markers,missing_threshold,output_message)
	
	if missing_imputation_method!="0":
		[data_transp,output_message]=perform_missing_values_imputation(data_transp,missing_imputation_method, folder_name,output_message) 
		
	if normalization_method!="0":
		[data_transp,output_message]=normalize_dataset(data_transp,folder_name,output_message,normalization_method)
		
	[data_transp,markers]=average_duplicate_measurements(data_transp,markers)	
	
	#print_data(data_transp,markers,tags,folder_name, 'preprocessed_data.csv',timepoints)
	data=list(map(list, zip(*data_transp)))
	tags=list(set(labels))
	
	
	
	#Step 6 Print Box Plots
	plotly.tools.set_credentials_file(username='theofilk', api_key='FvCv1RgpeQC6ng3rxNfR')
	#print_boxplots(markers, folded_data, timepoints, folder_name+'boxplots/folded/', 1,tags)
	
	print_boxplots(markers, data, labels, folder_name+'boxplots/raw/', 0,tags)
	
	if parametric_flag=='0':
		data_per_phenotype=list()
		
		shapiro_pvals=list()
		flag=0
		for k in range(len(tags)):
			data_per_phenotype.append([])
			for i in range(len(data_transp)):
				for j in range(len(data_transp[i])):
					if labels[j]==tags[k]:
						data_per_phenotype[k].append(data_transp[i][j])
			
			[shapiro,shapiro_pvals_temp]=st.shapiro(data_per_phenotype[k])
			shapiro_pvals.append(shapiro_pvals_temp)
			if shapiro_pvals_temp<0.05:
				flag=1
		if flag==1:
			[Hs,pvals] = kruskal_wallis_test(data, labels,markers)
		else:
			[Hs,pvals] = anova_test(data, labels,markers)
		if flag==1:
			output_file=open(folder_name+'kruskal_wallis_test.csv','w')
		else:
			output_file=open(folder_name+'anova_test.csv','w')
		output_file.write('ID\tPvalue\tAdjusted Pvalues\n')
		pvalues=smm.multipletests(pvals,method='fdr_bh')
		pvals2=pvalues[1]
		for i in range(len(pvals)):
			output_file.write(markers[i]+'\t'+str(pvals[i])+'\t'+str(pvals2[i])+'\n')
		output_file.close()
		if flag==1:
			output_message+='Kruskal Wallis test was used for multple phenotypes comparison because at least for one of the categories data are not normally distributed\n'
			parametric_flag=2
		else:
			output_message+='Anova test was used for multple phenotypes comparison because at least for one of the categories data are not normally distributed\n'
			parametric_flag=1
	elif parametric_flag=='1':
		[Hs,pvals] = anova_test(data, labels,markers)
		
		
		output_file=open(folder_name+'anova_test.csv','w')
		output_file.write('ID\tPvalue\n')
		for i in range(len(pvals)):
			output_file.write(markers[i]+'\t'+str(pvals[i])+'\n')
		output_file.close()
		output_message+='Anova test was used for multiple phenotypes comparison\n'
	else:
		[Hs,pvals] = kruskal_wallis_test(data, labels,markers)
		output_file=open(folder_name+'kruskal_wallis_test.csv','w')
		output_file.write('ID\tPvalue\tAdjusted Pvalues\n')
		pvalues=smm.multipletests(pvals,method='fdr_bh')
		pvals2=pvalues[1]
		for i in range(len(pvals)):
			output_file.write(markers[i]+'\t'+str(pvals[i])+'\t'+str(pvals2[i])+'\n')
		output_file.close()
		output_message+='Kruskal Wallis test was used for multple phenotypes comparison\n'
	
	#Step 8 Keep Significant Data Only
	#print_significant_boxplots(markers, folded_data, timepoints, folder_name+'boxplots/significant/folded/', 1, pvals,tags,p_value_threshold)
	print_significant_boxplots(markers, data, labels, folder_name+'boxplots/significant/raw/', 0, pvals,tags,p_value_threshold)
	
	
	print_heatmap(markers,data_transp,labels,tags, samples,0.4, folder_name+'heatmap_all.png')
	
	#Step 11 Perform Pairwise Statistical Analysis
	output_message=perform_pairwise_analysis(markers, data, folder_name,labels,tags,parametric_flag,p_value_threshold,output_message,paired_flag)
	
	out_file=open(folder_name+'info.txt','w')
	out_file.write(output_message)
	out_file.close()
	

	
	