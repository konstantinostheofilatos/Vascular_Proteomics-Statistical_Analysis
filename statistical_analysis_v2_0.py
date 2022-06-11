
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
import plotly.tools as tls
import statsmodels.stats.multitest as smm
import plotly
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import mstats
import copy
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler
from itertools import compress
from sklearn.decomposition import PCA
from sklearn.neighbors import LocalOutlierFactor

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
from sklearn.preprocessing import Imputer
from knnimpute import (
    knn_impute_few_observed,
    knn_impute_with_argpartition,
    knn_impute_optimistic,
    knn_impute_reference,
)

#import rpy2.robjects.lib.limma as limma


def outlier_detection(dataset,folder_name):
	pca = PCA(n_components=0.9,svd_solver = 'full')
	new_dataset=list()
	num_of_samples=0
	for j in range(len(dataset[0])):
		new_dataset.append([])
		for i in range(len(dataset)):
			new_dataset[num_of_samples].append(float(dataset[i][j]))
		num_of_samples+=1
	dataset_new=pca.fit_transform(new_dataset)
	print(len(dataset_new))
	print(len(dataset_new[0]))
	clf = LocalOutlierFactor(n_neighbors=20)
	y_pred = clf.fit_predict(dataset_new)
	write_one_dimensional_list_to_tab_delimited_file(y_pred,folder_name+'outlier_prediction.txt')
	return dataset_new
	


def parse_data(data_filename):
	num_of_lines=0
	proteins=list()
	data=list()
	samples=list()
	with open(data_filename) as data_fname:
		for line in csv.reader(data_fname, delimiter="\t"):
			if num_of_lines==0:
				for j in range(len(line)):
					if j>0:
						samples.append(line[j].strip())
			else:
				proteins.append(line[0])
				data.append([])
				for j in range(len(line)):
					if j>0:
						if line[j]!='':
							data[num_of_lines-1].append(float(line[j]))
						else:
							data[num_of_lines-1].append('')
			num_of_lines+=1
	print('Data were successfully parsed!')
	return [proteins,data,samples]
	
def parse_labels(labels_filename):
	symptomatic_asymptomatic=list()
	diabetes=list()
	num_of_lines=0
	with open(labels_filename) as labels_fname:
		for line in csv.reader(labels_fname, delimiter="\t"):
			if num_of_lines==0:
				for i in range(len(line)):
					if line[i].strip()=='Symptomatic':
						symptomatic_asymptomatic.append('Symptomatic')
					else:
						symptomatic_asymptomatic.append('Asymptomatic')
			else:
				for i in range(len(line)):
					diabetes.append((line[i].strip()))
			num_of_lines+=1
	print('Labels were successfully parsed!')
	return [symptomatic_asymptomatic, diabetes]
	
def parse_labels_new(labels_filename):
	labels=list()
	num_of_lines=0
	with open(labels_filename) as labels_fname:
		for line in csv.reader(labels_fname, delimiter="\t"):
			
			for i in range(len(line)):
				labels.append(line[i].strip())
			

	print('Labels were successfully parsed!')
	return labels
	
def parse_commorbidities(samples,commorbidities_filename):
	age=list()
	sex=list()
	statin=list()
	a_or_b=list()
	commorbidities=dict()
	commorbidities_flag=0
	num_of_lines=0
	commorbidities=list()
	with open(commorbidities_filename) as commorbidities_fname:
		for line in csv.reader(commorbidities_fname, delimiter="\t"):
			commorbidities_flag=1
			if len(commorbidities)==0:
				for i in range(len(line)):
					commorbidities.append([])
			for i in range(len(line)):
				commorbidities[i].append(line[i].strip())
	return commorbidities, commorbidities_flag
	



def filter_proteomics_dataset(dataset_initial,proteins,percentage,output_message):
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
					
		
		if ((missing/float(len(dataset_initial[0])))<(percentage)):
			#print(i)
			new_data.append([])
			for k in range(len(dataset_initial[i])):
				new_data[selected].append(dataset_initial[i][k])
			selected+=1
			new_proteins.append(proteins[i])
		else:
			
			missing_proteins+=1
	print('Data were successfully filtered!')
	print('Total Number of Proteins='+str(len(dataset_initial)))
	print('Total Number of Proteins with missing values less than predefined threshold='+str(selected))
	print('Percentage of Missing Values in all Proteins='+str(proteins_missing_values_percentage/float(len(dataset_initial)*len(dataset_initial[0]))))
	output_message+='Total Number of Molecules='+str(len(dataset_initial))+'\n'
	output_message+='Total Number of Molecules with missing values less than less than allowed threshold='+str(selected)+'\n'
	output_message+='Percentage of Missing Values in all molecules='+str(proteins_missing_values_percentage/float(len(dataset_initial)*len(dataset_initial[0])))+'\n'
	return [new_data,new_proteins,output_message]

def f(dataset_initial,missing_imputation_method, folder_name,output_message):
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
	
def write_list_to_tab_delimited_file(data,filename):
	file_id=open(filename,'w')
	for i in range(len(data)):
		for j in range(len(data[0])):
			file_id.write(str(data[i][j]))
			if j!=len(data[0])-1:
				file_id.write('\t')
			else:
				file_id.write('\n')
	file_id.close()			
	
def write_one_dimensional_list_to_tab_delimited_file(data,filename):
	file_id=open(filename,'w')
	for i in range(len(data)):
		file_id.write(str(data[i]))
		file_id.write('\n')
	file_id.close()
	
def differential_expression_analysis(proteins, data, age, sex, statin, labels, mylist,suffix):
	proteins_vector=ro.StrVector(proteins)
	
	print(mylist)
	
					
	data_array=np.asarray(data)
	
	data_robject=ro.r.matrix(data_array, nrow=len(data_array))

	data_robject.rownames=proteins_vector
	labels_array=np.asarray(labels)
	labels_robject=ro.r.matrix(labels_array,nrow=len(labels_array))
	commorbidity_robject1=ro.r.matrix(np.asarray(age),nrow=len(labels_array))
	commorbidity_robject2=ro.r.matrix(np.asarray(sex),nrow=len(labels_array))
	commorbidity_robject3=ro.r.matrix(np.asarray(statin),nrow=len(labels_array))
	#print(labels_robject)
	#design=model.matrix(ro.FactorVector(labels_robject))
	#print(design)
	ro.r.assign('data',data_robject)
	ro.r.assign('labels', labels_robject)
	ro.r.assign('age',commorbidity_robject1)
	ro.r.assign('sex',commorbidity_robject2)
	ro.r.assign('statin',commorbidity_robject3)
	
	ro.r('design<-model.matrix(~factor(labels)+age+factor(sex)+factor(statin))')
	
	
	
	
	ro.r('fit <-lmFit (data, design)')
	
	ro.r('fit2<-eBayes(fit)')
	fit2=ro.r('eBayes(fit)')
	table=limma.topTable(fit2, coef=2,  n=len(data))
	columns=list(np.asarray(table.colnames))
	proteins_ids=list(np.asarray(table.rownames))
	table_array=(np.asarray(table))
	
	output_file=open('diff_exp_results'+suffix+'.txt','w')
	output_file.write('\t')
	for title in columns:
		output_file.write(str(title)+'\t')
	output_file.write('\n')
	for i in range(len(table_array[0])):
		output_file.write(str(proteins_ids[i])+'\t')
		for j in range(len(table_array)):
			output_file.write(str(table_array[j][i])+'\t')
		output_file.write('\n')
	output_file.close()
	
	return [table, proteins_ids, columns]
	
	
	
def differential_expression_analysis_new(proteins, data, labels, commorbidities_flag, commorbidities,commorbidities_type, mylist,folder_name,suffix):
	proteins_vector=ro.StrVector(proteins)
	
	print(mylist)
	
					
	data_array=np.asarray(data)
	
	data_robject=ro.r.matrix(data_array, nrow=len(data_array))

	data_robject.rownames=proteins_vector
	labels_array=np.asarray(labels)
	labels_robject=ro.r.matrix(labels_array,nrow=len(labels_array))
	
	#print(labels_robject)
	#design=model.matrix(ro.FactorVector(labels_robject))
	#print(design)
	ro.r.assign('data',data_robject)
	ro.r.assign('labels', labels_robject)
	if commorbidities_flag==0:
		ro.r('design<-model.matrix(~factor(labels))')
	else:
		command='design<-model.matrix(~factor(labels)'
		for i in range( len(commorbidities)):
			ro.r.assign('commorbidity_'+str(i),ro.r.matrix(np.asarray(commorbidities[i]),nrow=len(labels_array)))
			if commorbidities_type[i]=='0':
				command+='+factor(commorbidity_'+str(i)+')'
			else:
				command+='+commorbidity_'+str(i)
		command+=')'
		ro.r(command)
	
	ro.r('fit <-lmFit (data, design)')
	
	ro.r('fit2<-eBayes(fit)')
	fit2=ro.r('eBayes(fit)')
	table=limma.topTable(fit2, coef=2,  n=len(data))
	columns=list(np.asarray(table.colnames))
	proteins_ids=list(np.asarray(table.rownames))
	table_array=(np.asarray(table))
	
	output_file=open(folder_name+'diff_exp_results'+suffix+'.txt','w')
	output_file.write('\t')
	for title in columns:
		output_file.write(str(title)+'\t')
	output_file.write('\n')
	for i in range(len(table_array[0])):
		output_file.write(str(proteins_ids[i])+'\t')
		for j in range(len(table_array)):
			output_file.write(str(table_array[j][i])+'\t')
		output_file.write('\n')
	output_file.close()
	
	return [table, proteins_ids, columns]

	
def print_volcano_plots(proteins, pvals, log_fold_changes, filename,filename2,p_value_threshold):
	thresholds=list()
	selected_proteins=list()
	for i in range (len(pvals)):
		if pvals[i] <p_value_threshold :
			thresholds.append(1)
			selected_proteins.append(proteins[i])
		else:
			selected_proteins.append('')
			thresholds.append(0)
	ro.FactorVector(thresholds)
	log_pvals=list()
	colors=list()
	for i in range(len(pvals)):
		if pvals[i]<0.001:
			colors.append('red2')
		elif pvals[i]<0.01:
			colors.append('orange1')
		elif pvals[i]<0.05:
			colors.append('darkgreen')
		else:
			colors.append('darkblue')
		log_pvals.append(-1*math.log10(pvals[i]))
	colormap_raw = [['red2', '#ff0000'],['orange', '#FFA500'],['darkgreen', '#006400'],['darkblue', '#003366'] ]
	colormap_labels = [['red2', 'P < 0.001'],['orange', 'P < 0.01'], ['darkgreen', 'P < 0.05'], ['darkblue', 'P > 0.05']]
	colormap = ro.StrVector([elt[1] for elt in colormap_raw])
	colormap.names = ro.StrVector([elt[0] for elt in colormap_raw])	
	ro.r('p1<-expression(paste(-log[10], \"(p-value)\" ))')
	ro.r('p2<-expression(paste(log[2], \"(FC)\" ))')
	df_dict = {'Ids':ro.StrVector(selected_proteins),'threshold':ro.IntVector(thresholds),'log2FC': ro.FloatVector(log_fold_changes), 'MinusLog10Pvals': ro.FloatVector(log_pvals),'colors':ro.StrVector(colors)}
	
	df= ro.DataFrame(df_dict)
	
	r.png(filename)
	#gp=ggplot2.ggplot(data=df) +ggplot2.aes_string(x='log2FC', y='MinusLog10Pvals', label='Ids', colour='colors') +ggplot2.geom_point() +ggplot2.geom_text(colour ='black', check_overlap='True')+ggplot2.scale_colour_manual("", 
    #                              values=colormap,
    #                             breaks=colormap.names,
    #                              labels=[elt[1] for elt in 
    #                                      colormap_labels]) 
	gp=ggplot2.ggplot(data=df) +ggplot2.aes_string(x='log2FC', y='MinusLog10Pvals', label='Ids', colour='colors') +ggplot2.geom_point() +ggplot2.scale_colour_manual("", 
                                  values=colormap,
                                  breaks=colormap.names,
                                  labels=[elt[1] for elt in 
                                          colormap_labels]) 
	gp.plot()
	
	ro.r('p1<-expression(paste(-log[10], \"(p-value)\" ))')
	ro.r('p2<-expression(paste(log[2], \"(FC)\" ))')
	
	r("dev.off()")
	
	
	
def print_boxplots(new_proteins,proteins, data, values_list, labels_1, folder, pvals,p_value_threshold):
	
		
	for j in range(len(proteins)):
		num_of_cat=list()
		num_of_cat.append(0)
		num_of_cat.append(0)
		if pvals[j]<p_value_threshold:
			
			splitted_data=list()
			for i in range(2):
				
				temp_table=[]
				
				for k in range(len(labels_1)):
						
					if labels_1[k]==values_list[i]:
						
						ind=new_proteins.index(proteins[j])
						if data[ind][k]!='' and data[ind][k]!=-1000:
							num_of_cat[i]+=1
							temp_table.append(data[ind][k])
				splitted_data.append(temp_table)	
			mpl_fig = plt.figure()
			mpl_fig.subplots_adjust(bottom=0.15)
			ax = mpl_fig.add_subplot(111)

			ax.boxplot(splitted_data)
			
			ax.set_title(proteins[j])
			ax.set_ylabel('Relative Quantities')
			ax.set_xlabel(values_list[0]+' vs '+values_list[1])
			labels = [item.get_text() for item in ax.get_xticklabels()]
			labels[0] = values_list[0]+'\n(n=' +str(num_of_cat[0])+')'
			labels[1] = values_list[1]+'\n(n=' +str(num_of_cat[1])+')'

			ax.set_xticklabels(labels)
			#mpl_fig.Title(proteins[j])
			if not os.path.exists(folder):
				os.makedirs(folder)
			mpl_fig.savefig(folder+proteins[j]+'_boxplot.png')
			mpl_fig.clf()
			plt.close('all')
			#plotly_fig = tls.mpl_to_plotly( mpl_fig )
			#plot_url = py.plot(plotly_fig, 'mpl-multiple-boxplot')	

def average_duplicate_measurements(dataset_initial,markers):
	dataset={}
	dict_of_occurences={}
	num_of_elements=0
	for i in range(len(dataset_initial)):
		if dataset_initial[i][0] not in dataset:
			dict_of_occurences[markers[i]]=1
			dataset[markers[i]]=list()
			for j in range(len(dataset_initial[0])):
				if dataset_initial[i][j]!=-1000 and dataset_initial[i][j]!='':
					dataset[markers[i]].append(float(dataset_initial[i][j]))
				else:
					dataset[markers[i]].append(0)
		else:
			
			dict_of_occurences[markers[i]]+=1
			for j in range(len(dataset_initial[0])):
				if dataset_initial[i][j]!=-1000 and dataset_initial[i][j]!='':
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
		write_one_dimensional_list_to_tab_delimited_file(averages, folder_name+'averages_for_missing_values_imputation.txt')
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
	
def print_data(data,markers,labels,folder_name, filename):
	file=open(folder_name+filename,'w')
	message=''
	for i in range(len(data[0])):
		message=message+'\t'+labels[i]
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
	for i in range(len(data[0])):
		if labels[i]==unique_labels[0]:
			samples_colors.append("#FF0000") # Red
		else:
			samples_colors.append("#0000FF") # Blue
	
	
	print('Heatmaps as PNG')
	r.png(filename, width = 10, height = 10,units = 'in',res = 1200)
	data_array=np.asarray(data)
	
	#r.heatmap(data_array)
	
	r.heatmap(data_array, cexRow=scaling, labCol=np.asarray(samples),labRow=np.asarray(proteins),ColSideColors = np.asarray(samples_colors), col = ro.r.redgreen(75), show_heatmap_legend = True, show_annotation_legend = True)
	r("dev.off()")	
	

def parse_data_carotid(data_filename):
	num_of_lines=0
	proteins=list()
	data=list()
	samples=list()
	with open(data_filename) as data_fname:
		for line in csv.reader(data_fname, delimiter="\t"):
			if num_of_lines==0:
				for j in range(len(line)):
					if j>0:
						samples.append(line[j].strip())
			else:
				proteins.append(line[0])
				data.append([])
				for j in range(len(line)):
					if j>0:
						if line[j]!='':
							data[num_of_lines-1].append(float(line[j]))
						else:
							data[num_of_lines-1].append(-1000)
			num_of_lines+=1
	print('Data were successfully parsed!')
	return [proteins,data,samples]
	
def parse_labels_carotid(labels_filename):
	symptomatic_asymptomatic=list()
	diabetes=list()
	num_of_lines=0
	with open(labels_filename) as labels_fname:
		for line in csv.reader(labels_fname, delimiter="\t"):
			if num_of_lines==0:
				for i in range(len(line)):
					if line[i].strip()=='Symptomatic':
						symptomatic_asymptomatic.append('Symptomatic')
					else:
						symptomatic_asymptomatic.append('Asymptomatic')
			else:
				for i in range(len(line)):
					diabetes.append((line[i].strip()))
			num_of_lines+=1
	print('Labels were successfully parsed!')
	return [symptomatic_asymptomatic, diabetes]
	
def parse_labels_new_carotid(labels_filename):
	labels=list()
	num_of_lines=0
	with open(labels_filename) as labels_fname:
		for line in csv.reader(labels_fname, delimiter="\t"):
			
			for i in range(len(line)):
				labels.append(line[i].strip())
			

	print('Labels were successfully parsed!')
	return labels
	
def parse_commorbidities_carotid(samples,commorbidities_filename):
	age=list()
	sex=list()
	statin=list()
	a_or_b=list()
	commorbidities=dict()
	
	num_of_lines=0
	with open(commorbidities_filename) as commorbidities_fname:
		for line in csv.reader(commorbidities_fname, delimiter="\t"):
			commorbidities[line[0]]=[int(line[1].strip()), (line[2].strip()), (line[3].strip())]
			num_of_lines+=1
			
	for i in range(len(samples)):
		if samples[i][-1]=='A':
			a_or_b.append('A')
		else:
			a_or_b.append('B')
		age.append(commorbidities[samples[i][0:-1]][0])
		sex.append(commorbidities[samples[i][0:-1]][1])
		statin.append(commorbidities[samples[i][0:-1]][2])
	print('Commorbidities were successfully parsed!')
	return [age,sex,statin,a_or_b]

def cluster_samples(samples,data):
	data = StandardScaler().fit_transform(data)
	print(data)
	#data should be [samples X PCAs]
	db = DBSCAN(eps=0.3, min_samples=3).fit(data)
	core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
	core_samples_mask[db.core_sample_indices_] = True
	labels = db.labels_
	print(labels)
	print(len(labels))
	print(samples)
	# Number of clusters in labels, ignoring noise if present.
	n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

	print('Estimated number of clusters: %d' % n_clusters_)
	
	#print("Silhouette Coefficient: %0.3f"		  % metrics.silhouette_score(data_for_clustering, labels))


	return labels

	# # Black removed and is used for noise instead.
	# unique_labels = set(labels)
	# print(unique_labels)
	# colors = [plt.cm.Spectral(each)
			  # for each in np.linspace(0, 1, len(unique_labels))]
	# resultFile= open(filename,'w') 
	# message=''
	# for k, col in zip(unique_labels, colors):
		

		# class_member_mask = (labels == k)
		# data_for_clustering_filtered=list()
		# mpl_fig = plt.figure()
		# ax = mpl_fig.add_subplot(111)
		# protein_list=''
		# for iter in range(len(data_for_clustering)):
			# #print(labels)
			# #print(iter)
			# if labels[iter]==k:
				# #data_for_clustering_filtered.append(copy.deepcopy(data_for_clustering_filtered))
				# protein_list=protein_list+'\t'+proteins[iter]
				# ax.plot(data_for_clustering[iter])
		
		
		
		
		# ax.set_title('Cluster_'+str(k))
		# ax.set_ylabel('Folded Concetrations')
		# ax.set_xlabel('Time-points')
		# #labels_n=list()
		# labels_n = [item.get_text() for item in ax.get_xticklabels()]
		
		# labels_n [1] = '0'
		# labels_n [3] ='60min'
		# labels_n [5] = '8h'
		# labels_n [7] ='24h'
		# ax.set_xticklabels(labels_n)
		# #mpl_fig.Title(proteins[j])
		# mpl_fig.savefig(name+'_Cluster_'+str(k)+'_figure.png')
		# mpl_fig.clf()
		
		# message+='Cluster_'+str(k)+':'+protein_list+'\n'
		# #
	# resultFile.write(message)
	# resultFile.close()	

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
		for j in range(len(control[i])):
			control_data_per_protein.append(control[i][j])
		for j in range(len(condition[i])):
			condition_data_per_protein.append(condition[i][j])
		
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
	

if __name__ == "__main__":
	
	#python statistical_analysis_v2_0.py example_dataset.txt example_labels.txt commorbidities2.txt commorbidities_types.txt 0 2 2 0.1 0.05 0 statistical_analysis_results2/
	
	proteomics_data_filename=sys.argv[1]
	labels_filename=sys.argv[2]
	commorbidities_filename=sys.argv[3]
	commorbidities_types_filename=sys.argv[4]
	parametric_flag=sys.argv[5]
	normalization_method=sys.argv[6]
	missing_imputation_method=sys.argv[7]
	missing_threshold=float(sys.argv[8])
	p_value_threshold=float(sys.argv[9])
	paired_flag=int(sys.argv[10])
	folder_name=sys.argv[11]
	if not os.path.exists(folder_name):
		os.makedirs(folder_name)
	commorbidities_types=list()
	commorbidities_flag=0
	with open(commorbidities_types_filename) as commorbidities_types_fname:
		for line in csv.reader(commorbidities_types_fname, delimiter="\t"):
			for i in range(len(line)):
				commorbidities_types.append(line[i])
	if len(commorbidities_types)>0:
		commorbidities_flag=1
	[markers,data,samples]=parse_data(proteomics_data_filename)
	labels=parse_labels_new(labels_filename)
	unique_labels=list(set(labels))
	
	print(unique_labels)
	
	if commorbidities_flag!=0:
		[commorbidities,commorbidities_flag]=parse_commorbidities(samples,commorbidities_filename)
	else:
		commorbidities=list()
	output_message=''
	if commorbidities_flag==0:
		output_message+='No commorbidities were provided\n'
	else:
		output_message+=str(len(commorbidities))+' commorbidities were provided.\n'
	[data_filtered,markers,output_message] = filter_proteomics_dataset(data,markers,missing_threshold,output_message)
	
	dataset_imputed=data_filtered
	if missing_imputation_method!="0":
		[dataset_imputed,output_message]=perform_missing_values_imputation(data_filtered,missing_imputation_method, folder_name,output_message) 
		
	if normalization_method!="0":
		[dataset_imputed,output_message]=normalize_dataset(dataset_imputed,folder_name,output_message,normalization_method)
	if missing_imputation_method!="0":	
		[dataset_imputed,markers]=average_duplicate_measurements(dataset_imputed,markers)	
		pca_data=outlier_detection(dataset_imputed, folder_name)	
	#data_transp=list(map(list, zip(*pca_data)))
	
	
	print_data(dataset_imputed,markers,samples,folder_name, 'preprocessed_data.txt')
	
	#write_list_to_tab_delimited_file(dataset_imputed, 'SMILR_Missing_values_imputed.txt')
	#write_one_dimensional_list_to_tab_delimited_file(new_proteins, 'SMIRL_Consistent_Proteins.txt')
	
	if parametric_flag=='0':
		data_a=list()
		data_b=list()
		for i in range(len(dataset_imputed)):
			for j in range(len(dataset_imputed[0])):
				if dataset_imputed[i][j]!='' and dataset_imputed[i][j]!=1000:
					if labels[j]==unique_labels[0]:
						data_a.append(dataset_imputed[i][j])
					else:
						data_b.append(dataset_imputed[i][j])
		[shapiro,shapiro_pvals_a]=st.shapiro(data_a)
		[shapiro,shapiro_pvals_b]=st.shapiro(data_b)
		if shapiro_pvals_a>0.05 and shapiro_pvals_b>0.05:
			test_flag=1
			output_message+='Data are normally distributed. Parametric testing will be done.\n'
		else:
			test_flag=2
			output_message+='Data are not normally distributed. Non parametric testing will be done.\n'
	elif parametric_flag=='1':
		test_flag=1
		output_message+='Parametric testing will be done.\n'
	else:
		test_flag=2
		output_message+='Non Parametric testing will be done.\n'
	if test_flag==1:
		if paired_flag==1:
			counter1=0
			counter2=0
			new_com=list()
			for j in range(len(dataset_imputed[0])):
				if labels[j]==unique_labels[0]:
					new_com.append(counter1)
					counter1+=1
				else:
					new_com.append(counter2)
					counter2+=2
			commorbidities_flag=1
			commorbidities_types=list()
			commorbidities_types.append('0')
		
		[diff_table, diff_proteins, diff_columns]=differential_expression_analysis_new(markers, dataset_imputed, labels,commorbidities_flag,commorbidities,commorbidities_types, unique_labels, folder_name,unique_labels[0]+'VS'+unique_labels[1])
		print_volcano_plots(diff_proteins, diff_table[4], diff_table[0], folder_name+unique_labels[0]+'VS'+unique_labels[1]+'_volcano_plot_corrected.png',folder_name+unique_labels[0]+'VS'+unique_labels[1]+'_volcano_plot_corrected_unlabeled.png',p_value_threshold)
		print_volcano_plots(diff_proteins, diff_table[3], diff_table[0], folder_name+unique_labels[0]+'VS'+unique_labels[1]+'_volcano_plot.png',folder_name+unique_labels[0]+'VS'+unique_labels[1]+'_volcano_plot_unlabeled.png',p_value_threshold)
		print_boxplots(markers,diff_proteins, dataset_imputed, unique_labels, labels, folder_name+unique_labels[0]+'VS'+unique_labels[1]+'_boxplots_corrected/',diff_table[4],p_value_threshold)
		print_boxplots(markers,diff_proteins, dataset_imputed, unique_labels, labels, folder_name+unique_labels[0]+'VS'+unique_labels[1]+'_boxplots/',diff_table[3],p_value_threshold)
		filtered_data=list()
		position=0
		diff_proteins_corrected=list()
		for i in range(len(diff_proteins)):
			prot=diff_proteins[i]
			if diff_table[3][i]<p_value_threshold:
				diff_proteins_corrected.append(prot)
				ind=markers.index(prot)
				filtered_data.append([])
				
				for j in range (len(dataset_imputed[ind])):
					filtered_data[position].append(dataset_imputed[ind][j])
				position+=1
		if position>1:
			print_heatmap(diff_proteins_corrected,filtered_data,labels,unique_labels, samples,0.4, folder_name+'heatmap_significant_not_corrected.png')	

		filtered_data=list()
		position=0
		diff_proteins_corrected=list()
		for i in range(len(diff_proteins)):
			prot=diff_proteins[i]
			if diff_table[4][i]<p_value_threshold:
				diff_proteins_corrected.append(prot)
				ind=markers.index(prot)
				filtered_data.append([])
				
				for j in range (len(dataset_imputed[ind])):
					filtered_data[position].append(dataset_imputed[ind][j])
				position+=1
		if position>1:
			print_heatmap(diff_proteins_corrected,filtered_data,labels,unique_labels, samples,0.4, folder_name+'heatmap_significant_corrected.png')	

		print_heatmap(markers,dataset_imputed,labels,unique_labels, samples,0.4, folder_name+'heatmap_all.png')
	
	else:
		category1=list()
		category2=list()
		
		for i in range(len(dataset_imputed)):
			category1.append([])
			category2.append([])
			for j in range(len(dataset_imputed[i])):
				if labels[j]==unique_labels[0]:
					category1[i].append(dataset_imputed[i][j])
				else:
					category2[i].append(dataset_imputed[i][j])
		[pvals,folds,stdevs]=wilcoxon_rank_sum_test(markers,category1, category2, paired_flag )
		for k in range(len(pvals)):
			#print(pvals[k])
			if 'nan' in str(pvals[k]):
				pvals[k]=1
		pvalues=smm.multipletests(pvals,method='fdr_bh')
		pvals2=pvalues[1]
		message=''
		message+='IDs\tInitial Pvalue\tAdjusted Pvalue\tFold Change\tStandard Deviation of Fold Changes\n'				
		for k in range(len(pvals)):
			message+=markers[k]+'\t'+str(pvals[k])+'\t'+str(pvals2[k])+'\t'+str(folds[k])+'\t'+str(stdevs[k])+'\n'
		resultFile= open(folder_name+'all_pvals_'+str(unique_labels[0])+'_vs_'+str(unique_labels[1])+'.tsv','w') 
		resultFile.write(message)
		resultFile.close()
		print_volcano_plots(markers, pvals2, folds, folder_name+unique_labels[0]+'VS'+unique_labels[1]+'_volcano_plot_corrected.png',folder_name+unique_labels[0]+'VS'+unique_labels[1]+'_volcano_plot_corrected_unlabeled.png',p_value_threshold)
		print_volcano_plots(markers, pvals, folds, folder_name+unique_labels[0]+'VS'+unique_labels[1]+'_volcano_plot.png',folder_name+unique_labels[0]+'VS'+unique_labels[1]+'_volcano_plot_unlabeled.png',p_value_threshold)
		print_boxplots(markers,markers, dataset_imputed, unique_labels, labels, folder_name+unique_labels[0]+'VS'+unique_labels[1]+'_boxplots_corrected/',pvals2,p_value_threshold)
		print_boxplots(markers,markers, dataset_imputed, unique_labels, labels, folder_name+unique_labels[0]+'VS'+unique_labels[1]+'_boxplots/',pvals,p_value_threshold)
		filtered_data=list()
		position=0
		diff_proteins_corrected=list()
		for i in range(len(markers)):
			prot=markers[i]
			if pvals[i]<p_value_threshold:
				diff_proteins_corrected.append(prot)
				ind=markers.index(prot)
				filtered_data.append([])
				
				for j in range (len(dataset_imputed[ind])):
					filtered_data[position].append(dataset_imputed[ind][j])
				position+=1
		if position>1:
			print_heatmap(diff_proteins_corrected,filtered_data,labels,unique_labels, samples,0.4, folder_name+'heatmap_significant_not_corrected.png')	

		filtered_data=list()
		position=0
		diff_proteins_corrected=list()
		for i in range(len(markers)):
			prot=markers[i]
			if pvals2[i]<p_value_threshold:
				diff_proteins_corrected.append(prot)
				ind=markers.index(prot)
				filtered_data.append([])
				
				for j in range (len(dataset_imputed[ind])):
					filtered_data[position].append(dataset_imputed[ind][j])
				position+=1
		if position>1:
			print_heatmap(diff_proteins_corrected,filtered_data,labels,unique_labels, samples,0.4, folder_name+'heatmap_significant_corrected.png')	

		print_heatmap(markers,dataset_imputed,labels,unique_labels, samples,0.4, folder_name+'heatmap_all.png')

	out_file=open(folder_name+'info.txt','w')
	out_file.write(output_message)
	out_file.close()
	
	