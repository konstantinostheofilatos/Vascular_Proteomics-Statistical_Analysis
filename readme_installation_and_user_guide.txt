Installation Details.

These statistical analysis scripts have been implemented in python 3.6. Thus, for their utilization they require this version of python to be installed. If they are to be used outside the python's folder, then the python.exe file should be added to path.

Additionally these scripts require the following python libraries to be installed (they can be installed with simply running the command >pip install library_name subjective that the command pip has been properly installed.):

-numpy
-scipy
-matplotlib
-plotly
-statsmodels
-sklearn
-knnimpute
-rpy2

Furthermore this installation also calls several R scripts. It has been tested with R version 3.4.2 so it is strongly suggested that this version is installed, and R executable added to path.

Several R libraries should be firstly installed before running the script including:
-gplots
-limma
-statmod
-ggrepel
-lattice
-beanplots

Moreover, within the codes before running it the users should put the correct paths in the following lines:
-os.environ['R_HOME'] = 'C:/Program Files/R/R-3.4.2'
-os.environ['R_USER'] = 'C:/Users/Konstantinos/AppData/Local/Programs/Python/Python36/Lib/site-packages/rpy2'

User Guides.
Two scripts are provided for simple statistical analysis. The first one (statistical_analysis_v2_0.py) is for comparing in between two different phenotypes (e.g. control vs disease) while the second (statistical_data_analysis_multiple_conditions_v3_0.py) is for analyzing data with more than 2 phenotypes/conditions.

statistical_data_analysis_multiple_conditions_v3_0.py:
To run this script from the command line users should execute the following command:
python statistical_data_analysis_multiple_conditions_v3_0.py [argv1] [argv2] [argv3] [argv4] [argv5] [argv6] [argv7] [argv8] [argv9] [argv10]
-[argv1] The name of the data file. It should be a tab delimited txt file with the first line corresponding to the markers name (e.g. gene ids, uniprot ids etc) and all other lines having the quantification values per sample
-[argv2] The name of a tab delimited file which will include the phenotypic description of the samples. One row of data with a number of columns equal to the number of samples
-[argv3] The name of a tab delimited file which will include the sample ids. One row of data with a number of columns equal to the number of samples.
-[argv4] 0 if you do not know if parametric or non-parametric analysis should be applied, 1 for parametric analysis and 2 for non parametric analysis
-[argv5] 0 for not performing normalization, 1 for performing arithmetic scaling to 0-1 interval and 2 for performing logarithmic normalization
-[argv6] 0 for not performing any imputation, 1 for performing imputation with the average value of the marker, and 2 for performing missing value imputation using the KNN-Impute algorithm
-[argv7] A float number indicating the maximum accepted percentage of missing values for a marker to be considered in the analysis. (e.g. 0.3 for 30% threshold)
-[argv8] A p-value threshold for considering considering a marker statistically significantly differentiated
-[argv9] 0 when unpaired analysis should be done and 1 for paired analysis
-[argv10] The name of a folder where the results will be stored

An example for running this script with the example data stored in the folder would be:
python statistical_data_analysis_multiple_conditions_v3_0.py example_dataset.txt example_labels.txt example_samples.txt 0 2 2 0.3 0.05 0 example_results/

If multiple versions of python are installed the python should be replaced with python3.6 in the above command to force this script to be executed with python 3.6 version.

statistical_analysis_v3_0.py:
To run this script from the command line users should execute the following command:
python statistical_analysis_v3_0.py [argv1] [argv2] [argv3] [argv4] [argv5] [argv6] [argv7] [argv8] [argv9] [argv10]
-[argv1] The name of the data file. It should be a tab delimited txt file with the first line corresponding to samples identifiers and the first column to the markers name (e.g. gene ids, uniprot ids etc). Every line correspond to the quantification values of a marker for all samples.
-[argv2] The name of a tab delimited file which will include the phenotypic description of the samples. One row of data with a number of columns equal to the number of samples
-[argv3] The name of a tab delimited file which will include the commorbidities information. One row of data will correspond to the values of one comorbidity over all samples.
-[argv4] The name of a tab delimited file which will include information about the types of comorbidities. One line with n values are expected with n beingthe number of comorbidities. These values should be 0 if this comorbidity will include discrete values (such as YES or NO) and 1 if it will inclue continues ordinal values(such as age). 
-[argv5] 0 if you do not know if parametric or non-parametric analysis should be applied, 1 for parametric analysis and 2 for non parametric analysis
-[argv6] 0 for not performing normalization, 1 for performing arithmetic scaling to 0-1 interval and 2 for performing logarithmic normalization
-[argv7] 0 for not performing any imputation, 1 for performing imputation with the average value of the marker, and 2 for performing missing value imputation using the KNN-Impute algorithm
-[argv8] A float number indicating the maximum accepted percentage of missing values for a marker to be considered in the analysis. (e.g. 0.3 for 30% threshold)
-[argv9] A p-value threshold for considering considering a marker statistically significantly differentiated
-[argv10] 0 when unpaired analysis should be done and 1 for paired analysis
-[argv11] The name of a folder where the results will be stored

An example for running this script with the example data stored in the folder would be:
python statistical_analysis_v2_0.py example_dataset.txt example_labels.txt commorbidities2.txt commorbidities_types.txt 0 2 2 0.1 0.05 0 statistical_analysis_results2/

If multiple versions of python are installed the python should be replaced with python3.6 in the above command to force this script to be executed with python 3.6 version.


