# site-effect-analysis
site-effect analysis is a Python script that predicts the functional effects of site-specific protein phosphorylation. site-effect analysis uses proteomic and phosphoproteomic datasets to create linear regression models that highlight when a site-specific protein phosphorylates or dephosphorylates its target. 

The residual values of these graphs are used as beta values for principal component analysis to graph each site-specific protein as a point on a 2D graph that is also referred to as a cluster graph. This graph colors each point respective to their known site effect using a dataset that stores known site effects. This can be used to analyze clusters or correlation on a specific axis to determine if unknown site-specific proteins are activatory or inhibitory. This process is repeated to make a cluster graph for each gene symbol in the data and all graphs are stored in local directories. 

## Getting Started 
To run this script on your system, there are a few steps you will need to follow to ensure that everything runs smoothly. 

### Dependencies
You will need to have Python 3.14.5 or later installed on your machine. If you need help with installing Python click [here](https://www.python.org/downloads/). If you need to check the current version of Python you have installed, try typing `python --version` into your terminal or `python3 --version` if the first command does not work. 

You will also need to have the following libraries installed by using the following commands:
```
pip install matplotlib
pip install scipy
pip install scikit-learn
```
Installing these libraries should also install numpy, but if numpy is not installed, run the following command: 
```
pip install numpy
```

### Installing 
All you need to do for installation is download siteEffectAnalysis.py into a Python workspace on your system. 

### Data Formatting
Before you execute the program make sure that your datasets are properly formatted to make sure that the Python script runs properly. 

#### Phosphoproteomic and Proteomic datasets
Your phosphoproteomic and proteomic datasets should be formatted as follows:
+ Make sure that the data is saved as a .csv file. Comma delimited is preferable. 
+ The first row in both datasets should consist of labels/headers for each column. Make sure that the label for each target is in the same order across both datasets and that both datasets have the same number of targets.
+ For phosphoproteomic data, the first label should be id, the second geneSymbol, the third variableSites, and the rest are IDs for your targets.
+ For the proteomic data, the first label should be id, the second geneSymbol, and the rest are IDs for your targets. 
+ The first column below the label in both datasets should consist of an ID to label each site-specific protein. 
+ The second column below the label in both datasets must consist of the gene symbol for that row's site-specific protein. (ex. "CDK1") 
+ The third column below the label in the phosphoproteomic dataset must consist of the site(s) for that row's site-specific protein. (ex. "T99t")
+ If there is more than one site for a site-specific protein, divide each site with a space (" "). Make sure there are no extra spaces anywhere. (ex. "S109 S110s")
+ Each subsequent column below the label should consist of the phosphoprotein/total protein value for each respective target.
+ If a phosphoprotein/total protein value is missing, mark the missing value with "NA". This is case-sensitive and must be typed in all caps.
After following these instructions, your phosphoproteomic and proteomic datasets can be used for this Python script.

This is roughly how your phosphoproteomic data should look like: 

<img width="688" height="46" alt="image" src="https://github.com/user-attachments/assets/c9b6ddb9-7c39-4bac-b389-66053c3ec106" />

This is roughly how your proteomic data should look like: 

<img width="563" height="44" alt="image" src="https://github.com/user-attachments/assets/17b2ea48-6a44-4668-ae17-8cc8c3f57d3b" />


#### Known site-effects dataset (site-effect dataset) 
Your site-effect data should be formatted as follows:
+ Make sure that the data is saved as a .csv file. Comma delimited is preferable.
+ The first column must consist of a gene symbol mentioned in the phosphoproteomic and proteomic dataset. (ex. "CDK1")
+ The second column must consist of a site for that gene symbol without the additional letter at the end. (ex. "T31") Do not put multiple sites in the same cell. If a gene symbol has multiple sites, provide it a row for each site.
+ The third column should say what that site does. This could be "phosphorylation" or "dephosphorylation". 
+ The fourth column must consist of the number "1" to imply that the site is an activatory site, "0" to imply that the site effect is unknown, and "-1" to imply that the site is an inhibitory site.
+ Make sure that there are no extra spaces anywhere.
+ After following these instructions, your site-effects dataset can be used for this Python script.

This is roughly what your site-effect data should look like.

<img width="385" height="29" alt="image" src="https://github.com/user-attachments/assets/88d557c6-4395-4f44-94d5-4fd5605ab489" />


### Executing the program
To run siteEffectAnalysis.py, type the following command in the terminal in your Python workspace directory:
```
python3 siteEffectAnalysis.py
```
When the program runs, it will open your file manager and ask you to select your phosphoproteomic data. Select your correctly formatted phosphoproteomic data. Your file manager will then open again and you will have to select your properly formatted proteomic data. Your file manager will open for a third time and you will have to select your properly formatted site-effect data. 

Afterward, the terminal will ask you to provide a name for two directories that will be stored in your current directory. One of these directories will store your cluster graphs and will have the name you provided and _Clusters at the end of the name. The other directory will store your residual graphs and will have the name you provided and _Residual at the end of the name. For example, if you type the name "test" then the two directories that will be created in your current working directory will be test_Clusters and test_Residual. If you make a directory name that already exists, the script will still run, but no new directory will be made. 

After this, if you see the current gene symbol on screen, then the program is running correctly and graphs should shortly start to generate in the clusters and residual directories. 

## Help
This is general advice to troubleshoot any issue you have. 
+ Make sure all dependencies are properly installed. 
+ If you are seeing an index error or any issue relating to the data, refer to Data Formatting and make sure that the data is properly formatted. 
+ If there is an issue with file or directory creation, make sure that your system has sufficient storage and that you have the proper permissions to create directories.
