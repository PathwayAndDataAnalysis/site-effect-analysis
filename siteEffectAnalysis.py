# Imported libraries.
import csv
import math
import matplotlib.patches as patches
import matplotlib.pyplot as pyplot
import numpy
import os
from scipy import stats
from sklearn.decomposition import PCA
from tkinter import filedialog

# File reading. 
phosphoproteome = filedialog.askopenfile("r" , title="Select your phosphoproteomic data.")
proteome = filedialog.askopenfile("r", title="Select your proteomic data.")
siteEffects = filedialog.askopenfile("r", title="Select your site effect data. ")

phosphoproteomeReader = csv.reader(phosphoproteome)
proteomeReader = csv.reader(proteome)
siteEffectsReader = csv.reader(siteEffects)

phosphoproteomeRows = []
proteomeRows = []
siteEffectsRows = []

for row in phosphoproteomeReader:
    phosphoproteomeRows.append(row.copy())
for row in proteomeReader:
    proteomeRows.append(row.copy())
for row in siteEffectsReader:
    siteEffectsRows.append(row.copy())

# Functions.

# Features scales the x data provided using normalization. Makes all values floats instead of strings and replaces all NA values with 0.0.
def featureScale(xData):
    scaledData = xData.copy()
    mean = 0.0
    std = 0.0
    n = 0
    for i in range (0, len(scaledData)):
        if scaledData[i] != "NA":
            scaledData[i] = float(scaledData[i])
            mean += scaledData[i]
            n += 1
    mean /= n
    for i in range(0, len(scaledData)):
        if scaledData[i] != "NA":
            temp = scaledData[i] - mean
            temp = math.fabs(temp)
            temp = temp * temp
            std += temp
    std = math.sqrt(std)
    for i in range(0, len(scaledData)):
        if scaledData[i] != "NA":
            scaledData[i] -= mean
            scaledData[i] /= std
        else:
            scaledData[i] = 0.0
    return scaledData

# Makes all values in the provided row of data floats instead of strings and repalces all NA values with 0.0.
def makeFloats(data):
    newData = data.copy()
    for i in range(0, len(newData)):
        if newData[i] != "NA":
            newData[i] = float(newData[i])
        else: 
            newData[i] = 0
    return newData

# Finds and returns the number of cells that contain data from start to the end of the row.
def filledCells(proteinRow, start):
    count = 0
    for i in range(start, len(proteinRow)):
        if proteinRow[i] != "NA":
            count += 1
    return count

# Finds and returns each protein with a matching gene symbol in an array of rows that represent proteins.
def findMatchingGS(geneSymbol, proteinData):
    matching = []
    for row in proteinData:
        if row[1] == geneSymbol:
            matching.append(row.copy())
    return matching

# Finds and returns each site effect with a matching gene symbol in an array of row that represents site effects.
def findSites(geneSymbol, siteEffectData):
    matching = []
    for row in siteEffectData:
        if row[0] == geneSymbol:
            matching.append(row.copy())
    return matching

# Locates the first row with a matching site in an array rows that represent site effects.
def findSite(site, siteEffectData):
    for row in siteEffectData:
        if row[1] == site:
            return row.copy()
    return []

# A function to predict values for linear regression models. 
slope = 0
intercept = 0
def f(x):
    return slope * x + intercept

phosphoproteomeStart = 3
proteomeStart = 2

# Patches used for the legend of cluster graphs. 
red = patches.Patch(color="red", label="Inhibatory site")
green = patches.Patch(color="green", label="Activatory site")
yellow = patches.Patch(color="yellow", label="Mixed")
gray = patches.Patch(color="gray", label="Unknown site")
handlesL = [red, green, yellow, gray]

# Accesses the path of the current working directory. 
path = os.getcwd()
# Asks the user to provide a name that will be used to make a directory that stores residual graphs 
# and another directory that stores cluster graphs. Existing directories can also be used. 
print("Provide the name of your data that will be used to name the directories where graphs are stored: ", end="")
name = input()
nameResidual = name + "_Residual"
try:
    os.mkdir(nameResidual)
    print(path + "\\" + nameResidual + " has been sucessfully created.")
except FileExistsError:
    print(path + "\\" + nameResidual + " already exists.")
nameCluster = name + "_Clusters"
try:
    os.mkdir(nameCluster)
    print(path + "\\" + nameCluster + " has been succesfully created.")
except FileExistsError:
    print(path + "\\" + nameCluster + " already exists.")

allBetaVals = [] # Keeps track of the beta values for each protein of a gene symbol.
colors = [] # Keeps track of the color of each point in the cluster graph based on its known functional effect.
allSites = [] # Keeps track of each site to label each point in the cluster graphs. 
phoRN = 1 # Keeps track of the current row number in the phosphoproteomic data.
geneSymbol = "" # Keeps track of the current gene symbol. 
while phoRN < len(phosphoproteomeRows):
    betaVals = [] # Keeps track of the beta values for the current protein. 
    # Checks if the current gene sybmol (if there is one) is different from the next to create a cluster graph for the current gene symbol. 
    if (phoRN != 1 and geneSymbol != phosphoproteomeRows[phoRN][1]):
        allBetaVals = numpy.array(allBetaVals)
        # If there are less than two proteins represented in allBetVals, PCA cannot be done for this gene symbol and no cluster graph is made.
        if len(allBetaVals) >= 2:
            pca = PCA(n_components=2)
            pca.fit(allBetaVals)
            pcaB = pca.transform(allBetaVals)
            pca1 = pcaB[:, 0] # First principal component. 
            pca2 = pcaB[:, 1] # Second principal component. 
            # Plots each point with the color from colors and the label from allSites. 
            for i in range(0, len(pca1)):
                pyplot.plot(pca1[i], pca2[i], 'o', color=colors[i])
                pyplot.annotate(allSites[i], (pca1[i], pca2[i]), alpha=0.6, color="gray")
            # Labels the graph's legend with the patches from handlesL. 
            pyplot.legend(handles=handlesL)
            # Titles the graph with the name provided by the user earlier and the current gene symbol in parentheses. 
            pyplot.title(name + " (" + geneSymbol + ")")
            # Labels the axes with the principal components. 
            pyplot.xlabel("First Principal Component")
            pyplot.ylabel("Second Principal Component")
            # Saves the cluster graph to the cluster directory created earlier. 
            pyplot.savefig(path + "\\" + nameCluster + "\\" + geneSymbol + ".png")
            # Closes the graph to not interfere with the next graph. 
            pyplot.close()
        # Makes allBetaVals, colors, and allSites empty since they now have to represent the next geneSymbol. 
        allBetaVals = allBetaVals.tolist()
        allBetaVals = []
        colors = []
        allSites = []
    geneSymbol = phosphoproteomeRows[phoRN][1] # Updates the current geneSymbol. 
    sites = phosphoproteomeRows[phoRN][2] # Keeps track of the site(s) of the current protein. 
    gsProteomeRows = findMatchingGS(geneSymbol, proteomeRows) # Accesses all of the protemic information regarding the current gene symbol.
    # If there is no proteomic information, the current protein is skipped. 
    if gsProteomeRows == []:
        phoRN += 1
        continue
    # Selects the proteomic information with the most information to generate the most accurate residual graph. 
    mostCells = 0
    maxIndex = 0
    for i in range(0, len(gsProteomeRows)):
        currentCells = filledCells(gsProteomeRows[i], 18)
        if currentCells > mostCells:
            maxIndex = i
            mostCells = currentCells
    x = gsProteomeRows[maxIndex][proteomeStart:].copy() # The best row of proteomic information is used as the x data for the residual graph. 
    y = phosphoproteomeRows[phoRN][phosphoproteomeStart:].copy() # The current phosphoproteomic row is used as the y data for the residual graph. 
    xCells = float(filledCells(x, 0)) / len(x) # Keeps track of the percentage of filled cells in the x data. 
    yCells = float(filledCells(y, 0)) / len(y)  # Keeps track of the percentage of the filled celts in the y data. 
    # If less than 60% of the y or x data is filled, then the current protein is skipped as the data in the graph would be too inaccurate. 
    if xCells < 0.60 or yCells < 0.60:
        phoRN += 1
        continue
    x = featureScale(x)
    y = makeFloats(y)
    # Creates a linear regression model for the x and y data that will be used to calculate residuals. 
    slope, intercept, r, p, std_err = stats.linregress(x, y)
    corr, p_val = stats.pearsonr(x, y)
    # Displays the current gene symbol, correlation coefficient, and p-value.  
    print("Current gene symbol: " + geneSymbol)
    print("Correlation: ", end="")
    print(corr)
    print("P-value: ", end="")
    print(p_val)
    # Plots the curernt linear regression model and creates the best-fit line. 
    model = list(map(f, x))
    pyplot.scatter(x, y)
    # Titles the graph with the current gene symbol and site(s) of the current protein. 
    pyplot.title(geneSymbol + " (" + sites + ")")
    # Labels the axes with the data's source. 
    pyplot.xlabel("Total Protein")
    pyplot.ylabel("Phosphoprotein Value")
    # Plots the best fit line. 
    pyplot.plot(x, model)
    # Saves the residual graph to the residual directory created earlier. 
    pyplot.savefig(path + "\\" + nameResidual + "\\" + geneSymbol + "_" + sites + ".png")
    # Closes the graph to not interfere with the next graph. 
    pyplot.close()
    # Calculates every residual value of the current linear regression model to be the beta values of the curernt protein. 
    for i in range(0, len(x)):
        val = y[i] - f(x[i])
        betaVals.append(val)
    # Saves the beta values of the current protein to allBetaVals. 
    allBetaVals.append(betaVals)
    # Finds the known site effects of the current protein and splits each site to be an indiviudal string. 
    gsSiteEffects = findSites(geneSymbol, siteEffectsRows)
    sites = sites.split(" ")
    incOrDec = "0" # Keeps track if the current site specific protein is activatory ("1"), unknown ("0"), inhibitory ("-1"), or mixed ("2").
    for site in sites:
        siteEff = findSite(site[0:len(site) - 1], gsSiteEffects)
        # If there is no known site effect, the protein is labeled as unknown. 
        if siteEff != []:
            # If the current protein is has an activatory and inhibitory site, then it is labeled as mixed. 
            if incOrDec != siteEff[3] and incOrDec != "0":
                incOrDec = "2"
                break
            else: 
                incOrDec = siteEff[3]
    # Appends the color correlated to the current protein's functional effect to colors to be used later. 
    if incOrDec == "-1":
        colors.append("red")
    if incOrDec == "0":
        colors.append("gray")
    if incOrDec == "1":
        colors.append("green")
    if incOrDec == "2":
        colors.append("yellow")
    # Appends the first site of the current protein to allSites to be used later. 
    allSites.append(sites[0][0:len(sites[0]) - 1])
    # Iterates to the next protein in the phosphoproteomic data. 
    phoRN += 1

# Checks if the last gene symbol can create a cluster graph. 
if len(allBetaVals) > 2:
    pca = PCA(n_components=2)
    pca.fit(allBetaVals)
    pcaB = pca.transform(allBetaVals)
    pca1 = pcaB[:, 0] # First principal component.
    pca2 = pcaB[:, 0] # Second principal component 
    # Plots each point with the color from colors and the label from allSites. 
    for i in range(0, len(pca1)):
        pyplot.plot(pca[i], pca[i], 'o', color=colors[i])
        pyplot.annotate(allSites[i], (pca1[i], pca2[i]), alpha=0.6, color="gray") 
    # Labels the graph's legend with the patches from handlesL. 
    pyplot.legend(handles=handlesL)
    # Titles the graph with the name provided by the user earlier and the current symbol in parentheses. 
    pyplot.title(name + " (" + geneSymbol + ")")
    # Labels the axes with the principal components. 
    pyplot.xlabel("First principal component")
    pyplot.ylabel("Second principal component")
    # Saves the cluster graph to the cluster directory created earlier. 
    pyplot.savefig(path + "\\" + nameCluster + "\\" + geneSymbol + ".png")
    # Closes the graph. 
    pyplot.close()

# File Closure. 
phosphoproteome.close()
proteome.close()
siteEffects.close()