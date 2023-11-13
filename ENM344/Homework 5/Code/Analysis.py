import csv
import numpy as np
import time
import matplotlib.pyplot as plt
from scipy.stats import t as ttt

# If varying specides, change this name
species = 'Xiphius gladius'
# If downloading to run, change below paths for to access data properly
species_path = "Homework 5/BioTIMEQuery_24_06_2021.csv"
climate_path = "Homework 5/Annual_Surface_Temperature_Change.csv"
# Analyze raw animal data to get yearly swordfish data.
speciesDict = {}
SigPueCount = {}
with open(species_path, 'r', encoding='utf-8') as csv_file:
    csv_reader = csv.reader(csv_file)
    for i, row in enumerate(csv_reader):
        if i == 0:
            for j in range(len(row)):
                print(str(j)+":"+row[j])
        else:
            if row[14] in speciesDict.keys():
                speciesDict[row[14]] = speciesDict[row[14]] + 1
            else:
                speciesDict[row[14]] = 1
            
            if(row[14] == species):
                if int(row[4]) in SigPueCount.keys():
                    SigPueCount[int(row[4])] = SigPueCount[int(row[4])] + 1
                else:
                    SigPueCount[int(row[4])] = 1

        if(i % 1000000 == 0):
            print(i)

# Get Average Temperature data by year
temperature = {}
yearCount = {}
# Make proper entries for dictionaries
for i in range(14):
    year = 1992 + i
    temperature[year] = float(0)
    yearCount[year] = 0
# Loop through and extract temperature data
with open(climate_path, 'r', encoding='utf-8') as csv_file:
    csv_reader = csv.reader(csv_file)
    for i, row in enumerate(csv_reader):
        if i == 0:
            for j in range(len(row)):
                print(str(j)+":"+row[j])
        else:
            goodRow = True
            for j in range(14):
                if(row[41+j] == ''):
                    goodRow = False
            if(goodRow):
                for j in range(14):
                    year = 1992+j
                    temperature[year] = temperature[year] + float(row[41+j])
                    yearCount[year] = yearCount[year] + 1
# Distill data into average temperature by year
averageTemps = {}
for i in range(14):
    year = 1992 + i
    averageTemps[year] = temperature[year]/yearCount[year]

# Combine temperature and population data by year
popTempData = {}
for i in range(14):
    year = 1992 + i 
    popTempData[averageTemps[year]] = SigPueCount[year]

popTempKeys = list(popTempData.keys())
popTempVals = list(popTempData.values())

# Now we have combined data, we need to find significance
# Pearson Correlation Coefficient
xmean = sum(popTempKeys)/len(popTempKeys)
ymean = sum(popTempVals)/len(popTempVals)

numer = 0
denom1 = 0
denom2 = 0
for i in range(len(popTempKeys)):
    numer += (popTempKeys[i]-xmean)*(popTempVals[i]-ymean)
    denom1 += (popTempKeys[i]-xmean)*(popTempKeys[i]-xmean)
    denom2 += (popTempVals[i]-ymean)*(popTempVals[i]-ymean)
r = numer/(np.sqrt(denom1)*np.sqrt(denom2))
# print(r)

# Now we need to calculate the t-value associated with this.
N = len(popTempKeys)
t = abs(r*np.sqrt(N-2)/np.sqrt(1-r**2)) # Took absolute value because we expect it to be negative, but only mag matters

# Find associated alpha that gives a critical t value below this.
a = 0.05
tcrit = ttt.ppf(1-a, N-2)
error = t-tcrit
iter = 0
while(abs(error) > 10e-6):
    a -= error/(100000*max(tcrit,t))
    da = error/(100000*max(tcrit,t))
    tcrit = ttt.ppf(1-a, N-2)
    error = t-tcrit
    if(iter % 2500 == 0):
        print(error)
    iter = iter+1
time.sleep(1)
print()
print("Final Results")
print("t = "+str(t))
print("tcrit = "+str(tcrit))
print("a = "+str(a))

# Now that we have our a, and it is below 0.05, we are done with NHST and have proven our idea.
# We can also try bootstrapping
Niter = 100000
bootstrap_distributions = []
r_original = r

popTempKeys = np.array(list(popTempData.keys()))
popTempVals = np.array(list(popTempData.values()))
for i in range(Niter):
    sample_indices = np.random.choice(range(len(popTempKeys)), size=len(popTempVals), replace=True)
    bootstrap_sample_x = [popTempKeys[k] for k in sample_indices]
    bootstrap_sample_y = [popTempVals[k] for k in sample_indices]
    
    xmean = sum(bootstrap_sample_x)/len(bootstrap_sample_x)
    ymean = sum(bootstrap_sample_y)/len(bootstrap_sample_y)
    
    numer = 0
    denom1 = 0
    denom2 = 0
    for j in range(len(bootstrap_sample_x)):
        numer += (bootstrap_sample_x[j]-xmean)*(bootstrap_sample_y[j]-ymean)
        denom1 += (bootstrap_sample_x[j]-xmean)**2
        denom2 += (bootstrap_sample_y[j]-ymean)**2
    bootstrap_r = numer/(np.sqrt(denom1)*np.sqrt(denom2))
    bootstrap_distributions.append(bootstrap_r)
    if(i % 10000 == 0):
        print(i)
    

extreme_values_count = sum(abs(rs) >= abs(r) for rs in bootstrap_distributions)
p_value = extreme_values_count / Niter

print(p_value)

plt.figure()
plt.scatter(SigPueCount.keys(),SigPueCount.values())
plt.xlabel("Year")
plt.ylabel("Population")
plt.title("Xiphius gladius Population")

plt.figure()
plt.scatter(averageTemps.keys(),averageTemps.values())
plt.xlabel("Year")
plt.ylabel("Temperature Increase [C]")
plt.title("Global Average Temperature Trends")

plt.figure()
plt.scatter(popTempData.keys(),popTempData.values())
plt.xlabel("Temperature Increase [C]")
plt.ylabel("Population")
plt.title("Xiphius gladius Population Trend vs Temperature Change")
plt.show()
