import csv
import numpy as np
import matplotlib.pyplot as plt

# Helper Functions
def timeStringToSeconds(time):
    stringSec = time[-2:]
    stringMin = time[-5:-3]
    stringHour = time[0:1]
    seconds = int(stringSec) + int(stringMin)*60 + int(stringHour)*3600
    return seconds

# Data Taken From https://web.chessdigits.com/data
# Move number vs Time Spent
timeSpent = []
moveNumberCount = []
averageTime = []
# Counting Variables
categoryCount = {"Bullet": 0,"Blitz":0, "Rapid":0, "Classical":0}
dRatingWins = {}
dRatingCount = {}
rateCount = {}
firstMoveCount = {"d4 Win": 0, "d4 Loss": 0,"d4 Draw":0, "e4 Win":0, "e4 Loss":0,"e4 Draw":0, "Other Win":0, "Other Loss": 0,"Other Draw":0}
# Start Reading Data
csvfile = open('200k_blitz_rapid_classical_bullet.csv')
csvreader = csv.reader(csvfile)
maxCentipawns = 327.65
MateIn = 0

rowCount = 0
for row in csvreader:
    rowCount = rowCount + 1
    print(rowCount)
    if(rowCount == 1): # Deal with Header Row
        print("")
    else:
        # This is just an early break for testing, comment out for full run
        if(rowCount > 100000000):
            break
        
        # First, let us count how many of each category is played
        if(row[621] in categoryCount):
            categoryCount[row[621]] += 1
        else:
            categoryCount[row[621]] = 1
        
        # Now, let us see what the starting time and the increment were for this game
        plusIndex = row[13].find("+")
        timeMax = int(row[13][:plusIndex])
        increment = int(row[13][plusIndex+1:])
        
        # First we get the dt vs move number data
        timePrevWhite = timeStringToSeconds(row[421])
        timePrevBlack = timeStringToSeconds(row[421])
        moveNumber = 0
        for i in range(0,200):
            # Move Number vs Time Spent
            moveNumber = i + 1
            if(len(row[421+i]) != 0):
                if(len(moveNumberCount) < moveNumber):
                    moveNumberCount.append(float(0))
                    timeSpent.append(float(0))
                if i % 2 == 0: # If White Move
                    moveNumberCount[i] += 1
                    timeSpent[i] += timePrevWhite - timeStringToSeconds(row[421+i])+increment
                    timePrevWhite = timeStringToSeconds(row[421+i])
                else:
                    moveNumberCount[i] += 1
                    timeSpent[i] += timePrevBlack - timeStringToSeconds(row[421+i])+increment
                    timePrevBlack = timeStringToSeconds(row[421+i])
        
        # Now we calculate the rating difference win rate
        window = float(20)
        dR = int(row[17]) - int(row[3]) # Positive for White Higher
        dR = round(dR/window)*window
        if(dR in dRatingCount):
            dRatingCount[dR] += 1
        else:
            dRatingCount[dR] = float(1)
        if(row[9][-1] != "*"):
            if(int(row[9][-1]) == 0):
                if(dR in dRatingWins):
                    dRatingWins[dR] += 1
                else:
                    print("Add")
                    dRatingWins[dR] = float(1)
        # Rating Histogram
        rating = int(row[17])
        if(rating in rateCount):
            rateCount[rating] += 1
        else:
            rateCount[rating] = 1
        rating = int(row[3])
        if(rating in rateCount):
            rateCount[rating] += 1
        else:
            rateCount[rating] = 1
        
        #d4 vs e4 vs other win rates
        
        condition = ''
        move = ''
        if(len(row[9]) == 3):
            if(int(row[9][-1]) == 0):
                condition = ' Win'
            else:
                condition = ' Loss'
        else:
            condition = ' Draw'
        if(row[21] == "d4" or row[21] == "e4"):
            move = row[21]
        else:
            move = "Other"
        phrase = move+condition
        firstMoveCount[phrase] += 1
        
# Post Process Time v Move Number
for i in range(0,len(moveNumberCount)-1):
    if(i%2 == 0):
        timeAverage = (timeSpent[i]+timeSpent[i+1])/(moveNumberCount[i]+moveNumberCount[i+1])
        averageTime.append(timeAverage)
# Post Processing dRating Win Rate
myKeys = list(dRatingWins.keys())
myKeys.sort()
dR_sorted = {i: dRatingWins[i] for i in myKeys}
dRCountSorted = {i: 100*dRatingCount[i]/2*rowCount for i in myKeys}
myValues = list(dR_sorted.values())
myOccurances = list(dRCountSorted.values())
for i in range(0,len(myKeys)):
    myValues[i] = 100 * myValues[i] / dRatingCount[myKeys[i]]
    
myRatings = list(rateCount.keys())
myRatings.sort()
sortedRating = {i: 100*rateCount[i]/rowCount for i in myRatings}
myRates = list(sortedRating.values())
# Display and Testing
print("===========TESTING=============")
# Move Number vs Time
plt.figure(1)
plt.plot(averageTime,'r-')
plt.xlabel("Move Number")
plt.ylabel("Average Time Taken (s)")
plt.title("Average Time vs Move Number Trends")
plt.grid(True)
# dRating vs Win Rate
plt.figure(2)
ax1 = plt.subplot()
l1 = ax1.plot(myKeys,myValues,'r-')
ax2 = ax1.twinx()
l2 = ax2.plot(myKeys,myOccurances,'b-')
plt.xlabel("Rating Difference (elo)")
ax1.set_ylabel("Win Rate (%)", color = 'r')
ax2.set_ylabel("Difference Frequency (%)", color = 'b')
ax2.set_ylim(ax1.get_ylim())
plt.title("Rating Difference Statistics")
plt.grid(True)
plt.xlim(-800,800)
plt.ylim(bottom = 0)
# Rating Histogram
plt.figure(3)
plt.plot(myRatings,myRates,'r-')
plt.xlabel("Rating")
plt.ylabel("Percentage of Players")
plt.title("Rating Distribution")
plt.grid(True)
# Winning First Move Pie Chart
plt.figure(4)
myMoves = list(firstMoveCount.keys())
myMovesCounts = list(firstMoveCount.values())
colorSet = ['#ADD8E6', '#6495ED', '#0000FF', '#98FB98', '#32CD32', '#006400', '#FFB6C1', '#FF69B4', '#FF1493']
plt.pie(myMovesCounts,explode = None,labels = myMoves, colors = colorSet)
plt.title("Opening Move Success")
# Categories Bar Graph
plt.figure(5)
myCats = list(categoryCount.keys())
categoryCount = {i: 100*categoryCount[i]/rowCount for i in myCats}
myNums = list(categoryCount.values())
plt.bar(myCats,myNums,color = 'blue')
plt.xlabel("Time Control")
plt.ylabel("Percent of all Games")
plt.title("Time Control Frequency")
plt.show()