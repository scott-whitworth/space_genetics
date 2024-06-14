import os
import time
import csv #found from: https://www.geeksforgeeks.org/reading-csv-files-in-python/
import copy #learned about deep copying here: https://stackoverflow.com/questions/17873384/how-to-deep-copy-a-list

class Parser:
    def __init__(self, fileName, objectiveNum):
        self.row = 0 #what row to look at for new data. (That way we don't start at the beginning of the csv file each time we update it)
        self.storedChange = 0 #This variable detects if the csv file has been updated. found solution for file change from: https://stackoverflow.com/questions/182197/how-do-i-watch-a-file-for-changes
        self.fileName = fileName
        self.objectiveNum = objectiveNum
        self.objectives = [] #lists the names of each objective so that we can label the data in plot.py
        self.legendDisplay = False #displays if we need to add a legend to the plot (since we should only do it once). Put it here so that animate could access it.
    #this function returns true if the file has been updated
    def checkIfUpdate(self):
        currentChange = os.stat(self.fileName).st_mtime
        if currentChange!=self.storedChange:
            self.storedChange = currentChange #our new stored change
            return True #there is an update!
        return False
    
    #goes through the first line of the csv to grab the objective names
    def initialRead(self):
        with open(self.fileName, mode = 'r') as file:   #for correct way to open files: https://stackoverflow.com/questions/11555468/how-should-i-read-a-file-line-by-line-in-python
            csvFile = csv.reader(file)
            for lines in csvFile:     
                for i in range(self.objectiveNum):
                    self.objectives.append(lines[i+1]) #the first one is just generation
                self.row=1
                
                return
    #for each line (index of row and beyond) grab the future generation data.
    def getNextData(self):
        with open(self.fileName, mode = 'r') as file:
            csvFile = csv.reader(file)
            row = 0 #what row we're checking now
            dataReturn = [] 
            for lines in csvFile: #for each line in the csv file...
                currentGenData = [] #data for this generation
                if (row>=self.row): #only start grabbing data we haven't grabbed before
                    for i in range(self.objectiveNum+1): #for each objective (skipping the first column since it's just generation number)...
                        currentGenData.append(lines[i]) #add this data to the objective data for this generation
                    dataReturn.append(copy.deepcopy(currentGenData)) #add this to total new data collected
                
                row+=1 

            self.row=row #since we've looked down to the row-th column, set our self.row to this so that we don't look at previously processed data
            return dataReturn


            

    

