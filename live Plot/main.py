import parse
import numpy
import plot
import argparse

#first, parse the arguments: genPerformance file location and the number of objectives

#------Argument Parsing ------------------------------------------------------------------------ Help from https://docs.python.org/dev/library/argparse.html
argParser = argparse.ArgumentParser(prog="live plot", description="Tracks generation performance in real time")
argParser.add_argument("filename")
argParser.add_argument('-o', '--objectivesnum', default=1, help="List the number of objectives the simulation is currently optimizing for")
args = argParser.parse_args()

NUM_OF_OBJ = int(args.objectivesnum)
#---------------------------------------------------------------------------



def main():
    #create a newFile object to mostly parse through the csv file
    newFile = parse.Parser(args.filename, NUM_OF_OBJ)
    
    #call plot to update and plot the data
    plot.plot(NUM_OF_OBJ, newFile)
 
if __name__ == "__main__":
    main()