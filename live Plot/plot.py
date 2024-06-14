import matplotlib.pyplot as plt
import matplotlib.animation
COLOR_LIST = ["blue", "green", "red", "purple", "black"]
def plot(num_of_objectives, csvFile):

    

    fig, ax = plt.subplots() #set up the plot
    x = [] #these two variables store each point data for each objective. 1-D is the objective, 2-D is the data
    y = []
    for i in range(num_of_objectives): #add a sublist of x and y cords for each objective
        x.append([])
        y.append([])

    csvFile.initialRead() #grab the objective names
    
    
    #scatterPlot = ax.scatter(x,y, label="test")
    ax.set(xlabel="Generations", ylabel="DiffObjectives") #set up the axises for our plot


    def animate(_frame): #for animation: https://stackoverflow.com/questions/42722691/python-matplotlib-update-scatter-plot-from-a-function and https://matplotlib.org/stable/users/explain/animations/animations.html
        if(csvFile.checkIfUpdate()==True): #if there's been a change in the file...
            #time.sleep(2)
            newData = csvFile.getNextData() #grab the new data...
            for currentGeneration in newData: #for each generation in the new data...
                for i in range(len(x)):
                    x[i].append(float(currentGeneration[0])) #add new generation data to x and y           converting from string to double: https://stackoverflow.com/questions/482410/how-do-i-convert-a-string-to-a-double-in-python
                for i in range(len(y)):
                    y[i].append(float(currentGeneration[i+1]))

            # print ("x data: ", x)
            # print ("y data: ", y)
            scatterList = [] #a list of the individual plots for each objective
            for i in range(num_of_objectives):
                scatter = ax.scatter(x[i], y[i], label=csvFile.objectives[i], c=COLOR_LIST[i])
                scatterList.append(scatter)

            if csvFile.legendDisplay==False: #display the legend only once
                csvFile.legendDisplay=True
                ax.legend()
            

        

        
        


    #plt.legend()
    ani = matplotlib.animation.FuncAnimation(fig, func=animate, frames=40, interval=500) #animate the plot

    plt.show()