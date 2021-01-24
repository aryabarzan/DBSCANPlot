import math;
import matplotlib.pyplot as plt
#----------------------------------------------------------------
INPUT_FILE_NAME="2d-2c-norm.dat";

marker=['s','o','d','p','*','H','D','h','x'];
color=['c','r','g','b','y'];
markerStyle=[];
for m in marker:
    for c in color:
        markerStyle.append(c+m);


#---------------------- Class Point ----------------------------
class Point:
    def __init__(self,x,y,classNumber):
        self.x=x;
        self.y=y;
        self.classNumber=classNumber;
        self.clusterNumber=0;
        self.reachability_distance=None;
        self.core_distance=None;
        self.processed=False;
    def setX(self,x):
        self.x=x;
    def setY(self,y):
        self.y=y;
    def setClusterNumber(self,clusterNumber):
        self.clusterNumber=clusterNumber;
    def setClassNumber(self,classNumber):
        self.classNumber=classNumber;
    def setReachability_distance(self,reachability_distance):
        self.reachability_distance=reachability_distance;
    def setCore_distance(self,core_distance):
        self.core_distance=core_distance;
    def setCoreDistance(self,dataset,epsilone,minpts):
        neighborsDistance=[];
        for p in dataset:
            d=self.distanceTo(p);
            if(d<=epsilone):
                neighborsDistance.append(d);

        if (len(neighborsDistance)>=minpts):
           neighborsDistance.sort();
           self.core_distance=neighborsDistance[minpts];


    def setProcessed(self,processed):
        self.processed=processed;
    def getX(self):
        return self.x;
    def getY(self):
        return self.y;
    def getClusterNumber(self):
        return self.clusterNumber;
    def getClassNumber(self):
        return self.classNumber;

    def getReachability_distance(self):
        return self.reachability_distance;
    def getCore_distance(self):
        return self.core_distance;
    def getProcessed(self):
        return self.processed;
    def distanceTo(self,to):
        return math.sqrt(math.pow(self.getX() -to.getX(), 2)+math.pow(self.getY() -to.getY(), 2));
    def getNeighboresEpsilone(self,epsilone,dataset):
    #Direct Density Reachible
        neighbors=[];
        for p in dataset:
            if(self.distanceTo(p)<=epsilone):
                neighbors.append(p);

        return neighbors;
    
    def clone(self):
        p=Point(self.x,self.y,self.classNumber)
        return p;

    def printPoint(self):
        print '(',self.getX(),',',self.getY(),')->',self.getReachability_distance(),'-->',self.getProcessed(),'cd',self.getCore_distance(),'ccclass',self.getClassNumber(),'cluster number',self.getClusterNumber();
#---------------------- End Class Point ---------------------------
def clone(d):
    dataSet=[]
    for point in d:
        dataSet.append(point.clone())
    return dataSet;

#====================== Reading data from input file ===========================
def readDataSet(inputFileName):
    fin=open(inputFileName,"r");
    dataSet=[];

    for line in fin: # read rest of lines
        v=line.split();
        x=float(v[0]);
        y=float(v[1]);
        classNumber=int(v[2])
        p=Point(x,y,classNumber);
        dataSet.append(p);


    fin.close();
    return dataSet;

#====================== Functions ========================================
def DBSCAN(dataSet,minPts,epsilone,figureNumber):
    dataSet=clone(dataSet)
    d_Unprocessed=list(dataSet);
    no_of_clusters=0;
    clusters=[];
    while(len(d_Unprocessed)>0):
        p=d_Unprocessed[0];
        neighborsOfP=findNeighboreEpsilone(p,epsilone,dataSet);
        if(len(neighborsOfP)<minPts):#p is a noise point
            d_Unprocessed.remove(p);
        else :#p is a core point
            no_of_clusters=no_of_clusters+1;
            p.setClusterNumber(no_of_clusters);
            markCluster(p,epsilone,minPts,dataSet);
            #d_Unprocessed=d_Unprocessed-new cluster


            for p1 in list(d_Unprocessed):
                if(p1.getClusterNumber()>0):
                    d_Unprocessed.remove(p1);
    plotClusters(dataSet, no_of_clusters, figureNumber)

def findNeighboreEpsilone(point,epsilone,dataSet):
    #Direct Density Reachible
    neighbors=[];
    for p in dataSet:
        if(point.distanceTo(p)<=epsilone):
            neighbors.append(p);


    return neighbors;

def markCluster(point,epsilone,minPts,dataSet):

    #Direct Density Reachible
    neighbors=findNeighboreEpsilone(point,epsilone,dataSet);
    if(len(neighbors)>=minPts):#//p is a core point
        for n in neighbors:
            if(n.getClusterNumber()<=0):
                n.setClusterNumber(point.getClusterNumber());
                markCluster(n,epsilone,minPts,dataSet);

def plotClusters(dataSet,no_of_clusters,figureNumber):
    plt.figure(figureNumber)
    clusters=[[]  for v in range(no_of_clusters+1)];

    min_x=dataSet[0].getX();
    max_x=dataSet[0].getX();
    min_y=dataSet[0].getY();
    max_y=dataSet[0].getY();
  
    for p in dataSet:
        clusters[p.getClusterNumber()].append(p);
        if(p.getX()<min_x):
            min_x=p.getX();
        if(p.getX()>max_x):
            max_x=p.getX();
        if(p.getY()<min_y):
            min_y=p.getY();
        if(p.getY()>max_y):
            max_y=p.getY();

    #-------------------- plot clusters -------------------------------
    for c in clusters:
        x=[];
        y=[];
        clusterNumber=0;
       
        for p in c:
            clusterNumber=p.getClusterNumber();
            x.append(p.getX());
            y.append(p.getY());
           # p.printPoint()
        if(clusterNumber==0):
            plt.plot(x,y,'k+');
        else:
            plt.plot(x,y,markerStyle[clusterNumber%len(markerStyle)]);
       # print clusterNumber+1

        
    plt.xlim(min_x-1,max_x+1);
    plt.ylim(min_y-1,max_y+1);

#====================== End of Functions =================================


minPts=4


inputFileNames=['2d-2c-norm','2d-4c-no9','2d-4c-norm']
threshold=0
figureNumber=0;
epsilones=[0.2, 0.4, 0.6,0.8,1,1.2,1.4,1.6]
for inputFileName in inputFileNames:
    for epsilone in epsilones:
        dataSet=readDataSet(inputFileName+".dat");
        figureNumber=figureNumber+1
        DBSCAN(dataSet, minPts, epsilone,figureNumber)

plt.show()