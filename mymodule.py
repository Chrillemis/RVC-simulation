import numpy as np

#Tried to use in ordered cleaning with obstacles
def SeekCoord(Room, RVCX, RVCY, t, N, j, i):
    seekX = N
    seekY = N
    if Room[N-int(np.floor(RVCY[t]))-1][int(np.floor(RVCX[t]))] != 0:
        if j%2 == 0:
            for k in range(N-i):
                if Room[N-int(np.floor(RVCY[t]))-1][int(np.floor(RVCX[t]))+k] == 0:
                    seekX = i+k
                    seekY = j
                    break
        if j%2 == 1:
            for k in range(i):
                if Room[N-int(np.floor(RVCY[t]))-1][int(np.floor(RVCX[t]))-k] == 0:
                    seekX = i-k
                    seekY = j
                    break
    return seekX, seekY


class Robot():
    def __init__(self, Grid, Tmax, eff, N, d, drmax, Room):
        self.Grid = Grid
        self.Tmax = Tmax
        self.eff = eff
        self.N = N
        self.d = d
        self.drmax = drmax
        self.L = self.N-0.001
        self.RVCX = [self.L/2]
        self.RVCY = [self.L/2]
        self.Frames = []
        if isinstance(Room, np.ndarray or list):
            self.Room = Room
        if isinstance(Room, int):
            self.Room = np.zeros((self.N,self.N))

    
class SnakeBot(Robot)  :  
    def SnakeRVC(self, t=0, i=0, j=0):
        while self.Grid.sum() >= self.eff:
            if t > self.Tmax:
                break
            t += 1     
            self.Grid[j][i] += -self.d*self.eff
            if j%2 == 0:
                i += 1
                if i == self.N:
                    j += 1
            if j%2 == 1:
                i += -1 
                if i == -1:
                    i = 0
                    j += 1
            if j == self.N:
                j = 0 
        for i in range(self.N):
            for j in range(self.N):
                if self.Grid[i][j] < self.eff:
                    self.Grid[i][j] = 0
        return self.Grid    
    def SnakeRVCObstacle(self, t=0, i=0, j=0, seekX=0, seekY=0): #I believe i might have overcomplicated this
        self.RVCX  = [0]
        self.RVCY = [0]
        def clean(j, i):
            # if j > self.N-1:
            #     j = self.N-1
            # if i > self.N-1:
            #     i = self.N-1
            self.Grid[j][i] += -self.d*self.eff
            
        def seek(j, i, seekX, seekY):
            if j%2 == 0:
                while i != 0:
                    i += -1
                    self.RVCX.append(i)
                    self.RVCY.append(j)                        
                    clean(j = j, i = i)
                while j != 0:
                    j += -1
                    self.RVCX.append(i)
                    self.RVCY.append(j)                        
                    clean(j = j, i = i)
                while i != seekX:
                    i += 1
                    self.RVCX.append(i)
                    self.RVCY.append(j)                        
                    clean(j = j, i = i)
                while j != seekY:
                    j += 1
                    self.RVCX.append(i)
                    self.RVCY.append(j)                        
                    clean(j = j, i = i)   
            if j%2 == 1:
                while i != self.N-1:
                    i += 1
                    self.RVCX.append(i)
                    self.RVCY.append(j)                        
                    clean(j = j, i = i)
                while j != 0:
                    j += -1
                    self.RVCX.append(i)
                    self.RVCY.append(j)                        
                    clean(j = j, i = i)
                while i != seekX:
                    i += -1
                    self.RVCX.append(i)
                    self.RVCY.append(j)                        
                    clean(j = j, i = i)
                while j != seekY:
                    j += 1
                    self.RVCX.append(i)
                    self.RVCY.append(j)                        
                    clean(j = j, i = i)
                
        
        while self.Grid.sum() >= self.eff:
            if t > self.Tmax:
                break
            self.RVCX.append(i)
            self.RVCY.append(j)   
            t += 1     
            if self.Room[j][i] != 0:
                seekX, seekY = SeekCoord(self.Room, self.RVCX, self.RVCY, t, self.N, j, i)
                j = self.RVCY[t]
                i = self.RVCX[t]
                print(j)
                print(i)
                seek(j = j, i = i, seekX = seekX, seekY = seekY)
                i = seekX
                j = seekY

                    
            clean(j = j, i = i)
            if j%2 == 0:
                i += 1
                if i == self.N:
                    j += 1
            if j%2 == 1:
                i += -1 
                if i == -1:
                    i = 0
                    j += 1
            if j == self.N:
                j = 0 
        for i in range(self.N):
            for j in range(self.N):
                if self.Grid[i][j] < self.eff:
                    self.Grid[i][j] = 0
        return self.Grid, self.RVCX, self.RVCY, seekX, seekY
    
class RandoBot(Robot):
    def RandomPath(self, t=0):
        
        UpBoundX = self.L
        UpBoundY = self.L
        LowBoundX = 0.001
        LowBoundY = 0.001

        while self.Tmax >= t:
            t+=1
            
            theta = np.random.uniform(0,2*np.pi)
            
            self.RVCX.append(np.cos(theta)*self.drmax+self.RVCX[t-1])
            self.RVCY.append(np.sin(theta)*self.drmax+self.RVCY[t-1])
            
            #If hits a wall, bounce off
            if self.RVCX[t] < LowBoundX:
                self.RVCX[t] = LowBoundX
            if self.RVCX[t] > UpBoundX:
                self.RVCX[t] = UpBoundX
            if self.RVCY[t] < LowBoundY:    
                self.RVCY[t] = LowBoundY
            if self.RVCY[t] > UpBoundY:
                self.RVCY[t] = UpBoundY
            
            #If hits an obstacle, bounce off
            if self.Room[self.N-int(np.floor(self.RVCY[t]))-1][int(np.floor(self.RVCX[t]))] != 0:
                self.RVCX[t] = self.RVCX[t-1]
                self.RVCY[t] = self.RVCY[t-1]
        return self.RVCX, self.RVCY
             
def RandomRVC(RVCX, RVCY, Grid, Tmax, d, N, eff):
    Frames = [Grid.copy()]
    for i in range(Tmax+2):
        A = Frames[i].copy()    
        if A[N-int(np.floor(RVCY[i]))-1][int(np.floor(RVCX[i]))] != 0: #If the new position for the robot isn't clean
            A[N-int(np.floor(RVCY[i]))-1][int(np.floor(RVCX[i]))] += -d*eff #Clean the element
            Frames.append(A.copy())
        else:
            Frames.append(A.copy())
    return Frames

 #For when i need to run it a bunch of times, and don't care about the path
def QuickRandomRVC(RVCX, RVCY, Grid, Tmax, d, N, eff):
    Dust = []
    for i in range(Tmax+1):
        if Grid[N-int(np.floor(RVCY[i]))-1][int(np.floor(RVCX[i]))] >= eff:
            Grid[N-int(np.floor(RVCY[i]))-1][int(np.floor(RVCX[i]))] += -d*eff
        if Grid[N-int(np.floor(RVCY[i]))-1][int(np.floor(RVCX[i]))] < eff:
            Grid[N-int(np.floor(RVCY[i]))-1][int(np.floor(RVCX[i]))] = 0
        Dust.append(Grid.sum())
    return Dust

def SEMfunc(eff):
    RandomData = np.genfromtxt("RandomRVCeff"+str(eff)+".csv", delimiter = ',', skip_header = 1)
    SnakeData = np.genfromtxt("SnakeRVCeff"+str(eff)+".csv", delimiter = ',', skip_header = 1)

    RandomErr, Random = [], []
    for i in range(len(RandomData[0])):
        SEM = []
        for j in range(len(RandomData)):
            SEM.append(RandomData[j][i])
        RandomErr.append(np.std(SEM)/np.sqrt(len(SEM)))
        Random.append(np.mean(SEM))
    return RandomErr, Random, SnakeData

def SEMfuncSpeed(drmax):
    RandomData = np.genfromtxt("RandomRVCdrmax"+str(drmax)+".csv", delimiter = ',', skip_header = 1)

    RandomErr, Random = [], []
    for i in range(len(RandomData[0])):
        SEM = []
        for j in range(len(RandomData)):
            SEM.append(RandomData[j][i])
        RandomErr.append(np.std(SEM)/np.sqrt(len(SEM)))
        Random.append(np.mean(SEM))
    return RandomErr, Random
        
def SEMfuncObstacle(eff, N):
    RandomData = np.genfromtxt("RandomRVCeff"+str(eff)+"N"+str(N)+".csv", delimiter = ',', skip_header = 1)
    SnakeData = np.genfromtxt("SnakeRVCeff"+str(eff)+"N"+str(N)+".csv", delimiter = ',', skip_header = 1)
    RandomDataObstacle = np.genfromtxt("RandomRVCObstacleeff"+str(eff)+"N"+str(N)+".csv", delimiter = ',', skip_header = 1)
    RandomErr, Random, RandomObstacle, RandomErrObstacle = [], [], [], []
    for i in range(len(RandomData[0])):
        SEM, SEMObstacle = [], []
        for j in range(len(RandomData)):
            SEM.append(RandomData[j][i])
            SEMObstacle.append(RandomDataObstacle[j][i])
        RandomErr.append(np.std(SEM)/np.sqrt(len(SEM)))
        Random.append(np.mean(SEM))
        RandomErrObstacle.append(np.std(SEMObstacle)/np.sqrt(len(SEMObstacle)))
        RandomObstacle.append(np.mean(SEMObstacle))
    return RandomErr, Random, RandomErrObstacle, RandomObstacle, SnakeData