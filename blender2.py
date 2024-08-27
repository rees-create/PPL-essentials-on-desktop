import numpy as np
import math
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

XSize = 1024
ZSize = 1024
class Module:
    def __init__(self, nOctaves, wavelength, persistence, lacunarity, numWaves) -> None:
        #self.dimensions = dimensions
        self.nOctaves = nOctaves
        self.wavelength = wavelength
        self.persistence = persistence
        self.lacunarity = lacunarity
        self.numWaves = numWaves

dimensions = (512, 512)

module1 = Module(    
    #dimensions = (128, 128),
    nOctaves = 4,
    wavelength = 0.23,
    persistence = 0.08,
    lacunarity = 3.14,
    numWaves = 4
)

module2 = Module(    
    #dimensions = (128, 128),
    nOctaves = 4,
    wavelength = 0.43,
    persistence = 0.015,
    lacunarity = 5.9,
    numWaves = 4
)

#first step is to plot a gridded terrain, next is to blend
def DirectionalXY(x, y, numDirections, waveIndex):
    axis = (waveIndex * math.pi) / numDirections
    xComponent = x * math.sin(axis)
    yComponent = y * math.cos(axis)
    return xComponent + yComponent

def Index1Dto2D(index, n_cols):
    col = index % n_cols
    row = index // n_cols
    return (col, row)

def SetHeightmapPartition(partitionDimensions, moduleGrid):
    terrainHeightmap = [[0 for x in range(XSize)] for z in range(ZSize)]
    partitionX, partitionY = partitionDimensions

    if ZSize % partitionY != 0:
        raise ValueError("Your partition dimensions don't factor correctly")
    
    for gridIndex1D in range(partitionX * partitionY):
        if gridIndex1D >= (partitionX * partitionY):
            raise ValueError("gridIndex1D out of range. Make sure gridIndex1D < partitionDimensions.x * partitionDimensions.y")
        
        XStep = XSize // partitionX
        ZStep = ZSize // partitionY
        moduleX, moduleY = Index1Dto2D(gridIndex1D, partitionX)
        module = moduleGrid[moduleX][moduleY]
        moduleHeightmap = ppl_heightmap(
            (XStep, ZStep),
            module.nOctaves,
            module.wavelength,
            module.persistence,
            module.lacunarity,
            module.numWaves
        )
        startingCoords = Index1Dto2D(gridIndex1D, partitionY)
        startingCoords = (startingCoords[0] * XStep, startingCoords[1] * ZStep)

        for z in range(startingCoords[1], startingCoords[1] + ZStep):
            for x in range(startingCoords[0], startingCoords[0] + XStep):
                zeroedIndex = (x - startingCoords[0], z - startingCoords[1])
                terrainHeightmap[z][x] = moduleHeightmap[zeroedIndex[1]][zeroedIndex[0]]

    return terrainHeightmap


def ppl_heightmap(dimensions, nOctaves, wavelength, persistence, lacunarity, numWaves):
    heightmap = [[0 for y in range(dimensions[1])] for x in range(dimensions[0])]

    for z in range(dimensions[1]):
        for x in range(dimensions[0]):
            pointValue = 0
            for waveIndex in range(1, numWaves + 1):
                ensembleValue = 0
                for _octave in range (1, nOctaves + 1):
                    alpha = persistence ** _octave
                    omega = (1/wavelength) * (lacunarity ** _octave)
            
                    directionalTheta = DirectionalXY(x, z, numWaves, waveIndex)
                    directionalSizeTheta = dimensions[0] + dimensions[1]
                    
                    if directionalSizeTheta == 0:
                        directionalSizeTheta = 1
                    
                    parameter = omega * directionalTheta / directionalSizeTheta
                    wave = math.sin(parameter)
                    octave = alpha * wave
                    ensembleValue += octave
                    
                pointValue += ensembleValue
                
            heightmap[z][x] = pointValue
            
    return heightmap #let's not normalize for now; it's definitely something to consider though

#ok let's plot to test
#initialize graph
fig = plt.figure()
ax = fig.add_subplot(projection = '3d')
#make data

x = np.linspace(0, 1, XSize) #np.array(range(XSize))/XSize
y = np.linspace(0, 1, ZSize)#np.array(range(ZSize))/module.dimensions[1]
X, Y = np.meshgrid(x,y)
#z = np.array([0 for i in range(dimensions[1])])
terrainGrid = [
    [module1, module1, module2, module1],
    [module1, module2, module2, module2],
    [module1, module1, module2, module1],
    [module1, module2, module2, module2]
]

#z = np.array([[0 for _x in range(dimensions[0])] for _y in range(dimensions[1])])
partitionDimensions = (len(terrainGrid), len(terrainGrid[0]))
z = SetHeightmapPartition(partitionDimensions, terrainGrid)
z = np.array(z)

ax.plot_surface(X, Y, z)
ax.set_aspect('equal')

plt.show()
