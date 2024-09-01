import numpy as np
import math
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

XSize = 1024
ZSize = 1024 #in normal situations (including matplotlib) this should by YSize, but Z is used just for 
             #Unity 3D space conventions. Don't get confused by this
class Module:
    def __init__(self, nOctaves, wavelength, persistence, lacunarity, numWaves) -> None:
        #self.dimensions = dimensions
        self.nOctaves = nOctaves
        self.wavelength = wavelength
        self.persistence = persistence
        self.lacunarity = lacunarity
        self.numWaves = numWaves

    def single_point(self, point):
        return ppl_heightmap(
            point,
            self.nOctaves,
            self.wavelength,
            self.persistence,
            self.lacunarity,
            self.numWaves,
            single_point = True
        )


#dimensions = (512, 512)

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


def ppl_heightmap(dimensions, nOctaves, wavelength, persistence, lacunarity, numWaves, single_point = False):
    if not single_point:
        heightmap = [[0 for y in range(dimensions[1])] for x in range(dimensions[0])]

    for z in range(dimensions[1]):
        for x in range(dimensions[0]):
            pointValue = 0
            for waveIndex in range(1, numWaves + 1):
                ensembleValue = 0
                for _octave in range (1, nOctaves + 1):
                    alpha = persistence ** _octave
                    omega = (1/wavelength) * (lacunarity ** _octave)
                    if single_point:
                        directionalTheta = DirectionalXY(dimensions[0], dimensions[1], numWaves, waveIndex)
                    else:
                        directionalTheta = DirectionalXY(x, z, numWaves, waveIndex)
                    directionalSizeTheta = dimensions[0] + dimensions[1]
                    
                    if directionalSizeTheta == 0:
                        directionalSizeTheta = 1
                    
                    parameter = omega * directionalTheta / directionalSizeTheta
                    wave = math.sin(parameter)
                    octave = alpha * wave
                    ensembleValue += octave
                    if single_point:
                        break
                pointValue += ensembleValue
                if single_point:
                    return pointValue
            heightmap[z][x] = pointValue
            
    return heightmap #let's not normalize for now; it's definitely something to consider though

#plotting
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

#blend here
# [REMINDER]: in normal situations (including matplotlib) the Zs in blending should be Ys, but Z is used just for 
# Unity 3D space conventions. Don't get confused by this
thiccness_ratio = 4
weightDividers = ([0.2, 0.4], [0.5, 0.8])
fullInterpDividers = ([0] + weightDividers[0] + [1], [0] + weightDividers[1] + [1]) 
nDividersX, nDividersZ = (len(weightDividers[0]) + 2, len(weightDividers[1]) + 2) # remember that outer default dividers 0 and 1 will be added 

biomesDimens = (len(terrainGrid[0]), len(terrainGrid))  #biomes: another name for terrainGrid, small b for biome :)
bDX, bDZ = biomesDimens
divIntervalX, divIntervalZ = XSize / bDX, ZSize / bDZ
thiccnessX, thiccnessZ = (divIntervalX / thiccness_ratio, divIntervalZ / thiccness_ratio)
dividerX, dividerZ = ([XSize / bDX * i for i in range(1, bDX + 1)], [ZSize / bDZ * i for i in range(1, bDZ + 1)])

blendLinesX = np.array([line for line in np.linspace(np.array(dividerX) - thiccnessX, np.array(dividerX) + thiccnessX, nDividersX)]).T
blendLinesZ = np.array([line for line in np.linspace(np.array(dividerZ) - thiccnessZ, np.array(dividerZ) + thiccnessZ, nDividersZ)]).T

class BlendDivider:
    def __init__(self, position, weight, zero_index):
        self.position = position
        self.weight = weight
        self.zero_index = zero_index
    def __repr__(self):
        return f'<BlendDivider> (position: {self.position}, weight: {self.weight}, zero_index: {self.zero_index})'
    def createBlendTable(blendLinesX, blendLinesZ):
        fullbDX, fullbDY = fullInterpDividers
        x_list = [zip(blendLineX, fullbDX, range(len(fullbDX))) for blendLineX in blendLinesX]
        z_list = [zip(blendLineZ, fullbDY, range(len(fullbDY))) for blendLineZ in blendLinesZ]
        blendWallsX = []
        blendWallsZ = []
        for wall_x, wall_z in zip(x_list, z_list):
            blendDivsX = []
            blendDivsZ = []
            for div_x, div_z in zip(wall_x, wall_z):
                position_x, weight_x, zeroIndex_x = div_x
                position_z, weight_z, zeroIndex_z = div_z
                blendDivsX.append(BlendDivider(position_x, weight_x, zeroIndex_x))
                blendDivsZ.append(BlendDivider(position_z, weight_z, zeroIndex_z))
            blendWallsX.append(blendDivsX)
            blendWallsZ.append(blendDivsZ)

        return blendWallsX, blendWallsZ
    
    def blended_point(self, start_divider, next_divider, points, interpolation_mode, t):
        if t < 0 or t > 1:
            return ValueError('t should be between 0 and 1')
        interpolation_modes = ['quadratic', 'cubic', 'sigmoid', 'linear']
        interpolation_functions = [lambda start, end: start + (end - start) * (t**2),
                                   lambda start, end: start + (end - start) * (t**3),
                                   lambda start, end: start + (end - start) * (1/(1 - np.exp(t)))] #this will break, get a normalized shifted function
        
        mode = interpolation_modes.index(interpolation_mode)
        erp = interpolation_functions[mode]

        weight = erp(start_divider.weight, next_divider.weight)
        local = points[0] * (1 - weight) 
        other = points[1] * weight 
        return (local + other) / 2

blendWallsX, blendWallsZ = BlendDivider.createBlendTable(blendLinesX, blendLinesZ)
for _y in ZSize: #I'm done with the Unity convention at this point
    for _x in XSize:
        heightmap_val = z[_y][_x]
        index = 0
        weight = 0
        for blendWallX, blendWallY in zip(blendWallsX, blendWallsZ): 
            # Yeah, the Unity convention is no longer in use here to keep the code readable.
        
            for divider_idx in range(len(blendWallX)):
                diff = _x - math.floor(blendWallX[divider_idx].position) #interpolant
                if diff > 0:
                    index = blendWallX[divider_idx].zero_index
                    weight = blendWallY[divider_idx].weight
                    break
            break
        points = (heightmap_val * (1 - weight), terrainGrid[index + 1].single_point((_x,_y)))
        


print(f'dividerX = {dividerX}')
print(f'blendLinesX = {blendLinesX}')
#print(blendTable)

ax.plot_surface(X, Y, z)
ax.set_aspect('equal')

plt.show()
