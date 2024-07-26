import math
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d, interp1d
import numpy as np


#ya we codin' blender version paperplanedev up here... nah just the blending algorithm
class Vector2:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __init__(self, tup):
        self.x = tup[0]
        self.y = tup[1]
    
    def veclist(self, tup_list):
        vec_list = []
        for tup in tup_list:
            vec_list.append(Vector2(tup[0], tup[1]))
        return vec_list
        
    def __repr__(self):
        return str((self.x, self.y))
    
    def __add__(self, other):
        return (self.x + other.x, self.y + other.y)
    
    def __sub__(self, other):
        return (self.x - other.x, self.y - other.y)
    
    def __truediv__(self, other):
        return (self.x / other.x, self.y / other.y)
    
    def to_tup(self):
        return (self.x, self.y)
    
    def remainders(self):
        return Vector2(self.x % 1, self.y % 1)
    
    def magnitude(self):
        pass

class Blender:
    def __init__(self, module, coords, start_percent, module_res, _9slice_res):
        self.module = module
        self.coords = coords
        self.start_percent = start_percent
        self.module_res = module_res
        self._9slice_res = _9slice_res
        # auto-calculated params
        size = Vector2((start_percent * module_res.x) / 2, (start_percent * module_res.y) / 2)
        self.discreteSize = size.remainders()
        self.position = Vector2()

def findAdjacents(index, diagonals = True):
    #vertical, diagonal, horizontal, negative diagonal
    directions = Vector2.veclist([(1,0), (1,-1), (0,1), (1,1)]) if diagonals else Vector2.veclist([(1,0), (0,1)])
    adjacent_indices = [(0,0) for i in range(2 * len(directions))]
    
    for i in range(len(directions)):
        idx_pos = (i, i+len(directions)) if i // 2 % 2 == 0 else (i+len(directions), i)
        adjacent_indices[idx_pos[0]] = index - directions[i]
        adjacent_indices[idx_pos[1]] = index + directions[i]
    return adjacent_indices

def DirectionalXY(x, y, numDirections, waveIndex):
    axis = (waveIndex * math.pi) / numDirections
    xComponent = x * math.sin(axis)
    yComponent = y * math.cos(axis)
    return xComponent + yComponent

def ppl_heightmap(dimensions, nOctaves, wavelength, persistence, lacunarity, numWaves):
    heightmap = [[0 for y in range(len(dimensions.y))] for x in range(len(dimensions.x))]

    for z in range(len(dimensions.y)):
        for x in range(len(dimensions.x)):
            pointValue = 0
            for waveIndex in range(1, numWaves + 1):
                ensembleValue = 0
                for _octave in range (1, nOctaves + 1):
                    alpha = persistence ** _octave
                    omega = (1/wavelength) * (lacunarity ** _octave)
            
                    directionalTheta = DirectionalXY(x, z, numWaves, waveIndex)
                    directionalSizeTheta = dimensions.x + dimensions.y
                    
                    if directionalSizeTheta == 0:
                        directionalSizeTheta = 1
                    
                    parameter = omega * directionalTheta / directionalSizeTheta
                    wave = math.sin(parameter)
                    octave = alpha * wave
                    ensembleValue += octave
                    
                pointValue += ensembleValue
                
            heightmap[z][x] = pointValue
            
    return heightmap #let's not normalize for now; it's definitely something to consider though

def weightInterpolation(terrainGrid, spline_type, _9slice_res, start_percent, direction):
    for module in terrainGrid:
        _9slice_size = (_9slice_res * 4 + 4, _9slice_res * 4 + 4)
        #Remember, while on a module, it's the center of the intra-modular 9-slice.
        #intra_modular_slice_coords = findAdjacents(module.moduleIndex) 
        standard_slice = findAdjacents((0,0))
        intra_modular_slice = [[[] for y in range(3)] for x in range(3)]
        #unpack standard_slice into intra_modular_slice
        for cell in standard_slice:
            intra_modular_slice[cell[0]][cell[1]].append #append inter-modular 9-slice information (class object)
                
    
   

    




                    
