import math
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d, interp1d
import numpy as np
from functools import total_ordering


#ya we codin' blender version paperplanedev up here... nah just the blending algorithm
@total_ordering
class Vector2:
    def __init__(self, *args):
        if len(args) == 2:
            self.x = args[0]
            self.y = args[1]
        elif len(args) == 1:
            self.x = args[0][0]
            self.y = args[0][1]
    
    def __eq__(self, other):
        return other.x == self.x and other.y == self.y
    
    def __lt__(self, other):
        return self.magnitude() < other.magnitude()

    def __call__(self, start, step = 1):
        """
        :param Vector2 start: Where to start iterating from
        :param int step: How much to step forward
        """
        self.current = Vector2(start)
        self.step = int(np.abs(step))
        return self
    
    def __iter__(self):
        self._diff_ = self - self.current
        self._started_ = False
        return self
    
    def __next__(self):
       
        if self._started_ == False:
            self._started_ = True
            return self.current
        
        propagation = Vector2(0,0)
        if self._diff_.x != 0:
            propagation.x = -self.step if self._diff_.x + np.abs(self._diff_.x) == 0 else self.step
        if self._diff_.y != 0:
            propagation.y = -self.step if self._diff_.y + np.abs(self._diff_.y) == 0 else self.step
        
        self.current = self.current + propagation
        self._diff_ = self._diff_ - propagation

        if self._diff_ == Vector2(0, 0): 
            # susceptible to infinite loops if step is uneven (addition of step doesn't exactly equal _diff_ 
            # at some point). To prevent this, create a destiny (direction to destination) variable on both sides, 
            # terminate on case where original + updated destinies = Vector2(0,0). However, blender might not
            # have to use uneven step. Just worth noting 
            raise StopIteration
        
        return self.current

    def veclist(self, tup_list):
        vec_list = []
        for tup in tup_list:
            vec_list.append(Vector2(tup[0], tup[1]))
        return vec_list
        
    def __repr__(self):
        return str((self.x, self.y))
    
    def __add__(self, other):
        if isinstance(other, tuple):
            return Vector2(self.x + other[0], self.y + other[1])
        else:
            return Vector2(self.x + other.x, self.y + other.y)
    
    def __sub__(self, other):
        if isinstance(other, tuple):
            return Vector2(self.x - other[0], self.y - other[1])
        else:
            return Vector2(self.x - other.x, self.y - other.y)

    def __mul__(self, other):
        return Vector2(self.x * other.x, self.y * other.y)
    
    def __truediv__(self, other):
        if isinstance(other, tuple):
            return Vector2(self.x / other[0], self.y / other[1])
        else:
            return Vector2(self.x / other.x, self.y / other.y)
    
    def to_tup(self):
        return (self.x, self.y)
    
    def remainders(self):
        return Vector2(self.x % 1, self.y % 1)
    
    def magnitude(self):
        return math.sqrt(self.x**2 + self.y**2)

vec = Vector2(0,0)
for v in vec((6, 0), step = 1):
    print(v) # make Vector2 comparable before this can work.
    

class Blender:
    class Cell:
        def __init__(self, grid, start_percent, module_res, coords):
            size = Vector2((start_percent * module_res.x) / 2, (start_percent * module_res.y) / 2)
            rem = super.add_remainder(size.remainders())
            self.discreteSize = size - rem
            self.coords = coords
            self.position = sum([grid[cell.y][cell.x].discreteSize for cell in Vector2(0,0)(coords.to_tup())])
            self.percentPosition = self.position / module_res      
    
    def __init__(self, module, coords, start_percent, module_res, _9slice_res):
        self.module = module
        self.coords = coords
        self.start_percent = start_percent
        self.module_res = module_res
        self._9slice_res = _9slice_res
        # auto-calculated params
        self.grid = []
        self.remainders = 0
        #
        #self.discreteSize = size - size.remainders()    
        self.position = Vector2()

    def add_remainder(self, remainder):
        self.remainders += remainder
        if self.remainders >= 1:
            self.remainders -= 1
            return 1
        return 0

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
            intra_modular_slice[cell[1]][cell[0]].append #append inter-modular 9-slice information (class object)
                    
    
   

    




                    
