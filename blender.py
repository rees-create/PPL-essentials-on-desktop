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

    def __call__(self, start, step):
        """
        :param Vector2 start: Where to start iterating from
        :param int step: How much to step forward
        """
        self.current = Vector2(start)
        self.step = step
        return self
    
    def __iter__(self):
        x = int(np.abs(self.x)) if self.x != 0 else 1
        y = int(np.abs(self.y)) if self.y != 0 else 1
        xd = int(self.x - self.current.x)
        xd = 1 if xd == 0 else xd

        yd = int(self.y - self.current.y)
        yd = 1 if yd == 0 else yd

        uxd = int(np.abs(self.x - self.current.x)) 
        uxd = 1 if uxd == 0 else uxd

        uyd = int(np.abs(self.y - self.current.y))
        uyd = 1 if uyd == 0 else uyd

        

        self._umin_ = min(uxd, uyd) # (unsigned minimum) smallest unsigned DIFFERENCE
        self._umax_ = max(uxd, uyd) # (unsigned maximum) biggest unsigned DIFFERENCE

        self._diagonal_moves_ = self._umin_
        
        self._diag_ = Vector2(self.step, self.step) * (Vector2(xd, yd) / Vector2(uxd, uyd)) 
        self.current = self.current - self._diag_ 

        #to avoid errors on length 1 iterations, use expression instead of {uxd, uyd}
        if int(np.abs(self.x - self.current.x)) > int(np.abs(self.y - self.current.y)):
            self._straight_ = Vector2(self.step, 0) * (Vector2(xd, yd) / Vector2(uxd, uyd))
        elif int(np.abs(self.y - self.current.y)) > int(np.abs(self.x - self.current.x)):
            self._straight_ = Vector2(0, self.step) * (Vector2(xd, yd) / Vector2(uxd, uyd))
        else:
            self._straight_ = self._diag_
        
        #destiny calculation (direction to destination vector), tracks overshoots on flip
        straight_distance = self.current + self._straight_  - self
        usdx = np.abs(straight_distance.x) if np.abs(straight_distance.x) > 0 else 1 #zero-guarded, unsigned straight_distance.x
        usdy = np.abs(straight_distance.y) if np.abs(straight_distance.y) > 0 else 1
        self._straight_destiny_ = straight_distance / Vector2(usdx, usdy) # 1-vector form of straight_distance

        diagonal_distance = self.current + self._diag_ - self
        uddx = np.abs(diagonal_distance.x) if np.abs(diagonal_distance.x) > 0 else 1
        uddy = np.abs(diagonal_distance.y) if np.abs(diagonal_distance.y) > 0 else 1
        self._diagonal_destiny_ = diagonal_distance / Vector2(uddx, uddy)

        return self
    
    def __next__(self):
        if self._umax_ == self._diagonal_moves_:
            self._diagonal_moves_ -= 1
        
        straight_distance = self.current + self._straight_  - self
        usdx = np.abs(straight_distance.x) if np.abs(straight_distance.x) > 0 else 1 #zero-guarded, unsigned straight_distance.x
        usdy = np.abs(straight_distance.y) if np.abs(straight_distance.y) > 0 else 1
        straight_destiny = straight_distance / Vector2(usdx, usdy) # 1-vector form of straight_distance

        diagonal_distance = self.current + self._diag_ - self
        uddx = np.abs(diagonal_distance.x) if np.abs(diagonal_distance.x) > 0 else 1
        uddy = np.abs(diagonal_distance.y) if np.abs(diagonal_distance.y) > 0 else 1
        diagonal_destiny = diagonal_distance / Vector2(uddx, uddy)
        #print(f'straight end? = {straight_distance != Vector2(0,0)} no overshoot? = {self._straight_ != straight_destiny}')
        #print(f'straight distance = {straight_distance}, straight_destiny = {straight_destiny}, diagonal_destiny = {diagonal_destiny}')
        diag_destiny_match = self._diagonal_destiny_ == diagonal_destiny
        if (self._diagonal_moves_ > 0): # \
                                             #or self.current + self._diag_ != self:
            self.current = self.current + self._diag_
            self._diagonal_moves_ -= 1
            return self.current
        
        elif self._straight_destiny_ == straight_destiny:
            self.current = self.current + self._straight_
            print(f'straight moves = {self._straight_}, diag = {self._diag_}')
            return self.current
        else:
            raise StopIteration

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

vec = Vector2(-2,0)
for v in vec((0, 0), 1):
    print(v) # make Vector2 comparable before this can work.
    

class Blender:
    class Cell:
        def __init__(self, grid, start_percent, module_res, coords, position, percentPosition):
            size = Vector2((start_percent * module_res.x) / 2, (start_percent * module_res.y) / 2)
            self.discreteSize = size - size.remainders()
            self.coords = coords        
    
    def __init__(self, module, coords, start_percent, module_res, _9slice_res):
        self.module = module
        self.coords = coords
        self.start_percent = start_percent
        self.module_res = module_res
        self._9slice_res = _9slice_res
        # auto-calculated params
        self.grid = []
        #
        #self.discreteSize = size - size.remainders()
        
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
                
    
   

    




                    
