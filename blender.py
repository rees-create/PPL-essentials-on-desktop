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

    def __call__(self, *start, step = 1):
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
        def __init__(self, blender, start_percent, module_res, coords):
            size = Vector2((start_percent * module_res.x) / 2, (start_percent * module_res.y) / 2)
            rem = blender.add_remainder(size.remainders())
            self.discreteSize = size - rem
            self.coords = coords
            self.position = sum([blender.grid[cell.y][cell.x].discreteSize for cell in Vector2(0,0)(coords.to_tup())])
            self.percentPosition = self.position / module_res
            self.spline_point = None
            
    
    def __init__(self, coords, start_percent, module_res, _9slice_res, spline_points, interpolation_type = "linear"):
        self.coords = coords
        self.start_percent = start_percent
        self.module_res = module_res
        self._9slice_res = _9slice_res
        self.interpolation_type = interpolation_type
        # force len of spline points to be correct
        while len(spline_points) < _9slice_res * 2:
            spline_points.append(0.0) 
        while len(spline_points) > _9slice_res * 2:
            spline_points.pop()

        self.spline_points = spline_points
        # auto-calculated params
        _9slice_size = (_9slice_res * 2 + 1, _9slice_res * 2 + 1)
        # to get x and y *coordinates*, subtract _9slice_res + 1 from *indices* x and y.
        self.center = Vector2(self._9slice_res, self._9slice_res)
        self.center_grid = [[self.Cell(self, start_percent, module_res, Vector2(x - _9slice_res - 1, y - _9slice_res - 1)) 
                                for y in range(_9slice_size[1])] for x in range(_9slice_size[0])]
        #this is just the inter-modular 9-slice tho, you have to make the intra-modular one using findAdjacents
        self.adjacents = findAdjacents(0,0)
        self.adjacents = np.array([Vector2(adjacent) for adjacent in self.adjacents]) * (2 * _9slice_res + 1)
        
        self.adjacent_grids = [[[self.Cell(self, start_percent, module_res, Vector2(x + adj[0] - _9slice_res - 1, y + adj[1] - _9slice_res - 1)) 
                                for y in range(_9slice_size[1])] for x in range(_9slice_size[0])] for adj in self.adjacents]
        
        self.remainders = 0

    def add_remainder(self, remainder):
        self.remainders += remainder
        if self.remainders >= 1:
            self.remainders -= 1
            return 1
        return 0
    
    def set_spline_points(self):
        # set adjacents, iterate up the center grid, then the adjacent grid
        # setting adjacents
        adjacents = np.array(findAdjacents((0,0))) * self._9slice_res #the plus 1 to iterate to the end
        adjacents = [Vector2(adj) for adj in adjacents]
        # iterating up center grid
        for center_edge in adjacents:
            for index, center_cell in enumerate(center_edge(0,0)):
                if center_cell != Vector2(0,0):
                    self.center_grid[center_cell.y][center_cell.x].spline_point = self.spline_points[index]
        # iterating up adjacent grid
        for adjacent_edge in range(len(adjacents)):
            center_facing_edge = Vector2(adjacents[(adjacent_edge + 4) % 8]) #find opposite adjacent
            for index, adjacent_cell in enumerate(adjacents[adjacent_edge](center_facing_edge.to_tup())):
                self.adjacent_grid[index][adjacent_cell.y][adjacent_cell.x].spline_point = self.spline_points[index + 2]
        # bad code alert: you're iterating a vector, so index will be the vector iteration step, not cyclic step

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
        
        #Remember, while on a module, it's the center of the intra-modular 9-slice.
        coords_list = [module.coords.x, module.coords.y]
        invalid_module = 0 in coords_list or terrainGrid.dimens in coords_list
        if not invalid_module:
            module.blender = Blender(module.coords, start_percent, module.resolution, _9slice_res)

   

    




                    
