class Grid:
    def __init__(self):
        self.grid = [
            ['M','L','P','M','M'],
            ['P','P','P','M','L'],
            ['P','M','P','M','L'],
            ['M','M','P','L','L'],
            ['M','M','L','P','L']
        ]
        self.memo = []
    
    def append_to_memo(self, id_, index):
        #Thus the first element is an integer and the rest are 2-tuples.
        for ch_index in range(len(self.memo)):
            if self.memo[ch_index][0] == id_:
                self.memo[ch_index].append(index)
                return
        self.memo.append([id_, index])
    def in_memo(self, index):
        found = False 
        found_count = 0
        idx = 0
        for el in self.memo:
            if index[0] == el[0] and index[1] == el[1]:
                found = True
                found_count += 1
            #if found_count > 1:
            #    del self.memo[idx]
            #    found_count -= 1
            #idx += 1
        return found
    

grid = Grid()

add = lambda tup1, tup2: (tup1[0] + tup2[0], tup1[1] + tup2[1])
neg = lambda tup: (-tup[0], -tup[1])
def findAdjacents(index, diagonals = True):
    #vertical, diagonal, horizontal, negative diagonal
    directions = [(1,0), (1,-1), (0,1), (1,1)] if diagonals else [(1,0), (0,1)]
    adjacent_indices = [(0,0) for i in range(2 * len(directions))]
    
    for i in range(len(directions)):
        idx_pos = (i, i+len(directions)) if i // 2 % 2 == 0 else (i+len(directions), i)
        adjacent_indices[idx_pos[0]] = add(index, neg(directions[i]))
        adjacent_indices[idx_pos[1]] = add(index, directions[i])
    return adjacent_indices

#print(findAdjacents((2,2), False))


#def in_memo(grid, index):
#    for chains in grid.memo:
#        for idx in chains:
#            if idx == index:
#                return True
#    return False

def grid_size(arr):
    size = 0
    for row in arr:
        for el in row:
            size += 1
    return size


def calcEqualityChains(grid, current_index = (0,0), prop_dir = (0,1), diagonals = True):
    #if memo reaches grid size kill all recursive calls
    if grid_size(grid.grid) == len(grid.memo):
        print(grid.memo)
        return
    
    #generate adjacent list for current index, add current index to memo
    adjacent_list = findAdjacents(current_index, diagonals)
    if not grid.in_memo(current_index):
        grid.memo.append(current_index)
    
    #list for equal adjacents, bool indicator of no equal adjacents
    equalAdjacents = []
    none_equal = True
    
    #find equal adjacents (that aren't memoized already)
    for adjacent in adjacent_list:
        if (-1 in adjacent) or (-1 in current_index) or (len(grid.grid) in adjacent) \
          or (len(grid.grid) in current_index): #for uneven grid use of len(grid) doesn't suffice.
            continue

        #check for equality
        if grid.grid[adjacent[0]][adjacent[1]] == grid.grid[current_index[0]][current_index[1]]:
            #collect equal adjacents
            equalAdjacents.append(adjacent)
            none_equal = False 
    
    #recursively visit non-memoized equal adjacents
    _else = 0 #Acting as an "else" for the loop, if it's equal to the loop length,
    #then the loop was skipped and none_equal should be True.
    for eq_adjacent in equalAdjacents:
        if grid.in_memo(eq_adjacent):
            _else += 1
            continue
        #smooth recursion
        calcEqualityChains(grid, eq_adjacent, prop_dir, diagonals)
    if _else == len(equalAdjacents):
        none_equal = True
        
    #propagate to nearest cell recursively if no indices are found (if loop above doesn't run)
    if none_equal:
        next_index = add(current_index, prop_dir)
        valid_used_adjacents = []
        for adjacent in adjacent_list: #rotate until indices are valid
            possible_prop_indices = [adjacent, current_index, next_index]
            invalid_possible_prop_index = (True in [-1 in idx for idx in possible_prop_indices]) \
                                          or (True in [len(grid.grid) in idx for idx in possible_prop_indices])
            if not(invalid_possible_prop_index):
                #rotate   
                prop_dir = add(adjacent, neg(current_index))
                next_index = add(current_index, prop_dir)
                if not grid.in_memo(next_index):
                    calcEqualityChains(grid, next_index, prop_dir, diagonals)
        #INTENTIONALLY look through grid for indices not in memo, recurse.
        if grid_size(grid.grid) > len(grid.memo):
            for row in range(len(grid.grid)):
                for col in range(len(grid.grid[0])):
                    index = (row, col)
                    last_index = grid.memo[len(grid.memo) - 1]
                    if not grid.in_memo(index) and \
                    not(grid.grid[row][col] == grid.grid[last_index[0]][last_index[1]]):
                        #jump
                        calcEqualityChains(grid, index, prop_dir, diagonals)
    
    
            
calcEqualityChains(grid, current_index = (4,2), prop_dir = (0,1), diagonals = True)
print(grid.memo)
