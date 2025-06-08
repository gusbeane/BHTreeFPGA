import numpy as np
import copy
MAX_DEPTH = 10

# Rotation table for 3D space
rottable3 = np.array([
    [36, 28, 25, 27, 10, 10, 25, 27], [29, 11, 24, 24, 37, 11, 26, 26], [8, 8, 25, 27, 30, 38, 25, 27],
    [9, 39, 24, 24, 9, 31, 26, 26], [40, 24, 44, 32, 40, 6, 44, 6], [25, 7, 33, 7, 41, 41, 45, 45],
    [4, 42, 4, 46, 26, 42, 34, 46], [43, 43, 47, 47, 5, 27, 5, 35], [33, 35, 36, 28, 33, 35, 2, 2],
    [32, 32, 29, 3, 34, 34, 37, 3], [33, 35, 0, 0, 33, 35, 30, 38], [32, 32, 1, 39, 34, 34, 1, 31],
    [24, 42, 32, 46, 14, 42, 14, 46], [43, 43, 47, 47, 25, 15, 33, 15], [40, 12, 44, 12, 40, 26, 44, 34],
    [13, 27, 13, 35, 41, 41, 45, 45], [28, 41, 28, 22, 38, 43, 38, 22], [42, 40, 23, 23, 29, 39, 29, 39],
    [41, 36, 20, 36, 43, 30, 20, 30], [37, 31, 37, 31, 42, 40, 21, 21], [28, 18, 28, 45, 38, 18, 38, 47],
    [19, 19, 46, 44, 29, 39, 29, 39], [16, 36, 45, 36, 16, 30, 47, 30], [37, 31, 37, 31, 17, 17, 46, 44],
    [12, 4, 1, 3, 34, 34, 1, 3], [5, 35, 0, 0, 13, 35, 2, 2], [32, 32, 1, 3, 6, 14, 1, 3],
    [33, 15, 0, 0, 33, 7, 2, 2], [16, 0, 20, 8, 16, 30, 20, 30], [1, 31, 9, 31, 17, 17, 21, 21],
    [28, 18, 28, 22, 2, 18, 10, 22], [19, 19, 23, 23, 29, 3, 29, 11], [9, 11, 12, 4, 9, 11, 26, 26],
    [8, 8, 5, 27, 10, 10, 13, 27], [9, 11, 24, 24, 9, 11, 6, 14], [8, 8, 25, 15, 10, 10, 25, 7],
    [0, 18, 8, 22, 38, 18, 38, 22], [19, 19, 23, 23, 1, 39, 9, 39], [16, 36, 20, 36, 16, 2, 20, 10],
    [37, 3, 37, 11, 17, 17, 21, 21], [4, 17, 4, 46, 14, 19, 14, 46], [18, 16, 47, 47, 5, 15, 5, 15],
    [17, 12, 44, 12, 19, 6, 44, 6], [13, 7, 13, 7, 18, 16, 45, 45], [4, 42, 4, 21, 14, 42, 14, 23],
    [43, 43, 22, 20, 5, 15, 5, 15], [40, 12, 21, 12, 40, 6, 23, 6], [13, 7, 13, 7, 41, 41, 22, 20]
], dtype=np.uint8)

# Subpixel lookup table for 3D space
subpix3 = np.array([
    [0, 7, 1, 6, 3, 4, 2, 5], [7, 4, 6, 5, 0, 3, 1, 2], [4, 3, 5, 2, 7, 0, 6, 1], [3, 0, 2, 1, 4, 7, 5, 6], [1, 0, 6, 7, 2, 3, 5, 4],
    [0, 3, 7, 4, 1, 2, 6, 5], [3, 2, 4, 5, 0, 1, 7, 6], [2, 1, 5, 6, 3, 0, 4, 7], [6, 1, 7, 0, 5, 2, 4, 3], [1, 2, 0, 3, 6, 5, 7, 4],
    [2, 5, 3, 4, 1, 6, 0, 7], [5, 6, 4, 7, 2, 1, 3, 0], [7, 6, 0, 1, 4, 5, 3, 2], [6, 5, 1, 2, 7, 4, 0, 3], [5, 4, 2, 3, 6, 7, 1, 0],
    [4, 7, 3, 0, 5, 6, 2, 1], [6, 7, 5, 4, 1, 0, 2, 3], [7, 0, 4, 3, 6, 1, 5, 2], [0, 1, 3, 2, 7, 6, 4, 5], [1, 6, 2, 5, 0, 7, 3, 4],
    [2, 3, 1, 0, 5, 4, 6, 7], [3, 4, 0, 7, 2, 5, 1, 6], [4, 5, 7, 6, 3, 2, 0, 1], [5, 2, 6, 1, 4, 3, 7, 0], [7, 0, 6, 1, 4, 3, 5, 2],
    [0, 3, 1, 2, 7, 4, 6, 5], [3, 4, 2, 5, 0, 7, 1, 6], [4, 7, 5, 6, 3, 0, 2, 1], [6, 7, 1, 0, 5, 4, 2, 3], [7, 4, 0, 3, 6, 5, 1, 2],
    [4, 5, 3, 2, 7, 6, 0, 1], [5, 6, 2, 1, 4, 7, 3, 0], [1, 6, 0, 7, 2, 5, 3, 4], [6, 5, 7, 4, 1, 2, 0, 3], [5, 2, 4, 3, 6, 1, 7, 0],
    [2, 1, 3, 0, 5, 6, 4, 7], [0, 1, 7, 6, 3, 2, 4, 5], [1, 2, 6, 5, 0, 3, 7, 4], [2, 3, 5, 4, 1, 0, 6, 7], [3, 0, 4, 7, 2, 1, 5, 6],
    [1, 0, 2, 3, 6, 7, 5, 4], [0, 7, 3, 4, 1, 6, 2, 5], [7, 6, 4, 5, 0, 1, 3, 2], [6, 1, 5, 2, 7, 0, 4, 3], [5, 4, 6, 7, 2, 3, 1, 0],
    [4, 3, 7, 0, 5, 2, 6, 1], [3, 2, 0, 1, 4, 5, 7, 6], [2, 5, 1, 6, 3, 4, 0, 7]
], dtype=np.uint8)

def generate_ph_key(pos, depth=MAX_DEPTH):
    pos = np.asarray(pos, dtype=np.uint32)
    x_in, y_in, z_in = pos

    # Assuming input coordinates are uint32 (32-bit).
    input_bit_width = pos.dtype.itemsize * 8

    # Calculate shift to extract the 'depth' most significant bits.
    # These shifted values will be treated as 'depth'-bit numbers.
    shift_amount = input_bit_width - depth

    if shift_amount > 0:
        x = x_in >> shift_amount
        y = y_in >> shift_amount
        z = z_in >> shift_amount
    elif shift_amount == 0: # depth == input_bit_width
        x, y, z = x_in, y_in, z_in
    else: # depth > input_bit_width, this case should ideally be handled by capping depth or ensuring valid inputs.
        raise ValueError(f"Depth {depth} is too large for the input type. Please cap depth or ensure valid inputs.")

    rotation = 0
    if depth <= 10:
        key = np.uint32(0)
    elif depth<=21:
        key = np.uint64(0)
    else:
        raise ValueError(f"Maximum depth supported is 21, requested depth of {depth}")

    # Loop 'depth' times. x, y, z are now effectively 'depth'-bit numbers.
    for i in range(depth):
        # Mask for the i-th MSB (from left, 0-indexed) of the 'depth'-bit numbers x,y,z.
        mask = 1 << (depth - 1 - i) 
        
        pix = (((x & mask) > 0) and 4 or 0) | \
              (((y & mask) > 0) and 2 or 0) | \
              (((z & mask) > 0) and 1 or 0)
        
        key = (key << 3) | int(subpix3[rotation][pix])
        rotation = int(rottable3[rotation][pix])

    return key

def generate_ph_key_from_list(pos_list, depth=MAX_DEPTH):
    keys = np.array([generate_ph_key(pos, depth) for pos in pos_list])
    return keys

class Node:
    def __init__(self, level, key, start_idx, Nleaf=1):
        if Nleaf != 1:
            raise ValueError('Nleaf must be 1')

        self.level = level
        self.key = key
        self.com = np.array([0.0, 0.0, 0.0])
        self.mass = 0.0
        self.start_idx = start_idx
        self.Npart = 0
        self.Nleaf = Nleaf

    def add_particle(self, pos_float, mass):
        self.com += pos_float * mass
        self.mass += mass
        self.Npart += 1

    def close(self, force_leaf=False):
        if self.Npart <= self.Nleaf or force_leaf:
            # del self.com
            # del self.mass
            self.is_leaf = True

        else:   
            self.com /= self.mass
            # del self.start_idx
            self.is_leaf = False
    
    def __repr__(self):
        return f"Node(level={self.level}, key={self.key}, start_idx={self.start_idx}, Npart={self.Npart}, Nleaf={self.Nleaf}, is_leaf={self.is_leaf})"

def generate_tree(pos_float, keys, mass, depth=MAX_DEPTH):
    tree = []
    node_temp = []

    # dummy for root node
    node_temp.append(None)

    # create empty nodes for first particle
    for j in range(1, depth+1):
        node_key = keys[0] >> 3*(depth - j)
        node = Node(j, node_key, 0)
        node.add_particle(pos_float[0], mass[0])
        node_temp.append(node)

    # create tree by looping through all particles
    for i, (posf, key) in enumerate(zip(pos_float, keys)):
        if i==0:
            continue

        for j in range(1, depth+1):
            node_key = key >> 3*(depth - j)

            if node_key == node_temp[j].key:
                node_temp[j].add_particle(posf, mass[i])

            else:
                # we need to close all the nodes below this level
                for k in range(j, depth+1):
                    force_leaf = k == depth
                    node_temp[k].close(force_leaf)

                    if node_temp[k].is_leaf:
                        k_leaf = k
                        # print(i, 'k_leaf', k_leaf)
                        break
                
                # now we add all of these nodes in reverse order and create new nodes
                for k in range(depth, j-1, -1):

                    if k <= k_leaf:
                        tree.append(copy.deepcopy(node_temp[k]))

                    key_for_level_k = keys[i] >> (3*(depth - k))

                    node_temp[k] = Node(k, key_for_level_k, i)
                    node_temp[k].add_particle(posf, mass[i])

                # now we can go to the next particle
                break
    

    # now finalize the tree
    # loop through all nodes and close them, stopping at the first leaf node
    for k in range(1, depth+1):
        force_leaf = k == depth-1
        node_temp[k].close(force_leaf)

        if node_temp[k].is_leaf:
            k_leaf = k
            break

    # now we add all of these nodes in reverse order and create new nodes
    for k in range(depth, 0, -1):
        if k <= k_leaf:
            tree.append(copy.deepcopy(node_temp[k]))

    # finally reverse the list
    tree.reverse()

    return tree

if __name__ == "__main__":
    data = np.load('nbody_golden.npz')
    pos = data['pos']
    pos_float = data['pos_float']
    pos_eval = data['pos_eval']
    pos_eval_float = data['pos_eval_float']
    acc = data['acc']
    M = np.full(pos.shape[0], data['M'])
    G = data['G']
    eps = data['eps']

    # sort in place
    keys = generate_ph_key_from_list(pos)
    sorted_indices = np.argsort(keys)
    pos = pos[sorted_indices]
    pos_float = pos_float[sorted_indices]
    acc = acc[sorted_indices]

    tree = generate_tree(pos_float, keys, M)