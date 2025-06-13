import numpy as np
import copy
MAX_DEPTH = 10

class PeanoHilbert:
    def __init__(self, dimensions=3, max_depth=MAX_DEPTH):
        self.dimensions = dimensions
        self.max_depth = max_depth
        
        if dimensions != 3:
            raise ValueError("Currently only 3D is supported")
            
        # Rotation table for 3D space
        self.rottable3 = np.array([
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
        self.subpix3 = np.array([
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

    def generate_key(self, pos, depth=None):
        """Generate Peano-Hilbert key for a single position"""
        if depth is None:
            depth = self.max_depth
            
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
            
            key = (key << 3) | int(self.subpix3[rotation][pix])
            rotation = int(self.rottable3[rotation][pix])

        return key

    def generate_keys(self, pos_list, depth=None):
        """Generate Peano-Hilbert keys for a list of positions"""
        if depth is None:
            depth = self.max_depth
            
        keys = np.array([self.generate_key(pos, depth) for pos in pos_list])
        return keys

# Legacy functions for backward compatibility
def generate_ph_key(pos, depth=MAX_DEPTH):
    """Legacy function - use PeanoHilbert class instead"""
    ph = PeanoHilbert(dimensions=3, max_depth=depth)
    return ph.generate_key(pos, depth)

def generate_ph_key_from_list(pos_list, depth=MAX_DEPTH):
    """Legacy function - use PeanoHilbert class instead"""
    ph = PeanoHilbert(dimensions=3, max_depth=depth)
    return ph.generate_keys(pos_list, depth)

class PointCloud:
    def __init__(self, pos, mass, eps, depth=MAX_DEPTH):
        # Store original integer positions
        self.pos_int = pos
        
        # Convert integer positions to floating point (normalized to [0,1])
        MAX_INT32 = np.iinfo(np.uint32).max
        self.pos_float = ((pos.astype(np.float64)) / MAX_INT32).astype(np.float32)
        
        # Create PeanoHilbert generator and compute keys
        self.ph_generator = PeanoHilbert(dimensions=3, max_depth=depth)
        self.keys = self.ph_generator.generate_keys(pos, depth)
        
        # Sort all attributes by the PH keys
        sorted_indices = np.argsort(self.keys)
        # sorted_indices = sorted_indices[::-1]
        
        # Apply sorting to all attributes
        self.pos = pos[sorted_indices]
        self.pos_float = self.pos_float[sorted_indices]
        self.keys = self.keys[sorted_indices]
        self.mass = mass[sorted_indices]
        self.eps = eps[sorted_indices] if hasattr(eps, '__len__') else eps

class NodeOrLeaf:
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
        self.is_leaf = False
        self.child_idx = -1
        self.Nchild = 0

    def add_particle(self, pos_float, mass, N=1):
        self.com += pos_float * mass
        self.mass += mass
        self.Npart += N

    def close(self, force_leaf=False):
        self.com /= self.mass
        self.is_leaf = self.Npart <= self.Nleaf or force_leaf

    def __repr__(self):
        return f"Node(level={self.level}, key={self.key}, start_idx={self.start_idx}, Npart={self.Npart}, Nleaf={self.Nleaf},\
 is_leaf={self.is_leaf}, child_idx={self.child_idx}, Nchild={self.Nchild})"

class BHTree:
    def __init__(self, point_cloud, depth=MAX_DEPTH):
        # Store reference to the point cloud
        self.point_cloud = point_cloud
        self.depth = depth
        
        # Generate the tree using the point cloud data
        self.tree = self._generate_tree(point_cloud.pos_float, point_cloud.keys, point_cloud.mass, depth)
    
    def _generate_tree(self, pos_float, keys, mass, depth=MAX_DEPTH, Nleaf=1):
        
        # --- Part 1 ---
        # create bottom level leaves. assumes keys are at the correct depth
        incoming_stream = []
        current_node = NodeOrLeaf(depth, keys[0], 0, Nleaf=Nleaf)
        current_node.add_particle(pos_float[0], mass[0])
        for i in range(1, len(keys)):
            if keys[i] == current_node.key:
                current_node.add_particle(pos_float[i], mass[i])
            else:
                # close current node and add to stream
                current_node.close(force_leaf=True)
                incoming_stream.append(current_node)
                
                # create new node
                current_node = NodeOrLeaf(depth, keys[i], i, Nleaf=Nleaf)
                current_node.add_particle(pos_float[i], mass[i])
        
        # close last node
        current_node.close(force_leaf=True)
        incoming_stream.append(current_node)

        # --- Part 2 ---
        # create tree bottom up, streaming nodes and tendrils to be merge sorted by the next kernel
        tree = []
        itree = 0
        for current_level in range(depth-1, -1, -1):
            print('processing level:', current_level, " incoming_stream len:", len(incoming_stream))
            print('incoming_stream:', incoming_stream[0])
            tendril_stream = []
            nodeleaf_stream = []

            current_children = []
            Nchild = 0

            # start with next node
            incoming_node = incoming_stream[0]
            current_key = incoming_node.key >> 3
            current_node = NodeOrLeaf(current_level, current_key, incoming_node.start_idx)
            current_node.add_particle(incoming_node.com*incoming_node.mass, incoming_node.mass, incoming_node.Npart)

            current_children.append(incoming_node)
            Nchild += 1

            # now process incoming stream
            for i in range(1, len(incoming_stream)):
                incoming_node = incoming_stream[i]
                current_key = incoming_node.key >> 3
                if current_key == current_node.key:
                    # we add the node to the parent node
                    current_node.add_particle(incoming_node.com * incoming_node.mass, incoming_node.mass, incoming_node.Npart)
                    current_children.append(incoming_node)
                    Nchild += 1
                else:
                    # emit children to final tree if they are not only children
                    # otherwise, promote child and place into tendril stream
                    if Nchild > 1:
                        current_node.child_idx = itree
                        current_node.Nchild = Nchild
                        for j in range(Nchild):
                            tree.append(current_children[j])
                            itree += 1
                        
                        # emit the node to the nodeleaf stream
                        current_node.close()
                        nodeleaf_stream.append(current_node)
                    else:
                        current_children[0].level -= 1
                        current_children[0].key >>= 3
                        tendril_stream.append(current_children[0])

                    

                    current_node = NodeOrLeaf(current_level, current_key, incoming_node.start_idx)
                    current_node.add_particle(incoming_node.com*incoming_node.mass, incoming_node.mass, incoming_node.Npart)

                    # reset children
                    current_children = [incoming_node]
                    Nchild = 1

            # now we close out the last node
            if Nchild > 1:
                current_node.child_idx = itree
                current_node.Nchild = Nchild
                for j in range(Nchild):
                    tree.append(current_children[j])
                    itree += 1

                # emit the node to the nodeleaf stream
                current_node.close()
                nodeleaf_stream.append(current_node)
                
            else:
                current_children[0].level -= 1
                current_children[0].key >>= 3
                tendril_stream.append(current_children[0])

            # now we need to merge the nodeleaf 
            incoming_stream = []
        
            if len(nodeleaf_stream) == 0:
                incoming_stream = tendril_stream
                continue

            if len(tendril_stream) == 0:
                incoming_stream = nodeleaf_stream
                continue

            i,j = 0,0
            while True:
                current_nodeleaf = nodeleaf_stream[i]
                current_tendril = tendril_stream[j]
                
                if current_nodeleaf.key < current_tendril.key:
                    incoming_stream.append(current_nodeleaf)
                    i += 1
                    if i == len(nodeleaf_stream):
                        while j < len(tendril_stream):
                            incoming_stream.append(tendril_stream[j])
                            j += 1
                        break

                elif current_nodeleaf.key > current_tendril.key:
                    incoming_stream.append(current_tendril)
                    j += 1
                    if j == len(tendril_stream):
                        while i < len(nodeleaf_stream):
                            incoming_stream.append(nodeleaf_stream[i])
                            i += 1
                        break
                else:
                    print(current_nodeleaf)
                    print(current_tendril)
                    raise ValueError("This should not occur!")
        
        # end
        # We don't add this because it is the root node, which we will always open
        # print('incoming_stream=', incoming_stream)
        # for node in incoming_stream:
            # node.child_idx = itree
            # tree.append(node)

        tree.reverse()
        for node in tree:
            node.child_idx = len(tree) - node.child_idx - node.Nchild

        return tree
    
    def walk_tree(self, pos_eval, G, theta=0.5):
        """
        Walk the tree to compute acceleration at pos_eval using Barnes-Hut algorithm
        
        Args:
            pos_eval: position where to evaluate acceleration
            G: gravitational constant
            theta: opening criterion parameter
            
        Returns:
            acc: computed acceleration
            Nint_leaf: number of leaf interactions
            Nint_node: number of node interactions
        """

        acc = np.zeros(3)
        Nint_leaf = 0
        Nint_node = 0
        
        
        


    # def walk_tree(self, pos_eval, G, theta=0.5):
    #     """
    #     Walk the tree to compute acceleration at pos_eval using Barnes-Hut algorithm
        
    #     Args:
    #         pos_eval: position where to evaluate acceleration
    #         G: gravitational constant
    #         theta: opening criterion parameter
            
    #     Returns:
    #         acc: computed acceleration
    #         Nint_leaf: number of leaf interactions
    #         Nint_node: number of node interactions
    #     """
    #     Nint_leaf = 0
    #     Nint_node = 0
    #     acc = np.zeros(3)
    #     pos0 = pos_eval
        
    #     # we start by opening all the top level nodes
    #     current_level = 1
    #     for tr in self.tree:
    #         # we are ignoring all nodes that are finer (higher) than the current level
    #         if tr.level > current_level:
    #             continue

    #         # if tr.level is smaller than current_level, it means that we've exhausted all the nodes at current level
    #         # and we are moving on to the next auncle/grand-auncle/... node
    #         current_level = tr.level

    #         # okay, so if we are at a leaf, we loop through all particles and add to our acc
    #         if tr.is_leaf:
    #             for i in range(tr.Npart):
    #                 rsq = np.sum((self.point_cloud.pos_float[tr.start_idx+i] - pos0)**2)
    #                 rcubed = rsq**(3./2.)
    #                 acc += G * self.point_cloud.mass[tr.start_idx+i] * (self.point_cloud.pos_float[tr.start_idx+i] - pos0) / rcubed
    #                 Nint_leaf += 1
    #             continue
            
    #         # otherwise, we are at a node, so we need to decide whether or not to open it
    #         # using geometric opening criterion
    #         Lnode = 1./(2**tr.level)
    #         d = np.linalg.norm(tr.com - pos0)
    #         # print(d, Lnode, theta)
    #         if d < Lnode/theta:
    #             current_level += 1
    #         else:
    #             # we are not opening this node, just use its com to compute acc
    #             acc += G * tr.mass * (tr.com - pos0) / (np.linalg.norm(tr.com - pos0)**3)
    #             Nint_node += 1
        
    #     return acc, Nint_leaf, Nint_node

class DirectGravity:
    def __init__(self, point_cloud):
        """
        Direct N^2 gravity computation class
        
        Args:
            point_cloud: PointCloud object containing particle data
        """
        self.point_cloud = point_cloud
    
    def compute_acceleration(self, pos_eval, G):
        """
        Compute acceleration at pos_eval using direct N^2 summation
        
        Args:
            pos_eval: position where to evaluate acceleration
            G: gravitational constant
            
        Returns:
            acc: computed acceleration
            Nint: number of interactions (always N for direct method)
        """
        pos_eval = np.asarray(pos_eval, dtype=np.float32)
        
        # Compute squared distances from all particles to evaluation point
        rsq = np.sum((self.point_cloud.pos_float - pos_eval)**2, axis=1)
        rcubed = rsq**(3./2.)
        
        # Compute acceleration contributions from all particles
        dacc = np.zeros_like(self.point_cloud.pos_float)
        np.divide(self.point_cloud.pos_float - pos_eval, rcubed[:, None], 
                 out=dacc, where=rcubed[:, None]!=0)
        
        # Sum all contributions
        acc = G * np.sum(dacc * self.point_cloud.mass[:, None], axis=0)
        
        # Number of interactions is always N for direct method
        Nint = len(self.point_cloud.mass)
        
        return acc, Nint
    
    def compute_accelerations_batch(self, pos_eval_list, G):
        """
        Compute accelerations for multiple evaluation points
        
        Args:
            pos_eval_list: array of positions where to evaluate acceleration
            G: gravitational constant
            
        Returns:
            acc_list: array of computed accelerations
            Nint_total: total number of interactions
        """
        pos_eval_list = np.asarray(pos_eval_list, dtype=np.float32)
        acc_list = np.zeros_like(pos_eval_list)
        
        for i, pos_eval in enumerate(pos_eval_list):
            acc_list[i], _ = self.compute_acceleration(pos_eval, G)
        
        Nint_total = len(pos_eval_list) * len(self.point_cloud.mass)
        
        return acc_list, Nint_total
    
if __name__ == '__main__':
    np.random.seed(9996)

    Npart = 10000

    # Maximum value for 32-bit unsigned integer
    MAX_INT32 = np.iinfo(np.uint32).max

    # Generate random 3D positions for particle field
    pos = np.random.randint(0, MAX_INT32, size=(Npart, 3), dtype=np.uint32)
    pos_float = ((pos.astype(np.float64)) / MAX_INT32).astype(np.float32)
    M = np.full(Npart, 1. / Npart)
    G = 1.

    # Softening length equal to interparticle spacing
    eps = np.full(Npart, np.sqrt(3.) / Npart)

    pcloud = PointCloud(pos, M, eps)
    bhtree = BHTree(pcloud)