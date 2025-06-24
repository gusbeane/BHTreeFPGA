#include <iostream>
#include "../include/peano_hilbert.h"
#include "../include/bhtree_types.h"
#include "../include/bhtree_config.h"
#include <random>
#include "../include/pointcloud.h"

void add_particle_to_node(NodeOrLeaf& node, RealPosition3D pos, double mass, int index = -1) {
    node.com_x += pos[0] * mass;
    node.com_y += pos[1] * mass;
    node.com_z += pos[2] * mass;
    node.mass += mass;
    node.Nparticles++;

    if (index != -1) {
        node.start_idx = index;
    }
}

std::vector<NodeOrLeaf> add_particle_to_tree(Tree& tree, Tree& temp_nodes, RealPosition3D pos, double mass, int key, int index) {
    // we loop through the levels and find the node that has the first discrepancy
    std::vector<NodeOrLeaf> new_nodes;
    new_nodes.reserve(MAX_DEPTH);
    bool flush_mode = false;
    bool quiet_mode = false;

    for (int i = 1; i < MAX_DEPTH+1; i++) {
        int node_key = key >> 3 * (MAX_DEPTH - i);
        if (node_key == temp_nodes[i-1].key && !flush_mode) {
            add_particle_to_node(temp_nodes[i-1], pos, mass);
        }
        else {
            // we need to emit all the nodes below this level and replace them with empty nodes
            flush_mode = true;

            if (!quiet_mode) {
                new_nodes.push_back(temp_nodes[i-1]);
            }

            temp_nodes[i-1] = NodeOrLeaf{};
            temp_nodes[i-1].level = i;
            temp_nodes[i-1].key = node_key;
            add_particle_to_node(temp_nodes[i-1], pos, mass);

            if (temp_nodes[i-1].Nparticles <= NLEAF) {
                quiet_mode = true;
            }
        }
    }

    return new_nodes;
}

Tree build_tree(PointCloud pc) {
    // tree is just a vector of NodeOrLeaf
    Tree tree;
    tree.reserve(static_cast<size_t>(1.5 * pc.size()));

    // Get the real positions (not integer ones)
    std::vector<RealPosition3D> pos = pc.get_real_positions();

    // We start with an array of temporary nodes with size of the max depth
    // We don't make the root node
    std::vector<NodeOrLeaf> temp_nodes(MAX_DEPTH);
    uint32_t key0 = pc.get_ph_keys()[0];
    for (int i = 1; i < MAX_DEPTH+1; i++) {
        temp_nodes[i-1].level = i;
        temp_nodes[i-1].key = key0 >> (3 * (MAX_DEPTH - i));
        add_particle_to_node(temp_nodes[i-1], pos[0], 1.0);
    }

    // Now we loop through all the particles and add them to the tree
    for (int i = 1; i < pc.size(); i++) {
        std::vector<NodeOrLeaf> new_nodes = add_particle_to_tree(tree, temp_nodes, pos[i], 1.0, pc.get_ph_keys()[i], i);
        tree.insert(tree.end(), new_nodes.begin(), new_nodes.end());
    }

    

    return tree;
}

int main(int argc, char* argv[]) {
    PointCloud pc = PointCloud::random(1000);
    pc.sort_by_ph_keys();
    
    return 0;
} 