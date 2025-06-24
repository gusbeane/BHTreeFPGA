#include <iostream>
#include "../include/bhtree_types.h"
#include "../include/bhtree_config.h"
#include "../include/pointcloud.h"
#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip>

std::string format_ph_key(uint32_t key, int level) {
    std::ostringstream oss;
    
    // Extract 3 bits per level, starting from the most significant bits
    for (int i = level - 1; i >= 0; i--) {
        if (i < level - 1) {
            oss << " ";  // Add space between levels
        }
        
        // Extract 3 bits for this level
        uint32_t level_bits = (key >> (3 * i)) & 0x7;  // 0x7 = 0b111
        
        // Convert to 3-bit binary string
        oss << ((level_bits & 4) ? '1' : '0')
            << ((level_bits & 2) ? '1' : '0') 
            << ((level_bits & 1) ? '1' : '0');
    }
    
    return oss.str();
}

std::string format_ph_key_padded(uint32_t key, int level, int pad_to_level = 5) {
    std::ostringstream oss;
    
    // Add the actual key bits first
    for (int i = level - 1; i >= 0; i--) {
        if (i < level - 1) {
            oss << " ";  // Add space between levels
        }
        
        // Extract 3 bits for this level
        uint32_t level_bits = (key >> (3 * i)) & 0x7;  // 0x7 = 0b111
        
        // Convert to 3-bit binary string
        oss << ((level_bits & 4) ? '1' : '0')
            << ((level_bits & 2) ? '1' : '0') 
            << ((level_bits & 1) ? '1' : '0');
    }
    
    // Add trailing padding spaces
    for (int i = level; i < pad_to_level; i++) {
        oss << "    ";  // 4 spaces (3 for digits + 1 for separator)
    }
    
    return oss.str();
}

void print_node(const NodeOrLeaf& node) {
    std::cout << "Level: " << node.level 
              << " Key: " << format_ph_key_padded(node.key, node.level) 
              << " Npart: " << std::left << std::setw(3) << node.Nparticles 
              << " StartIdx: " << node.start_idx
              << " COM: (" << std::fixed << std::setprecision(3) 
              << node.com_x << ", " << node.com_y << ", " << node.com_z << ")"
              << " Mass: " << node.mass << std::endl;
}

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

std::vector<NodeOrLeaf> add_particle_to_tree(Tree& temp_nodes, RealPosition3D pos, double mass, int key, int index) {
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

    std::reverse(new_nodes.begin(), new_nodes.end());

    return new_nodes;
}

std::vector<NodeOrLeaf> flush_out_nodes(Tree& temp_nodes) {
    std::vector<NodeOrLeaf> new_nodes;
    new_nodes.reserve(MAX_DEPTH);
    for (int i = 1; i < MAX_DEPTH+1; i++) {
        new_nodes.push_back(temp_nodes[i-1]);

        if(temp_nodes[i-1].Nparticles <= NLEAF) {
            break;
        }
    }
    std::reverse(new_nodes.begin(), new_nodes.end());
    return new_nodes;
}

Tree build_tree(PointCloud pc) {
    // tree is just a vector of NodeOrLeaf
    Tree tree;
    tree.reserve(static_cast<size_t>(1.5 * pc.size()));

    // Get the real positions (not integer ones)
    std::vector<RealPosition3D> pos = pc.get_real_positions();
    std::vector<double> mass = pc.get_masses();

    // We start with an array of temporary nodes with size of the max depth
    // We don't make the root node
    std::vector<NodeOrLeaf> temp_nodes(MAX_DEPTH);
    uint32_t key0 = pc.get_ph_keys()[0];
    for (int i = 1; i < MAX_DEPTH+1; i++) {
        temp_nodes[i-1].level = i;
        temp_nodes[i-1].key = key0 >> (3 * (MAX_DEPTH - i));
        add_particle_to_node(temp_nodes[i-1], pos[0], mass[0]);
    }

    // Now we loop through all the particles and add them to the tree
    for (int i = 1; i < pc.size(); i++) {
        std::vector<NodeOrLeaf> new_nodes = add_particle_to_tree(temp_nodes, pos[i], mass[i], pc.get_ph_keys()[i], i);
        tree.insert(tree.end(), new_nodes.begin(), new_nodes.end());
    }

    // Flush out the remaining nodes
    std::vector<NodeOrLeaf> new_nodes = flush_out_nodes(temp_nodes);
    tree.insert(tree.end(), new_nodes.begin(), new_nodes.end());

    // Reverse the order of the tree
    // Commented out for now for HLS
    // std::reverse(tree.begin(), tree.end());

    // Divide com by mass
    for (int i = 0; i < tree.size(); i++) {
        tree[i].com_x /= tree[i].mass;
        tree[i].com_y /= tree[i].mass;
        tree[i].com_z /= tree[i].mass;
    }

    return tree;
}

int main(int argc, char* argv[]) {
    PointCloud pc = PointCloud::random(1000);
    pc.sort_by_ph_keys();

    Tree tree = build_tree(pc);
    
    // Print first 10 nodes
    std::cout << "First 10 nodes:" << std::endl;
    for (int i = 0; i < 10; i++) {
        print_node(tree[tree.size() - i - 1]);
    }
    
    return 0;
} 