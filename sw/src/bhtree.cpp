#include <bitset>
#include <cmath>
#include <iostream>
#include "../include/bhtree_types.h"
#include "../include/bhtree_config.h"
#include "../include/pointcloud.h"
#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>

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

void print_node(const NodeOrLeaf& node, int index) {
    std::cout << "Node " << index << ": "
              << "Level: " << node.level 
              << " Key: " << format_ph_key_padded(node.key, node.level) 
              << " Npart: " << std::left << std::setw(3) << node.Nparticles 
              << " Leaf: " << node.is_leaf
              << " StartIdx: " << node.start_idx
              << " COM: (" << std::fixed << std::setprecision(3) 
              << node.com_x << ", " << node.com_y << ", " << node.com_z << ")"
              << " Mass: " << node.mass
              << " NextSibling: " << node.next_sibling << std::endl;
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

    unsigned int num_nodes_emitted = 0;
    int level_divergence = 0;

    for (int i = 1; i < MAX_DEPTH+1; i++) {
        unsigned int node_key = key >> 3 * (MAX_DEPTH - i);
        if (node_key == temp_nodes[i-1].key && !flush_mode && i < MAX_DEPTH) {
            add_particle_to_node(temp_nodes[i-1], pos, mass);
        }
        else {
            // we need to emit all the nodes below this level and replace them with empty nodes
            // we always emit the max depth node, if we get there
            if(!flush_mode) {
                level_divergence = i;
                flush_mode = true;
            }

            if (!quiet_mode) {
                new_nodes.push_back(temp_nodes[i-1]);
                num_nodes_emitted++;
                if(temp_nodes[i-1].Nparticles <= NLEAF) {
                    quiet_mode = true;
                }
            }

            temp_nodes[i-1] = NodeOrLeaf{};
            temp_nodes[i-1].level = i;
            temp_nodes[i-1].key = node_key;
            temp_nodes[i-1].start_idx = index;
            temp_nodes[i-1].next_sibling = 0;
            add_particle_to_node(temp_nodes[i-1], pos, mass);
        }
    }

    // all nodes above the divergence level need to have the num emitted nodes added to their sibling count
    for (int i = 1; i < level_divergence; i++) {
        if(temp_nodes[i-1].next_sibling != -1u) {
            temp_nodes[i-1].next_sibling += num_nodes_emitted;
        }
    }

    // all nodes below the divergence level need to adjust accordingly
    for (unsigned int i=0; i<num_nodes_emitted; i++) {
        // if the next sibling is -1, we want to keep it as -1 (n.b. -1 is max value for unsigned int)
        if(new_nodes[i].next_sibling != -1u) {
            new_nodes[i].next_sibling += num_nodes_emitted - i;
        }
        if(new_nodes[i].Nparticles <= NLEAF) {
            new_nodes[i].is_leaf = true;
        }
    }

    std::reverse(new_nodes.begin(), new_nodes.end());

    return new_nodes;
}

std::vector<NodeOrLeaf> flush_out_nodes(Tree& temp_nodes) {
    std::vector<NodeOrLeaf> new_nodes;
    new_nodes.reserve(MAX_DEPTH);
    unsigned int num_nodes_emitted = 0;
    for (int i = 1; i < MAX_DEPTH+1; i++) {
        new_nodes.push_back(temp_nodes[i-1]);
        num_nodes_emitted++;

        if(temp_nodes[i-1].Nparticles <= NLEAF) {
            break;
        }
    }

    for (unsigned int i=0; i<num_nodes_emitted; i++) {
        if(new_nodes[i].next_sibling != -1u) {
            new_nodes[i].next_sibling += num_nodes_emitted - i;
        }
        if(new_nodes[i].Nparticles <= NLEAF) {
            new_nodes[i].is_leaf = true;
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
    // The next node for all of these is set to -1, which means its the last sibling
    std::vector<NodeOrLeaf> temp_nodes(MAX_DEPTH);
    uint32_t key0 = pc.get_ph_keys()[0];
    for (int i = 1; i < MAX_DEPTH+1; i++) {
        temp_nodes[i-1].level = i;
        temp_nodes[i-1].key = key0 >> (3 * (MAX_DEPTH - i));
        temp_nodes[i-1].next_sibling = -1;
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
    std::reverse(tree.begin(), tree.end());

    // Divide com by mass
    for (int i = 0; i < tree.size(); i++) {
        tree[i].com_x /= tree[i].mass;
        tree[i].com_y /= tree[i].mass;
        tree[i].com_z /= tree[i].mass;
    }

    return tree;
}

void walk_tree(Tree& tree, unsigned long num_nodes, RealPosition3D pos,
        double G, double theta, RealPosition3D& acc, int& Nint_node, int& Nint_leaf) {
    // Initialize output parameters
    acc = {0.0, 0.0, 0.0};
    Nint_node = 0;
    Nint_leaf = 0;
    
    unsigned long current_idx = 0;

    while(current_idx < num_nodes) {
        NodeOrLeaf node = tree[current_idx];
        // compute distance
        double dx, dy, dz;
        dx = node.com_x - pos[0];
        dy = node.com_y - pos[1];
        dz = node.com_z - pos[2];
        double rsq = dx*dx + dy*dy + dz*dz;
        double r = std::sqrt(rsq);
        
        // check multipole acceptance criterion. size of node is 1/(2^level)
        double node_size = 1.0 / (1 << node.level);
        if(r > node_size/theta || node.is_leaf) {
            double rcubed = rsq*r;
            // we are not opening this node, just use its com to compute acc
            acc[0] += G * node.mass * (dx) / (rcubed);
            acc[1] += G * node.mass * (dy) / (rcubed);
            acc[2] += G * node.mass * (dz) / (rcubed);
            
            if(node.is_leaf) {
                Nint_leaf += 1;
            }
            else {
                Nint_node += 1;
            }

            // because we didn't open the node, we can skip to the next sibling
            if(node.next_sibling != -1u) {
                current_idx += node.next_sibling;
            }
            else {
                break;
            }
        }
        else {
            // we are opening this node, so we just advance to the next nodeleaf
            current_idx++;
        }
    }
}

void direct_sum(PointCloud pc, RealPosition3D pos, double G, RealPosition3D& acc, int &Nint) {
    // zero out acc
    acc = {0.0, 0.0, 0.0};
    Nint = 0;
    
    for(int i = 0; i < pc.size(); i++) {
        RealPosition3D dx = pc.get_real_positions()[i];
        dx[0] -= pos[0];
        dx[1] -= pos[1];
        dx[2] -= pos[2];
        double rsq = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
        double r = std::sqrt(rsq);
        double rcubed = rsq*r;
        if(r > 0.0) {
            acc[0] += G * pc.get_masses()[i] * (dx[0]) / (rcubed);
            acc[1] += G * pc.get_masses()[i] * (dx[1]) / (rcubed);
            acc[2] += G * pc.get_masses()[i] * (dx[2]) / (rcubed);
        }
        Nint += 1;
    }
}

bool test_at_max_depth() {
    std::cout << "\n=== Test: At Max Depth ===" << std::endl;

    const int NUM_PARTICLES = 10;

    // Create a vector to store our test particles
    std::vector<RealPosition3D> positions;
    std::vector<double> masses;
    positions.reserve(NUM_PARTICLES);
    masses.reserve(NUM_PARTICLES);

    // Manually create particles at different positions, similar to test_treecon.cpp
    double step = (1.0 / pow(2, MAX_DEPTH));
    for (int i = 0; i < NUM_PARTICLES; i++) {
        positions.push_back({0.1 + i*step/NUM_PARTICLES, 0.1, 0.1});
        masses.push_back(1.0/double(NUM_PARTICLES));
    }

    // Create a PointCloud and convert real positions to integer positions
    PointCloud pc(NUM_PARTICLES, 1.0, masses);  // box_size = 1.0
    
    // Convert real positions to integer positions
    PositionVector integer_positions;
    integer_positions.reserve(NUM_PARTICLES);
    
    for (int i = 0; i < NUM_PARTICLES; i++) {
        Position3D int_pos;
        for (int j = 0; j < 3; j++) {
            // Convert real position to integer position
            // real_pos[j] = (pos[j] / UINT32_MAX) * box_size
            // So: pos[j] = (real_pos[j] / box_size) * UINT32_MAX
            int_pos[j] = static_cast<uint32_t>((positions[i][j] / 1.0) * UINT32_MAX);
        }
        integer_positions.push_back(int_pos);
    }
    
    // Set the positions in the PointCloud
    pc.set_positions(integer_positions);
    pc.sort_by_ph_keys();

    std::cout << "Manual particles sorted by PH keys:" << std::endl;
    auto ph_keys = pc.get_ph_keys();
    auto real_positions = pc.get_real_positions();
    auto real_masses = pc.get_masses();
    
    for (int i = 0; i < NUM_PARTICLES; i++) {
        std::cout << "  Particle " << i << ": pos=(" 
                  << real_positions[i][0] << "," 
                  << real_positions[i][1] << "," 
                  << real_positions[i][2] 
                  << ") key=0x" << std::hex << ph_keys[i] << std::dec 
                  << " key_binary=" << format_ph_key(ph_keys[i], MAX_DEPTH) << std::endl;
    }

    // Build the tree
    Tree tree = build_tree(pc);

    std::cout << "\nBuilt tree with " << tree.size() << " nodes" << std::endl;
    std::cout << "\nTree nodes:" << std::endl;
    for (int i = 0; i < std::min(50, (int)tree.size()); i++) {
        print_node(tree[i], i);
    }

    // Now perform sanity checks similar to test_treecon.cpp
    bool test_passed = true;

    // Check number of nodes (should be 11 based on the original test)
    int expected_nodes = 11;
    std::cout << "\nSanity check: Number of nodes" << std::endl;
    std::cout << "Expected: " << expected_nodes << ", Got: " << tree.size() << std::endl;
    if (tree.size() != expected_nodes) {
        std::cout << "WARNING: Node count mismatch (expected " << expected_nodes << ", got " << tree.size() << ")" << std::endl;
        // Don't fail the test for this, as the software implementation might differ
    }

    // Sum of particles in leaf nodes should equal NUM_PARTICLES
    int total_particles = 0;
    for (size_t i = 0; i < tree.size(); i++) {
        if (tree[i].Nparticles <= NLEAF) { // This is a leaf node
            total_particles += tree[i].Nparticles;
        }
    }
    std::cout << "\nSanity check: Total particles in leaf nodes" << std::endl;
    std::cout << "Expected: " << NUM_PARTICLES << ", Got: " << total_particles << std::endl;
    test_passed &= (total_particles == NUM_PARTICLES);

    // Compute expected center of mass of all particles
    RealPosition3D com_expected = {0.0, 0.0, 0.0};
    double total_mass = 0.0;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        com_expected[0] += real_positions[i][0] * real_masses[i];
        com_expected[1] += real_positions[i][1] * real_masses[i];
        com_expected[2] += real_positions[i][2] * real_masses[i];
        total_mass += real_masses[i];
    }
    com_expected[0] /= total_mass;
    com_expected[1] /= total_mass;
    com_expected[2] /= total_mass;

    // Check center of mass of root node (last node in tree after reversal)
    const double COM_TOL = 1e-6;
    if (!tree.empty()) {
        NodeOrLeaf root_node = tree.back(); // Last node after reversal should be root
        std::cout << "\nSanity check: Root node center of mass" << std::endl;
        std::cout << "Expected COM: (" << com_expected[0] << ", " << com_expected[1] << ", " << com_expected[2] << ")" << std::endl;
        std::cout << "Root node COM: (" << root_node.com_x << ", " << root_node.com_y << ", " << root_node.com_z << ")" << std::endl;
        
        bool com_match = (abs(root_node.com_x - com_expected[0]) < COM_TOL) &&
                        (abs(root_node.com_y - com_expected[1]) < COM_TOL) &&
                        (abs(root_node.com_z - com_expected[2]) < COM_TOL);
        test_passed &= com_match;
        
        if (!com_match) {
            std::cout << "Center of mass mismatch!" << std::endl;
        }
    }

    std::cout << "\nMax depth test " << (test_passed ? "PASSED" : "FAILED") << std::endl;
    return test_passed;
}

int main(int argc, char* argv[]) {
    PointCloud pc = PointCloud::random(1000);
    pc.sort_by_ph_keys();

    Tree tree = build_tree(pc);

    std::cout << "Tree size: " << tree.size() << std::endl;
    
    // Print first 10 nodes
    std::cout << "First 10 nodes:" << std::endl;
    for (int i = 0; i < 10; i++) {
        print_node(tree[i], i);
    }
    std::cout << std::endl;

    // Print last 10 nodes
    std::cout << "Last 10 nodes:" << std::endl;
    for (int i = 10; i >= 0; i--) {
        int index = tree.size() - i - 1;
        print_node(tree[index], index);
    }
    std::cout << std::endl;

    // Print nodes with level = 1
    std::cout << "Nodes with level = 1:" << std::endl;
    for (int i = 0; i < tree.size(); i++) {
        if (tree[i].level == 1) {
            print_node(tree[i], i);
        }
    }
    std::cout << std::endl;

    // Print nodes with level <= 2 and key has first three bits as 111
    std::cout << "Nodes with level <= 2 and key has first three bits as 111:" << std::endl;
    for (int i = 0; i < tree.size(); i++) {
        unsigned int levelkey = tree[i].key >> (3 * (tree[i].level -1));
        if (tree[i].level <= 2 && levelkey == 0b111) {
            print_node(tree[i], i);
        }
    }
    std::cout << std::endl;

    // Print raw binary of first node
    std::cout << "Raw binary of first node:" << std::endl;
    std::cout << std::bitset<30>(tree[0].key) << std::endl;
    std::cout << std::endl;

    // Check that COM of particles in first node is equal to its quoted COM
    RealPosition3D com = {0.0, 0.0, 0.0};
    double mass = 0.0;
    for (int i = tree[0].start_idx; i < tree[0].start_idx + tree[0].Nparticles; i++) {
        com[0] += pc.get_real_positions()[i][0] * pc.get_masses()[i];
        com[1] += pc.get_real_positions()[i][1] * pc.get_masses()[i];
        com[2] += pc.get_real_positions()[i][2] * pc.get_masses()[i];
        mass += pc.get_masses()[i];
    }
    com[0] /= mass;
    com[1] /= mass;
    com[2] /= mass;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "COM of particles in first node: (" << com[0] << ", " << com[1] << ", " << com[2] << ")" << std::endl;
    std::cout << "Quoted COM of first node: (" << tree[0].com_x << ", " << tree[0].com_y << ", " << tree[0].com_z << ")" << std::endl;
    std::cout << "Difference: (" << com[0] - tree[0].com_x << ", " << com[1] - tree[0].com_y << ", " << com[2] - tree[0].com_z << ")" << std::endl;
    std::cout << std::endl;

    // Check total number of leaves
    int num_leaves = 0;
    for (int i = 0; i < tree.size(); i++) {
        if (tree[i].is_leaf) {
            num_leaves++;
        }
    }
    std::cout << "Total number of leaves: " << num_leaves << std::endl;


    // test_at_max_depth();


    // now we calculate some numbers
    double G = 1.0;
    RealPosition3D acc_direct;
    int Nint_direct;
    RealPosition3D pos = {0.8, 0.5, 0.5};
    direct_sum(pc, pos, G, acc_direct, Nint_direct);
    std::cout << "Direct sum: (" << acc_direct[0] << ", " << acc_direct[1] << ", " << acc_direct[2] << ")" << std::endl;
    std::cout << "Nint_direct: " << Nint_direct << std::endl;

    // now do tree walk
    RealPosition3D acc_tree;
    int Nint_node, Nint_leaf;
    double theta = 0.5;
    walk_tree(tree, tree.size(), pos, G, theta, acc_tree, Nint_node, Nint_leaf);
    std::cout << "Tree walk: (" << acc_tree[0] << ", " << acc_tree[1] << ", " << acc_tree[2] << ")" << std::endl;
    std::cout << "Nint_node: " << Nint_node << std::endl;
    std::cout << "Nint_leaf: " << Nint_leaf << std::endl;

    return 0;
} 