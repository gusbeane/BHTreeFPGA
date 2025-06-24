#include "../include/pointcloud.h"
#include <algorithm>
#include <random>
#include <climits>

// Constructor for random point cloud (default box_size = 1.0)
PointCloud::PointCloud(size_t num_points) : box_size(1.0) {
    generate_random_positions(num_points);
    update_real_positions();
}

// Constructor for random point cloud with specified box_size
PointCloud::PointCloud(size_t num_points, double box_size) : box_size(box_size) {
    generate_random_positions(num_points);
    update_real_positions();
}

// Constructor for random point cloud with specified box_size and masses
PointCloud::PointCloud(size_t num_points, double box_size, const std::vector<double>& masses) : box_size(box_size), masses(masses) {
    generate_random_positions(num_points);
    update_real_positions();
}

// Add a single position
void PointCloud::add_position(const Position3D& pos) {
    positions.push_back(pos);
    // Convert to real position and add to real_positions
    RealPosition3D real_pos;
    for (int i = 0; i < 3; i++) {
        real_pos[i] = (static_cast<double>(pos[i]) / UINT32_MAX) * box_size;
    }
    real_positions.push_back(real_pos);
}

// Add multiple positions
void PointCloud::add_positions(const PositionVector& pos_vec) {
    positions.insert(positions.end(), pos_vec.begin(), pos_vec.end());
    update_real_positions();
}

// Set positions (replaces existing ones)
void PointCloud::set_positions(const PositionVector& pos_vec) {
    positions = pos_vec;
    update_real_positions();
}

// Set box size and update real positions
void PointCloud::set_box_size(double new_box_size) {
    box_size = new_box_size;
    update_real_positions();
}

// Get box size
double PointCloud::get_box_size() const {
    return box_size;
}

// Convert integer positions to real positions
void PointCloud::update_real_positions() {
    real_positions.clear();
    real_positions.reserve(positions.size());
    
    for (const auto& pos : positions) {
        RealPosition3D real_pos;
        for (int i = 0; i < 3; i++) {
            real_pos[i] = (static_cast<double>(pos[i]) / UINT32_MAX) * box_size;
        }
        real_positions.push_back(real_pos);
    }
}

// Generate PH keys and sort positions by keys
void PointCloud::sort_by_ph_keys() {
    if (positions.empty()) {
        return;
    }
    
    // Generate Peano-Hilbert keys for all positions
    ph_keys = ph.generate_keys(positions);
    
    // Create pairs of keys and indices for sorting
    std::vector<std::pair<uint32_t, size_t> > key_index_pairs;
    key_index_pairs.reserve(ph_keys.size());
    
    for (size_t i = 0; i < ph_keys.size(); i++) {
        key_index_pairs.push_back(std::make_pair(ph_keys[i], i));
    }
    
    // Sort pairs by key values
    std::sort(key_index_pairs.begin(), key_index_pairs.end());
    
    // Create sorted vectors using the sorted indices
    std::vector<uint32_t> sorted_keys;
    PositionVector sorted_positions;
    std::vector<RealPosition3D> sorted_real_positions;
    
    sorted_keys.reserve(ph_keys.size());
    sorted_positions.reserve(positions.size());
    sorted_real_positions.reserve(real_positions.size());
    
    for (const auto& pair : key_index_pairs) {
        sorted_keys.push_back(pair.first);
        sorted_positions.push_back(positions[pair.second]);
        sorted_real_positions.push_back(real_positions[pair.second]);
    }
    
    // Replace original vectors with sorted ones
    ph_keys = std::move(sorted_keys);
    positions = std::move(sorted_positions);
    real_positions = std::move(sorted_real_positions);
}

// Get the sorted PH keys
const std::vector<uint32_t>& PointCloud::get_ph_keys() const {
    return ph_keys;
}

// Get the sorted integer positions
const PositionVector& PointCloud::get_positions() const {
    return positions;
}

// Get the sorted real positions
const std::vector<RealPosition3D>& PointCloud::get_real_positions() const {
    return real_positions;
}

// Get the masses
const std::vector<double>& PointCloud::get_masses() const {
    return masses;
}

// Get both sorted keys and integer positions
std::pair<std::vector<uint32_t>, PositionVector> PointCloud::get_sorted_data() const {
    return std::make_pair(ph_keys, positions);
}

// Get both sorted keys and real positions
std::pair<std::vector<uint32_t>, std::vector<RealPosition3D> > PointCloud::get_sorted_real_data() const {
    return std::make_pair(ph_keys, real_positions);
}

// Clear all data
void PointCloud::clear() {
    positions.clear();
    real_positions.clear();
    ph_keys.clear();
}

// Get the number of points
size_t PointCloud::size() const {
    return positions.size();
}

// Check if empty
bool PointCloud::empty() const {
    return positions.empty();
}

// Static factory method for creating random point clouds (box_size = 1.0)
PointCloud PointCloud::random(size_t num_points) {
    double mass = 1.0 / static_cast<double>(num_points);
    std::vector<double> masses(num_points, mass);
    return PointCloud(num_points, 1.0, masses);
}

// Generate random positions (0 to UINT32_MAX for each coordinate)
void PointCloud::generate_random_positions(size_t num_points) {
    positions.clear();
    positions.reserve(num_points);
    
    // Create random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> dis(0, UINT32_MAX);
    
    for (size_t i = 0; i < num_points; i++) {
        Position3D pos = {dis(gen), dis(gen), dis(gen)};
        positions.push_back(pos);
    }
} 