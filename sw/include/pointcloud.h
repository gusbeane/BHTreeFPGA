#ifndef POINTCLOUD_H
#define POINTCLOUD_H

#include <vector>
#include <utility>
#include "peano_hilbert.h"
#include "bhtree_types.h"

class PointCloud {
private:
    PositionVector positions;           // Integer positions
    std::vector<RealPosition3D> real_positions;  // Real (floating-point) positions
    std::vector<uint32_t> ph_keys;
    PeanoHilbert ph;
    double box_size;                    // Box size for coordinate conversion
    
    // Private helper methods
    void generate_random_positions(size_t num_points);
    void update_real_positions();       // Convert integer to real positions

public:
    // Constructors
    PointCloud() : box_size(1.0) {}
    PointCloud(size_t num_points);      // Uses default box_size = 1.0
    PointCloud(size_t num_points, double box_size);
    
    // Position management
    void add_position(const Position3D& pos);
    void add_positions(const PositionVector& pos_vec);
    void set_positions(const PositionVector& pos_vec);
    
    // Box size management
    void set_box_size(double new_box_size);
    double get_box_size() const;
    
    // Core functionality
    void sort_by_ph_keys();
    
    // Data access - integer positions
    const std::vector<uint32_t>& get_ph_keys() const;
    const PositionVector& get_positions() const;
    std::pair<std::vector<uint32_t>, PositionVector> get_sorted_data() const;
    
    // Data access - real positions  
    const std::vector<RealPosition3D>& get_real_positions() const;
    std::pair<std::vector<uint32_t>, std::vector<RealPosition3D> > get_sorted_real_data() const;
    
    // Utility methods
    void clear();
    size_t size() const;
    bool empty() const;
    
    // Static factory methods
    static PointCloud random(size_t num_points);                    // Uses box_size = 1.0
    static PointCloud random(size_t num_points, double box_size);
};

#endif // POINTCLOUD_H 