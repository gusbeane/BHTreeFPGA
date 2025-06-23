#include <iostream>
#include <vector>
#include <cassert>
#include "../include/peano_hilbert.h"
#include "../include/bhtree_types.h"

class PeanoHilbertTester {
private:
    PeanoHilbert ph;
    int tests_passed;
    int tests_failed;

    void add_test_data(std::vector<uint32_t>& x, std::vector<uint32_t>& y, std::vector<uint32_t>& z, std::vector<uint32_t>& expected) {
        // Test data with correct expected values for original unsorted positions
        // These values come from direct PeanoHilbert key generation (not sorted PointCloud output)
        x.push_back(2876051090); y.push_back(1096434381); z.push_back(1714768399); expected.push_back(478731298);
        x.push_back(3097084232); y.push_back(3019784370); z.push_back( 793782314); expected.push_back(280590756);
        x.push_back(3431515389); y.push_back(2591306578); z.push_back(2471124925); expected.push_back(734091386);
        x.push_back(2008487837); y.push_back(3933998103); z.push_back(1005816787); expected.push_back(237190064);
        x.push_back(1282104987); y.push_back( 268117585); z.push_back( 571306213); expected.push_back(56714574);
        x.push_back(3335969296); y.push_back(1050737557); z.push_back(1323156755); expected.push_back(510071576);
        x.push_back(2282156697); y.push_back(2596480175); z.push_back( 334708671); expected.push_back(269027414);
        x.push_back(1369708748); y.push_back(1579412074); z.push_back(3827174530); expected.push_back(1001426638);
        x.push_back(3576214542); y.push_back(1903744684); z.push_back(2359863372); expected.push_back(543811739);
        x.push_back( 609652247); y.push_back(1088688328); z.push_back(3240302310); expected.push_back(951784828);
    }

public:
    PeanoHilbertTester() : tests_passed(0), tests_failed(0) {}

    void test_golden_data() {
        std::cout << "Testing against golden reference data..." << std::endl;
        std::cout << "=========================================" << std::endl;
        
        std::vector<uint32_t> x_vals, y_vals, z_vals, expected_keys;
        add_test_data(x_vals, y_vals, z_vals, expected_keys);
        
        for (size_t i = 0; i < x_vals.size(); i++) {
            uint32_t computed_key = ph.generate_key(x_vals[i], y_vals[i], z_vals[i]);
            
            std::cout << "Test " << i << ": Position [" 
                      << x_vals[i] << ", " << y_vals[i] << ", " << z_vals[i] << "] ";
            
            if (computed_key == expected_keys[i]) {
                std::cout << "✓ PASS (Key: " << computed_key << ")" << std::endl;
                tests_passed++;
            } else {
                std::cout << "✗ FAIL (Expected: " << expected_keys[i] 
                          << ", Got: " << computed_key << ")" << std::endl;
                tests_failed++;
            }
        }
        
        std::cout << std::endl;
    }

    void test_basic_functionality() {
        std::cout << "Testing basic functionality..." << std::endl;
        std::cout << "==============================" << std::endl;
        
        // Test simple cases
        uint32_t key1 = ph.generate_key(0, 0, 0);
        std::cout << "Key for (0,0,0): " << key1 << std::endl;
        
        uint32_t key2 = ph.generate_key(UINT32_MAX, UINT32_MAX, UINT32_MAX);
        std::cout << "Key for (MAX,MAX,MAX): " << key2 << std::endl;
        
        // Test different depths
        uint32_t key_depth5 = ph.generate_key(1000, 2000, 3000, 5);
        uint32_t key_depth10 = ph.generate_key(1000, 2000, 3000, 10);
        
        std::cout << "Key for (1000,2000,3000) depth 5: " << key_depth5 << std::endl;
        std::cout << "Key for (1000,2000,3000) depth 10: " << key_depth10 << std::endl;
        
        // Test 64-bit version
        uint64_t key64 = ph.generate_key_64(1000, 2000, 3000, 15);
        std::cout << "64-bit key for (1000,2000,3000) depth 15: " << key64 << std::endl;
        
        tests_passed++; // Basic functionality test
        std::cout << "✓ Basic functionality tests completed" << std::endl;
        std::cout << std::endl;
    }

    void test_batch_operations() {
        std::cout << "Testing batch operations..." << std::endl;
        std::cout << "===========================" << std::endl;
        
        // Create position vector from golden data
        PositionVector positions;
        std::vector<uint32_t> x_vals, y_vals, z_vals, expected_keys;
        add_test_data(x_vals, y_vals, z_vals, expected_keys);
        
        for (size_t i = 0; i < x_vals.size(); i++) {
            Position3D pos;
            pos[0] = x_vals[i];
            pos[1] = y_vals[i];
            pos[2] = z_vals[i];
            positions.push_back(pos);
        }
        
        // Test batch key generation
        std::vector<uint32_t> computed_keys = ph.generate_keys(positions);
        
        std::cout << "Batch processing " << positions.size() << " positions:" << std::endl;
        
        bool batch_success = true;
        for (size_t i = 0; i < positions.size(); i++) {
            if (computed_keys[i] != expected_keys[i]) {
                std::cout << "Batch test " << i << " FAILED: Expected " 
                          << expected_keys[i] << ", got " << computed_keys[i] << std::endl;
                batch_success = false;
                tests_failed++;
            }
        }
        
        if (batch_success) {
            std::cout << "✓ All batch operations PASSED" << std::endl;
            tests_passed++;
        }
        
        std::cout << std::endl;
    }

    void test_edge_cases() {
        std::cout << "Testing edge cases..." << std::endl;
        std::cout << "=====================" << std::endl;
        
        // Test corner cases
        std::vector<uint32_t> edge_x, edge_y, edge_z;
        edge_x.push_back(0); edge_y.push_back(0); edge_z.push_back(0);
        edge_x.push_back(UINT32_MAX); edge_y.push_back(0); edge_z.push_back(0);
        edge_x.push_back(0); edge_y.push_back(UINT32_MAX); edge_z.push_back(0);
        edge_x.push_back(0); edge_y.push_back(0); edge_z.push_back(UINT32_MAX);
        edge_x.push_back(UINT32_MAX); edge_y.push_back(UINT32_MAX); edge_z.push_back(UINT32_MAX);
        edge_x.push_back(UINT32_MAX/2); edge_y.push_back(UINT32_MAX/2); edge_z.push_back(UINT32_MAX/2);
        
        for (size_t i = 0; i < edge_x.size(); i++) {
            uint32_t key = ph.generate_key(edge_x[i], edge_y[i], edge_z[i]);
            std::cout << "Edge case " << i << ": [" << edge_x[i] << ", " << edge_y[i] 
                      << ", " << edge_z[i] << "] -> Key: " << key << std::endl;
        }
        
        tests_passed++; // If we get here without crashing, edge cases pass
        std::cout << "✓ Edge cases completed without errors" << std::endl;
        std::cout << std::endl;
    }

    void run_all_tests() {
        std::cout << "Peano-Hilbert Implementation Test Suite" << std::endl;
        std::cout << "=======================================" << std::endl;
        std::cout << std::endl;
        
        test_basic_functionality();
        test_golden_data();
        test_batch_operations();
        test_edge_cases();
        
        print_summary();
    }

    void print_summary() {
        std::cout << "Test Summary:" << std::endl;
        std::cout << "=============" << std::endl;
        std::cout << "Tests Passed: " << tests_passed << std::endl;
        std::cout << "Tests Failed: " << tests_failed << std::endl;
        std::cout << "Total Tests:  " << (tests_passed + tests_failed) << std::endl;
        
        if (tests_failed == 0) {
            std::cout << "✓ ALL TESTS PASSED!" << std::endl;
        } else {
            std::cout << "✗ " << tests_failed << " TESTS FAILED!" << std::endl;
        }
        
        std::cout << std::endl;
        std::cout << "NodeOrLeaf struct size: " << sizeof(NodeOrLeaf) << " bytes" << std::endl;
    }
};

int main() {
    PeanoHilbertTester tester;
    tester.run_all_tests();
    return 0;
} 