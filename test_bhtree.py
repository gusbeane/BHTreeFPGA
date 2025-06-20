import numpy as np
from tqdm import tqdm
from prototype import PointCloud, DirectGravity, BHTree

# Set random seed for reproducibility
np.random.seed(9996)

# Parameters (same as gen_golden.py)
Npart = 10000
Neval = 1000  # Fewer evaluation points for faster testing

# Maximum value for 32-bit unsigned integer
MAX_INT32 = np.iinfo(np.uint32).max

print(f"Generating {Npart} particles and {Neval} evaluation points...")

# Generate random 3D positions for particle field
pos = np.random.randint(0, MAX_INT32, size=(Npart, 3), dtype=np.uint32)
M = np.full(Npart, 1. / Npart)
G = 1.

# Softening length equal to interparticle spacing
eps = np.full(Npart, np.sqrt(3.) / Npart)

# Generate random 3D positions to evaluate acc at
pos_eval = np.random.randint(0, MAX_INT32, size=(Neval, 3), dtype=np.uint32)
pos_eval_float = ((pos_eval.astype(np.float64)) / MAX_INT32).astype(np.float32)

print("Creating PointCloud...")
# Create PointCloud object (handles PH key generation and sorting)
point_cloud = PointCloud(pos, M, eps)

print("Setting up gravity computation methods...")
# Create DirectGravity and BHTree instances
direct_gravity = DirectGravity(point_cloud)
bh_tree = BHTree(point_cloud)

print("Computing accelerations with DirectGravity (exact N^2 method)...")
# Compute accelerations using direct method
acc_direct = np.zeros((Neval, 3), dtype=np.float32)
total_interactions_direct = 0

for i in tqdm(range(Neval), desc="Direct computation"):
    acc_direct[i], n_int = direct_gravity.compute_acceleration(pos_eval_float[i], G)
    total_interactions_direct += n_int

print("Computing accelerations with BHTree (Barnes-Hut approximation)...")
# Compute accelerations using Barnes-Hut tree
acc_bh = np.zeros((Neval, 3), dtype=np.float32)
total_interactions_leaf = 0
total_interactions_node = 0

theta_values = [0.3, 0.5, 0.7]  # Different opening criteria
max_rel_err = [8e-3, 5e-2, 8e-2]

for idx, theta in enumerate(theta_values):
    print(f"\n--- Testing with theta = {theta} ---")
    
    acc_bh = np.zeros((Neval, 3), dtype=np.float32)
    total_interactions_leaf = 0
    total_interactions_node = 0
    
    for i in tqdm(range(Neval), desc=f"BH tree (θ={theta})"):
        acc_bh[i], n_leaf, n_node = bh_tree.walk_tree(pos_eval_float[i], G, theta=theta)
        total_interactions_leaf += n_leaf
        total_interactions_node += n_node
    
    # Compute relative differences
    acc_diff = acc_bh - acc_direct
    acc_magnitude = np.linalg.norm(acc_direct, axis=1)
    
    # Avoid division by zero
    mask = acc_magnitude > 1e-12
    relative_error = np.zeros(Neval)
    relative_error[mask] = np.linalg.norm(acc_diff[mask], axis=1) / acc_magnitude[mask]
    
    # Compute statistics
    max_relative_error = np.max(relative_error)
    mean_relative_error = np.mean(relative_error[mask])
    median_relative_error = np.median(relative_error[mask])
    
    total_interactions_bh = total_interactions_leaf + total_interactions_node
    speedup_factor = total_interactions_direct / total_interactions_bh
    
    # Print results
    print(f"\n=== RESULTS for θ = {theta} ===")
    print(f"Total particles: {Npart}")
    print(f"Evaluation points: {Neval}")
    print(f"Direct method interactions: {total_interactions_direct:,}")
    print(f"BH tree interactions: {total_interactions_bh:,} (leaf: {total_interactions_leaf:,}, node: {total_interactions_node:,})")
    print(f"Speedup factor: {speedup_factor:.1f}x")
    print(f"\nAccuracy:")
    print(f"Maximum relative error: {max_relative_error:.2e}")
    print(f"Mean relative error: {mean_relative_error:.2e}")
    print(f"Median relative error: {median_relative_error:.2e}")
    
    # Find worst case for detailed analysis
    worst_idx = np.argmax(relative_error)
    print(f"\nWorst case (evaluation point {worst_idx}):")
    print(f"  Direct acceleration: {acc_direct[worst_idx]}")
    print(f"  BH tree acceleration: {acc_bh[worst_idx]}")
    print(f"  Relative error: {relative_error[worst_idx]:.2e}")
    
    # Assert that the maximum relative error is within expected bounds
    assert max_relative_error < max_rel_err[idx], f"Maximum relative error {max_relative_error:.2e} exceeds threshold {max_rel_err[idx]:.2e} for θ={theta}"

print("\n=== TEST COMPLETE ===")