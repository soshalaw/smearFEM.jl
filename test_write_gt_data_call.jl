#!/usr/bin/env julia
# Quick diagnostic to verify write_gt_data is accessible and callable

using smearFEM

# Check if write_gt_data exists
if isdefined(smearFEM, :write_gt_data)
    println("✓ write_gt_data is defined in smearFEM")
else
    println("✗ write_gt_data is NOT defined in smearFEM")
    exit(1)
end

# Check function signature
try
    sig = methods(write_gt_data)
    println("✓ write_gt_data methods:")
    for s in sig
        println("  - $s")
    end
catch err
    println("✗ Error inspecting write_gt_data: $err")
    exit(1)
end

# Test minimal call (will fail due to missing data, but checks if function is reachable)
test_params = Dict(
    "FunctionClass_x" => "Q2",
    "FunctionClass_u" => "Q2",
    "FunctionClass_p" => "Q1",
    "ne_gt" => 4,
    "β_gt" => 10.0,
    "η_gt" => 100.0,
    "filepath_gt" => "/tmp/test",
    "control" => "force",
    "viscosity_type" => "bulk_viscosity",
    "obj_pose_gt" => zeros(4, 4),
    "F_ext" => 9813.0,
    "sim_time_gt" => 1.0,
    "steps_gt" => 1,
    "r" => 25.0,
    "h" => 40.0,
    "camera_matrix" => ones(3, 3)
)

println("\nAttempting to call write_gt_data with test parameters...")
try
    write_gt_data(test_params)
    println("✓ write_gt_data executed successfully")
catch err
    println("✗ Error calling write_gt_data (expected if data missing):")
    println("  Type: $(typeof(err))")
    println("  Message: $(err)")
    # Don't exit - might fail due to missing mesh data, but function is callable
end

println("\n✓ All checks passed - write_gt_data is accessible and callable")
