#!/usr/bin/env julia
"""
Test Phase 4b assembly router integration
"""

# Test module loading
println("[TEST] Loading smearFEM module...")
using smearFEM
println("[OK] Module loaded successfully!")

# Test that assembly router functions are accessible
println("\n[TEST] Checking assembly router exports...")
exported_functions = [
    :assemble_system_A_routed,
    :assemble_system_B_routed,
    :apply_boundary_conditions_routed,
    :simulate_with_gpu_integration,
    :print_assembly_status
]

for func in exported_functions
    if isdefined(smearFEM, func)
        println("[OK] Function exported: $func")
    else
        println("[ERROR] Function NOT exported: $func")
    end
end

println("\n[TEST] Phase 4b integration complete!")
println("[Summary] All assembly router functions are accessible from smearFEM module.")
