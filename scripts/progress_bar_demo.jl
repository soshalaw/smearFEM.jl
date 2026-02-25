using ProgressMeter

println("=== Progress Bar Demo in Julia ===\n")

# Example 1: Thin, clean progress bar (default style)
println("Example 1: Thin Progress Bar")
println("-" ^ 40)
p = Progress(100, dt=0.1, desc="Processing: ", barlen=50, color=:green)
for i in 1:100
    sleep(0.05)  # Simulate work
    next!(p)
end
println()

# Example 2: Compact style (shorter bar)
println("\nExample 2: Compact Style")
println("-" ^ 40)
p = Progress(50, dt=0.1, desc="Computing: ", barlen=30, showspeed=false)
for i in 1:50
    sleep(0.05)
    next!(p)
end
println()

# Example 3: Progress bar with output below
println("\nExample 3: Progress with Output Below")
println("-" ^ 40)
p = Progress(50, dt=0.1, desc="Loading: ", barlen=50, output=stderr)
for i in 1:50
    sleep(0.05)
    # You can print to stdout while progress bar is on stderr
    if i % 10 == 0
        println("Checkpoint reached at iteration $i")
    end
    next!(p)
end
println()

# Example 4: Progress with percentage only (very clean)
println("\nExample 4: Minimal Style with Percentage")
println("-" ^ 40)
p = Progress(50, dt=0.1, desc="Computing: ", barlen=40)
for i in 1:50
    sleep(0.05)
    next!(p)
end
println()

# Example 5: Progress with status updates
println("\nExample 5: Progress with Status Updates")
println("-" ^ 40)
p = Progress(100, dt=0.1, desc="Status: ", barlen=50)
for i in 1:100
    sleep(0.03)
    
    # Update with status message
    if i == 25
        next!(p; showvalues = [(:Phase, "Initialization")])
    elseif i == 50
        next!(p; showvalues = [(:Phase, "Processing")])
    elseif i == 75
        next!(p; showvalues = [(:Phase, "Finalizing")])
    else
        next!(p)
    end
end
println()

println("\n✓ All examples completed!")
