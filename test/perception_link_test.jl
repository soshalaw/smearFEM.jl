using smearFEM
using HDF5
using Sockets
using Test

# Fixture written by smearPerception's own code (data/perception_link/make_fixture.py, which lists
# the record values checked here): its driver's HDF5 file, the bytes its publisher sent for the same
# records, and its sampler's values. Guards the byte layout, the axis reversal, the unit
# conversions and the sampling rule across languages.
const LINK_FIXTURE = joinpath(@__DIR__, "data", "perception_link")

# Probe pose of fixture step s [mm].
fixture_pose(s) = [cosd(30) -sind(30) 0 10; sind(30) cosd(30) 0 20; 0 0 1 30 + (s - 1); 0 0 0 1]

@testset "perception link" begin
    run = read_perception_steps(joinpath(LINK_FIXTURE, "steps.h5"))
    grid, steps = run.grid, run.steps
    points, samples = h5open(f -> (read(f["points"]), read(f["samples"])), joinpath(LINK_FIXTURE, "expected.h5"))

    @testset "offline file and sampling" begin
        @test grid.shape == [5, 7, 9]
        @test grid.x0 ≈ [12.0, -7.0, 0.0] && grid.delta ≈ 1.0 && grid.offset == [2, 3, 0]
        @test run.t0_ns == 300 && run.dt_ns == 100
        @test length(steps) == 3
        for (s, st) in enumerate(steps)
            @test st.n == 6 + s && st.t_ns == run.t0_ns + st.n * run.dt_ns && st.n_views == 5 - s
            @test size(st.phi) == (5, 7, 9) && size(st.count) == size(st.phi)
            @test size(st.points) == (3, 10 + 7 * (s - 1))
            @test all(st.count[isnan.(st.phi)] .== 0)
            v = sample_field(st.phi, grid, 1e3 .* points)
            @test isequal(isnan.(v), isnan.(samples[:, s]))
            ok = .!isnan.(v)
            @test count(ok) > 50    # not a vacuous comparison (97–151 of 302 points are defined)
            @test maximum(abs.(v[ok] .- 1e3 .* samples[ok, s])) < 1e-5
        end
        # N → kg·mm/s², N·m → kg·mm²/s², the pose back in Julia's orientation and in mm.
        @test steps[1].wrench ≈ [100, 0, 0, 1e4, 2e4, 3e4] && steps[1].contact == 0
        @test steps[2].wrench ≈ [1e3, 2e3, 2e3, 1e5, 2e5, 3e5] && steps[2].contact == 1
        @test all(isnan, steps[3].wrench) && steps[3].contact == -1
        @test steps[1].probe_pose ≈ fixture_pose(1) && steps[2].probe_pose ≈ fixture_pose(2)
        @test all(isnan, steps[3].probe_pose)
    end

    @testset "stream" begin
        link = PerceptionLink(IOBuffer(read(joinpath(LINK_FIXTURE, "stream.bin"))))
        @test link.grid.shape == grid.shape && link.grid.x0 ≈ grid.x0 && link.grid.offset == grid.offset
        @test link.t0_ns == run.t0_ns && link.dt_ns == run.dt_ns && link.meta["contact_threshold"] == 0.5
        for st in steps
            got = next_step!(link)
            @test (got.n, got.t_ns, got.n_views, got.contact) == (st.n, st.t_ns, st.n_views, st.contact)
            @test isequal(got.phi, st.phi) && got.count == st.count
            @test got.points == st.points && isequal(got.normals, st.normals)
            @test isequal(got.wrench, st.wrench) && isequal(got.probe_pose, st.probe_pose)
        end
        @test next_step!(link) === nothing
    end

    Sys.isunix() && @testset "socket" begin
        path = "/tmp/smear_test_$(getpid()).sock"
        bytes = read(joinpath(LINK_FIXTURE, "stream.bin"))
        # Sends `data` to the first client, as perception's publisher would.
        function serve(data)
            rm(path; force=true)
            server = listen(path)
            return @async try
                sock = accept(server)
                write(sock, data)
                close(sock)
            finally
                close(server)
                rm(path; force=true)
            end
        end
        try
            touch(path)    # a stale file, as a crashed run leaves: retried, not an error
            @test_throws ErrorException open_perception_link(path; timeout=0.2)
            task = serve(bytes)
            link = open_perception_link(path; timeout=5)
            @test [next_step!(link).n for _ in steps] == [7, 8, 9]
            @test next_step!(link) === nothing
            close_link(link)
            wait(task)
            # Perception gone without END (a crash), here in the middle of the last step.
            task = serve(bytes[1:end-100])
            link = open_perception_link(path; timeout=5)
            @test next_step!(link).n == 7 && next_step!(link).n == 8
            @test_throws EOFError next_step!(link)
            close_link(link)
            wait(task)
        finally
            rm(path; force=true)
        end
    end
end
