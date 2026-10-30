"""Writes the fixture of test/perception_link_test.jl with smearPerception's own code.

steps.h5: three step records in the format of smearPerception's `scripts/reconstruct_bag.py
--save-field` (its `ReconstructionWriter`, repacked without chunking to keep the file small), on a
grid of distinct sizes per axis, so a swapped axis cannot pass. stream.bin: the bytes its
`fe_link.StepPublisher` sends for the same records (hello, three steps, END). expected.h5:
`voxels.sample_field` of each step's float32 field at query points [m].

The records (the Julia test checks these values): n = 7, 8, 9 with dt_ns = 100 and t0_ns = 300;
wrench [0.1, 0, 0, 0.01, 0.02, 0.03], [1, 2, 2, 0.1, 0.2, 0.3], NaN [N, N m]; contact threshold
0.5 N; probe pose: rotation of 30° about z, translation (10, 20, 30 + s) mm, NaN at the last step.

Run in smearPerception's env (smear_ros), from this folder:
    python make_fixture.py <path to smearPerception>
"""
import importlib.util
import os
import socket
import sys
from pathlib import Path

import h5py
import numpy as np
from smear_perception.fe_link import StepPublisher, StepRecord, contact_state
from smear_perception.reconstruction.pipeline import StepResult
from smear_perception.reconstruction.voxels import compute_voxel_grid, sample_field

here = Path(__file__).parent
spec = importlib.util.spec_from_file_location("rb", Path(sys.argv[1]) / "scripts" / "reconstruct_bag.py")
rb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rb)

rng = np.random.default_rng(0)
grid = compute_voxel_grid((0.012, -0.007), [[-0.002, -0.0025, 0], [0.002, 0.0025, 0.008]], 0.001)
assert len(set(grid["shape"])) == 3
writer = rb.ReconstructionWriter(here / "chunked.h5", grid, {"dt_ns": 100}, save_field=True, with_gt=False)
path = f"/tmp/smear_fixture_{os.getpid()}.sock"
publisher = StepPublisher(grid, dt_ns=100, contact_threshold=0.5, path=path)
fe = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
fe.connect(path)  # the three steps fit the socket buffer

points = np.concatenate([grid["x0"] + rng.uniform([-0.0025, -0.0035, -0.0005], [0.0025, 0.0035, 0.0085], (300, 3)),
                         grid["centers"][-1, -1, -1][None], grid["centers"][0, 0, 0][None]])
wrenches = [[0.1, 0, 0, 0.01, 0.02, 0.03], [1, 2, 2, 0.1, 0.2, 0.3], [np.nan] * 6]
c, s30 = np.cos(np.pi / 6), np.sin(np.pi / 6)
expected = []
for s in range(3):
    phi = rng.normal(0, 0.002, grid["shape"])
    phi[rng.random(grid["shape"]) < 0.05] = np.nan
    count = np.where(np.isnan(phi), 0, rng.integers(1, 5, grid["shape"])).astype(np.int32)
    n = 10 + 7 * s
    pose = np.array([[c, -s30, 0, 0.010], [s30, c, 0, 0.020], [0, 0, 1, 0.030 + 0.001 * s], [0, 0, 0, 1]])
    wrench = np.array(wrenches[s], dtype=float)
    record = StepRecord(7 + s, 300 + (7 + s) * 100, 4 - s,
                        StepResult(phi, count, rng.normal(size=(n, 3)), rng.normal(size=(n, 3))),
                        wrench, contact_state(wrench, 0.5), np.full((4, 4), np.nan) if s == 2 else pose)
    writer.add(record)
    publisher.publish(record)
    expected.append(sample_field(phi.astype(np.float32).astype(float), points, grid))
writer.close()
publisher.close()
with fe, fe.makefile("rb") as f:
    (here / "stream.bin").write_bytes(f.read())

with h5py.File(here / "chunked.h5") as src, h5py.File(here / "steps.h5", "w") as dst:
    dst.attrs.update(src.attrs)
    for name, ds in src.items():
        dst[name] = ds[()]
(here / "chunked.h5").unlink()
with h5py.File(here / "expected.h5", "w") as f:
    f["points"] = points
    f["samples"] = np.array(expected)
