# Test Failures Triage Report

**Summary:** 19 passed, 2 failed, 4 errored (total 25 tests)

## Failures Identified

### 1. **File Not Found: `test/optimization_test stokes.jl`** [CRITICAL]
   - **Error:** `SystemError: opening file "...test/optimization_test stokes.jl": No such file or directory`
   - **Location:** `test/runtests.jl:7`
   - **Root Cause:** Typo in filename—spaces in name should likely be underscores
   - **Fix:** Correct the filename in `test/runtests.jl` from `optimization_test stokes.jl` to `stokes_optimization_test.jl` (or skip if file doesn't exist)
   - **Priority:** HIGH—blocks entire test run

### 2. **Missing HDF5 Fixture: `cylindergen/cylinder.h5`** [ERROR]
   - **Error:** `unable to determine if .../cylindergen/cylinder.h5 is accessible in the HDF5 format (file may not exist)`
   - **Location:** `test/IGA_tests/extraction_test.jl:9`
   - **Root Cause:** Test file requires external HDF5 fixture not present in repo
   - **Fix:** Either skip the test (conditional on fixture existence) or create a minimal fixture generator
   - **Priority:** MEDIUM—affects IGA tests only

### 3. **Linear Elasticity Test Error** [ERROR]
   - **Error:** Likely import/API mismatch in `test/linear_elasticity_test.jl`
   - **Location:** `test/runtests.jl:5`
   - **Root Cause:** Unclear from output; likely missing function or signature mismatch
   - **Fix:** Run test directly to see detailed error
   - **Priority:** MEDIUM

### 4. **Stokes Test Error** [ERROR]
   - **Error:** Stokes test initialization fails
   - **Location:** `test/stokes_test.jl:9`
   - **Root Cause:** Likely API or import issue
   - **Fix:** Check for missing boundary condition or model setup function
   - **Priority:** MEDIUM

### 5. **Aqua.jl Code Quality Issues** [FAILED - 2/6 tests]

   a. **Undefined Exports** [FAILED]
      - 34 symbols are used internally but not exported in `__init__.jl`
      - Examples: `apply_boundary_conditions_stokes`, `assemble_system`, `simulate_stokes`, etc.
      - **Fix Options:**
        - (A) Add missing exports to module definition (preferred if API is stable)
        - (B) Add `Aqua.test_undefined_exports(smearFEM; broken=true)` to mark as known issue
      - **Priority:** LOW (code quality, not functional)

   b. **Stale Dependencies** [FAILED]
      - Unused dependencies detected in `Project.toml`:
        - BenchmarkTools, Images, Aqua, ProfileView, PlotlyJS
      - **Fix:** Remove from `Project.toml` or move to `extras` if conditionally used
      - **Priority:** LOW (code quality)

---

## Categorized Action Plan

### BLOCKING (Fix First)
1. **Fix filename typo in `test/runtests.jl`**
   - Change line with bad filename reference
   - Verify all referenced files exist

### HIGH (Impact on Core Tests)
2. **Resolve linear elasticity and stokes test errors**
   - Run each test file directly to capture detailed error
   - Apply minimal API shims if needed

### MEDIUM (Impact on Specific Suites)
3. **Handle missing HDF5 fixture**
   - Conditional skip in `extraction_test.jl` if fixture not found

### LOW (Code Quality, Non-Blocking)
4. **Address Aqua failures (optional)**
   - Either export symbols or mark as known issues
   - Clean up unused dependencies from `Project.toml`

---

## Next Steps
1. Fix the filename typo immediately (should unblock many tests)
2. Run core tests individually to collect specific error details
3. Apply minimal fixes per category
4. Re-run full test suite
