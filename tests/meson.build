# Tests build file for FARM simulation

# Common module source (needed by most tests)
common_mod_src = files('../src/modules/common_module.f90')

# ── Unit test: Lorentz boost (TL) ──────────────────────────────────────────
test_tl_exe = executable('test_tl',
    sources: ['test_tl.f90', 'tl_only.f90'],
    fortran_args: fortran_args,
    install: false)
test('unit_lorentz_boost', test_tl_exe)

# ── Unit test: Nuclear initialization (InitNucl) ──────────────────────────
test_initnucl_exe = executable('test_initnucl',
    sources: ['test_initnucl.f90', '../src/physics/fermimotion.f90', common_mod_src,
              '../src/rng_bridge.cpp'],
    dependencies: all_deps,
    fortran_args: fortran_args,
    cpp_args: cpp_args,
    install: false)
test('unit_init_nuclear', test_initnucl_exe,
    workdir: meson.project_source_root())

# ── Unit test: Vertex z-position ──────────────────────────────────────────
test_vz_exe = executable('test_vz',
    sources: ['test_vz.f90', '../src/physics/vz_position.f90'],
    fortran_args: fortran_args,
    install: false)
test('unit_vertex_z', test_vz_exe)

# ── Golden file tests (require built simulation) ──────────────────────────
bash = find_program('bash')
sim_path = sim_exe.full_path()

foreach config : ['proton',
                  'deuterium', 'tritium', 'he3', 'he4', 'li6', 'li7',
                  'carbon', 'aluminum', 'iron', 'tin', 'lead',
                  'neon', 'krypton', 'xenon', 'copper',
                  'iron_cs', 'carbon_hs']
    test('golden_' + config,
        bash,
        args: [meson.project_source_root() / 'tests/test_golden.sh',
               config, sim_path],
        workdir: meson.project_source_root(),
        timeout: 300,
        is_parallel: false,
        depends: sim_exe)
endforeach

# ── Statistical distribution tests ────────────────────────────────────────
python = find_program('python3', required: false)
if python.found()
    test('statistical_distributions',
        python,
        args: [meson.project_source_root() / 'tests/test_distributions.py',
               '--sim', sim_path,
               '--workdir', meson.project_source_root()],
        timeout: 600,
        is_parallel: false,
        depends: sim_exe)
endif
