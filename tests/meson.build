# Tests build file for FARM simulation

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
