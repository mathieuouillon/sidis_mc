project('simulation', ['c', 'cpp', 'fortran'],
        version : '1.0',
        default_options : ['warning_level=0',
                          'optimization=3',
                          'debug=false'])

# Compiler flags
fortran_args = ['-O3', '-w', '-g']
cpp_args = ['-g', '-O3']

#── ROOT via root-config ─────────────────────────────────────────────────────#
root_incdir   = run_command('root-config', '--incdir', check: true).stdout().strip()
root_libs     = run_command('root-config', '--libs',   check: true).stdout().strip().split()

root_dep = declare_dependency(
  include_directories: include_directories(root_incdir),
  link_args:           root_libs,
)

# Find gfortran library
gfortran_dep = meson.get_compiler('cpp').find_library('gfortran', required : true)

# CERNLIB path and libraries
cernlib_libdir = '/cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/cernlib/2023/lib'

# Find CERNLIB libraries
mathlib_dep = meson.get_compiler('cpp').find_library('mathlib', dirs : cernlib_libdir, required : true)
pawlib_dep = meson.get_compiler('cpp').find_library('pawlib', dirs : cernlib_libdir, required : true)
packlib_dep = meson.get_compiler('cpp').find_library('packlib', dirs : cernlib_libdir, required : true)
kernlib_dep = meson.get_compiler('cpp').find_library('kernlib', dirs : cernlib_libdir, required : true)

# Source files
fortran_sources = [
  'simulation.f90',
  'fermimotion.f90',
  'pythia6428.f',
  'pythiaconfig.f90',
  'book.f90',
  'transfo.f90',
  'common_module.f90',
  'nucdens.f',
  'density.f90',
  'qweight.f90',
  'accep_fun.f',
  'clas_at12g.f',
  'clas12_accept.f',
  'read_par_clas12g.f',
  'smear_fun.f',
  'vz_position.f90',
  'simulation_helpers.f90',
  'cli.f90'
]

# All dependencies
all_deps = [
  root_dep,
  gfortran_dep,
  mathlib_dep,
  pawlib_dep,
  packlib_dep,
  kernlib_dep
]

# Create executable
sim_exe = executable('simulation',
          sources : [fortran_sources],
          dependencies : all_deps,
          fortran_args : fortran_args,
          cpp_args : cpp_args,
          install : true)

subdir('tests')