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

# Fortran physics sources (no main program)
fortran_sources = [
  # Modules
  'src/modules/common_module.f90',
  'src/modules/simulation_helpers.f90',

  # Physics
  'src/physics/fermimotion.f90',
  'src/physics/pythiaconfig.f90',
  'src/physics/book.f90',
  'src/physics/transfo.f90',
  'src/physics/nucdens.f90',
  'src/physics/density.f90',
  'src/physics/qweight.f90',
  'src/physics/vz_position.f90',

  # ISO_C_BINDING interface
  'src/bindings/c_interface.f90',

  # External libraries
  'src/external/pythia6428.f',
]

# C++ main program
cpp_sources = [
  'src/main.cpp',
  'src/rng_bridge.cpp',
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

# Include directories
inc_dirs = include_directories('src')

# Create executable (C++ main + Fortran physics)
sim_exe = executable('simulation',
          sources : [fortran_sources, cpp_sources],
          include_directories : inc_dirs,
          dependencies : all_deps,
          fortran_args : fortran_args,
          cpp_args : cpp_args,
          install : true)

subdir('tests')