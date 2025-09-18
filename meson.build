project('simulation', ['c', 'cpp', 'fortran'],
        version : '1.0',
        default_options : ['warning_level=0',
                          'optimization=3',
                          'debug=true'])

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
cernlib_path = '/cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/cernlib'
cernlib_libdir = cernlib_path / '2023/lib'

# Find CERNLIB libraries
mathlib_dep = meson.get_compiler('cpp').find_library('mathlib', 
                                                    dirs : cernlib_libdir, 
                                                    required : true)
pawlib_dep = meson.get_compiler('cpp').find_library('pawlib', 
                                                   dirs : cernlib_libdir, 
                                                   required : true)
packlib_dep = meson.get_compiler('cpp').find_library('packlib', 
                                                    dirs : cernlib_libdir, 
                                                    required : true)
kernlib_dep = meson.get_compiler('cpp').find_library('kernlib', 
                                                    dirs : cernlib_libdir, 
                                                    required : true)

# Source files
fortran_sources = [
  'simulation.f90',
  'fermimotion.f',
  'pythia6428.f',
  'pythiaconfig.f',
  'book.f',
  'transfo.f',
  'fermimotion2.f',
  'nucdens.f',
  'density.f',
  'qweight.f',
  'accep_fun.f',
  'clas_at12g.f',
  'clas12_accept.f',
  'read_par_clas12g.f',
  'smear_fun.f'
]

cpp_sources = [
  'ALERT_fastMC.cc',
 # 'main.cc'
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
executable('simulation',
          sources : [fortran_sources, cpp_sources],
          dependencies : all_deps,
          fortran_args : fortran_args,
          cpp_args : cpp_args,
          install : true)

# Custom target to show compiler version (equivalent to 'gfortran --version' in makefile)
run_target('show-version',
          command : [find_program('gfortran'), '--version'])