project('simulation', ['c', 'cpp', 'fortran'],
        version : '1.0',
        default_options : ['warning_level=2',
                          'optimization=3',
                          'debug=true',
                          'cpp_std=c++20'])

# Compiler flags
fortran_args = ['-O3', '-w', '-g']
cpp_args = ['-g', '-O3', '-Wall', '-Wextra']

# Find gfortran library (needed to link PYTHIA6 Fortran object)
gfortran_dep = meson.get_compiler('cpp').find_library('gfortran', required : true)

# Source files
fortran_sources = [
  'src/external/pythia6428.f',   # PYTHIA 6.4.28 (only Fortran)
]

cpp_sources = [
  'src/main.cpp',                # C++ main program
]

# Dependencies
all_deps = [
  gfortran_dep,
]

# Include directories
inc_dirs = include_directories('src')

# Create executable
sim_exe = executable('simulation',
          sources : [fortran_sources, cpp_sources],
          include_directories : inc_dirs,
          dependencies : all_deps,
          fortran_args : fortran_args,
          cpp_args : cpp_args,
          install : true)

subdir('tests')
