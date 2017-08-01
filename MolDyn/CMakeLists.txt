cmake_minimum_required (VERSION 2.8.11)
project (MOLDYN)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2 -Wall")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lm")
add_executable (moldyn
                Andersen.c
                RDF.c
                acf.c
                atomic_numbers.c
                export.c
                jmol_func.c
                linked_cells.c
                lj_forces.c
                lj_forces_cells.c
                main.c
                parameters.c
                paraview_func.c
                random_gen.c
                read_elements.c
                sim_aneal.c
                tersoff_forces_cells.c
                tersoffopt.c
                vrescale.c
                xi_define.c)
