cmake_minimum_required(VERSION 2.6)

#Configuration du projet
project(Klib)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BUILD_TYPE})

#Inclusion de Boost
include_directories(.)
link_directories(.)

#Configuration de l'exécutable
file(
	GLOB
	source_files
	*.cpp *.h
)
add_executable(
	Klib
	${source_files}
)

#Configuration de l'édition de liens
target_link_libraries(
	Klib
)





