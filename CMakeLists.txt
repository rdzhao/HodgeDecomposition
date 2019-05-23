cmake_minimum_required(VERSION 3.1)

project(DecomAll)

find_package(OpenGL)

find_package(Qt5Widgets)
find_package(Qt5OpenGL)
set(CMAKE_AUTOMOC ON)

find_package(CGAL)

include(${CGAL_USE_FILE})

include_directories(.)
include_directories(eigen)
include_directories(glm)

file(GLOB Headers "src/*.h")
file(GLOB Sources "src/*.cpp")
file(GLOB Resources "resources/*.png")

add_executable(DecomAll ${Headers} ${Sources} ${Resources})

qt5_use_modules( DecomAll Widgets)
qt5_use_modules( DecomAll OpenGL)

target_link_libraries( DecomAll ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
target_link_libraries( DecomAll ${QT_LIBRARIES} ${OPENGL_LIBRARIES} )