
check_PROGRAMS = \
  test_particle_storage \
  test_weighting
  
TESTS = \
  test_particle_storage

AM_CXXFLAGS = -I../src


test_particle_storage_SOURCES = test_particle_storage.cpp ../src/particles.cpp
test_particle_storage_LDADD = -L../src 

test_weighting_SOURCES = test_weighting.cpp