ddMd_tests_potentials_=ddMd/tests/potentials/Test.cpp

ddMd_tests_potentials_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_potentials_))
ddMd_tests_potentials_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_tests_potentials_:.cpp=.o))

