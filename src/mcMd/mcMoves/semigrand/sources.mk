mcMd_mcMoves_semigrand_=

ifdef INTER_BOND
mcMd_mcMoves_semigrand_+= \
    mcMd/mcMoves/semigrand/HomopolymerSemiGrandMove.cpp \
    mcMd/mcMoves/semigrand/GeneralpolymerSemiGrandMove.cpp \
    mcMd/mcMoves/semigrand/WangLandauMove.cpp \
    mcMd/mcMoves/semigrand/WangLandauAdaptiveStepMove.cpp
endif

mcMd_mcMoves_semigrand_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_mcMoves_semigrand_))
mcMd_mcMoves_semigrand_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_mcMoves_semigrand_:.cpp=.o))

