mcMd_mcMoves_micelle_=\
    mcMd/mcMoves/micelle/AggregatorMove.cpp \
    mcMd/mcMoves/micelle/RgUmbrellaMove.cpp

mcMd_mcMoves_micelle_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_mcMoves_micelle_))
mcMd_mcMoves_micelle_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_mcMoves_micelle_:.cpp=.o))

