mcMd_mcMoves_grandCanon_=\
    mcMd/mcMoves/grandCanon/PolymerPointExchangeMove.cpp \
    mcMd/mcMoves/grandCanon/PointInsertionMove.cpp

mcMd_mcMoves_grandCanon_SRCS=\
    $(addprefix $(SRC_DIR)/, $(mcMc_mcMoves_grandCanon_))
mcMd_mcMoves_grandCanonc_OBJS=\
    $(addprefix $(BLD_DIR)/, $(mcMd_mcMoves_grandCanon_:.cpp=.0))
