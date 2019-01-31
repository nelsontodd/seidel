#
# == Paths ==
#
BIN_DIR   := bin
LIB_DIR   := lib
BUILD_DIR := build

#
# == Files ==
#
LIBS := $(wildcard $(LIB_DIR)/*.cpp)
OBJS = $(LIBS:.cpp=.o)
#TODO: Do this more intelligently
BUILD_OBJS = $(foreach obj,$(OBJS),$(BUILD_DIR)/$(notdir $(obj)))

#
# == CC Flags ==
#
CC      := c++

#
# == Targets ==
#

default: iterativeLA
objs: $(OBJS)


clean:
	$(RM) $(BUILD_DIR)/*.o $(BIN_DIR)/*

%.o: %.cpp
	$(CC) -o $(BUILD_DIR)/$(notdir $@) -c $<

iterativeLA: $(OBJS) 
	$(CC) -o $(BIN_DIR)/iterativeLA $(BUILD_OBJS)
