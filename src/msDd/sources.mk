include $(SRC_DIR)/msDd/simulation/sources.mk
include $(SRC_DIR)/msDd/param/sources.mk

msDd_SRCS=$(msDd_param_SRCS) $(msDd_simulation_SRCS) 

msDd_OBJS=$(msDd_SRCS:.cpp=.o)

#$(msDd_LIB): $(msDd_OBJS)
#	$(AR) rcs $(msDd_LIB) $(msDd_OBJS)

