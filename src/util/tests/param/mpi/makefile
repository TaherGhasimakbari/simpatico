BLD_DIR_REL =../../../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/util/tests/param/mpi/sources.mk

all: $(util_tests_param_mpi_OBJS)

clean:
	rm -f $(util_tests_param_mpi_OBJS) 
	rm -f $(util_tests_param_mpi_OBJS:.o=.d)
	rm -f $(util_tests_param_mpi_OBJS:.o=)

-include $(util_tests_param_mpi_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

