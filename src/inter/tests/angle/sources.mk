inter_tests_angle_=inter/tests/angle/Test.cpp


inter_tests_angle_SRCS=\
     $(addprefix $(SRC_DIR)/, $(inter_tests_angle_))
inter_tests_angle_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(inter_tests_angle_:.cpp=.o))
