
SRC1 = gg_hex_v2.f90
TARGET1 = structMeshGen
CC = gfortran
FLAGS = 

$(TARGET1) : $(SRC1)
	$(CC) $(FLAGS) $(SRC1) -o $(TARGET1)

all : $(TARGET1) 

clean :
	-rm -rf $(TARGET1)  *.o core
 
