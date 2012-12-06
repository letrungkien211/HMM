CC = g++
RM = rm -f
LIBS = -lm
CFLAGS = -Wall -g 
INCLUDES = -I/usr/include/eigen3

TARGET = main
OBJS = hmm.o main.o
all: $(TARGET)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS) 

clean:
	$(RM) $(TARGET) *#* *~ *.o