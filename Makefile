CC = g++
RM = rm -f
LIBS = -lm -llapack -lblas 
CFLAGS = -Wall

TARGET = main
OBJS = main.o hmm.o
all: $(TARGET)

.c.o:
	$(CC) -g -c $(CFLAGS)  $<
$(TARGET): $(OBJS)
	$(CC) -g -o $@ $^ $(LIBS) 

clean:
	$(RM) $(TARGET) *#* *~ *.o