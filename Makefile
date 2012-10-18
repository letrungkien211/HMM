CC = g++
RM = rm -f
LIBS = -lm -llapack -lblas
CFLAGS = -Wall

TARGET = main
OBJS = main.o hmm.o
all: $(TARGET)

.c.o:
	$(CC) -c $(CFLAGS) $<
$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(LIBS)

clean:
	$(RM) $(TARGET) *#* *~ *.o