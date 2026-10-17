CHPL = chpl
CFLAGS = --fast
SRCCOM = main.chpl constants.chpl

TARGET = cavity

all: $(TARGET)

cavity: cavity.chpl $(SRCCOM)
	$(CHPL) $(CFLAGS) -o $@ $^

clean:
	rm -f $(TARGET)
