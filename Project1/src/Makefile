LINK_TARGET = ns.exe
OBJS =                    \
  lib/lib.o               \
  lib/NucleonScattering.o \
  lib/main.o

REBUILDABLES = $(OBJS) $(LINK_TARGET)

all: $(LINK_TARGET)

$(LINK_TARGET): $(OBJS)
	g++ -o $@ $^


lib/%.o: %.cpp
	g++ -o $@ -c $<


clean:
	rm -f $(REBUILDABLES)
