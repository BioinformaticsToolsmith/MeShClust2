all: bin/Red.o bin/meshclust2

bin/Red.o:
	mkdir -p bin
	mkdir -p bin/exception
	mkdir -p bin/nonltr
	mkdir -p bin/utility
	$(MAKE) -C src
bin/meshclust2: bin/Red.o
	$(MAKE) -C src/cluster
	cp src/cluster/meshclust2 bin

clean:
	$(MAKE) clean -C src
	$(MAKE) clean -C src/cluster
	$(RM) -r bin

rebuild: clean all
.PHONY: all clean
