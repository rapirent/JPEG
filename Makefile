all: encoder decoder

.PHONY : clean

encoder: encoder.c util.h bitmap_image.hpp
	g++ -o encoder encoder.c -O3

decoder: decoder.c util.h bitmap_image.hpp
	g++ -o decoder decoder.c -O3

clean:
	rm encoder decoder
