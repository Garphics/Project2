#My makefile
all: part1 raytracer
part1:
	g++ -o part1.exe part1.cpp -O2 -L/usr/X11R6/lib -lm -lpthread -lX11

raytracer:
	g++ -o part2.exe raytracer.cpp -O2 -L/usr/X11R6/lib -lm -lpthread -lX11

