gen: HBeam.h pion_generator.cpp
	g++  -Wall -L/Users/tujuba/panda/svn/ext-apr13/lib/root -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lMinuit -lPhysics -lMathCore -lThread -lpthread -lm -ldl -L/Users/tujuba/panda/svn/ext-apr13/lib -lPluto -I/Users/tujuba/panda/svn/ext-apr13/include/root -I/Users/tujuba/panda/svn/ext-apr13/include/pluto  pion_generator.cpp -o gen

gen_gprof: HBeam.h pion_generator.cpp
	g++ -Wall -pg -o gen_gprof -L/Users/tujuba/panda/svn/ext-apr13/lib/root -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lMinuit -lPhysics -lMathCore -lThread -lpthread -lm -ldl -L/Users/tujuba/panda/svn/ext-apr13/lib -lPluto -I/Users/tujuba/panda/svn/ext-apr13/include/root -I/Users/tujuba/panda/svn/ext-apr13/include/pluto  pion_generator.cpp 
