c++ -o lhe_conversion_HZZ `root-config --glibs --cflags` -lm lhe_conversion_HZZ.cpp  
#gcc -std=c++14 -g  trial.cpp TUtil.cc -I$ROOTSYS/include `root-config --libs ` -lMinuit -o trial
./lhe_conversion_HZZ outfile.lhe

