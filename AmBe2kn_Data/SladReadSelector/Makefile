TARGET = HistFromSlad
SRCS = HistFromSlad.C
DICTNAME = SladDict
SOURCES = SladReadSelector/SladReadSelector.C 
HEADERS = SladReadSelector/SladReadSelector.h 
OBJS = $(addsuffix .o, $(notdir $(basename $(SRCS))))
SOBJS = $(addsuffix .o, $(notdir $(basename $(SOURCES))))
#ROOTSYS = /Applications/root
ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS   = $(shell $(ROOTSYS)/bin/root-config --libs) -lProof -lMinuit #-lMemStat
ROOTLDFLAGS = $(shell $(ROOTSYS)/bin/root-config --ldflags)
ROOTCINT = $(ROOTSYS)/bin/rootcint
CXXFLAGS   = $(ROOTCFLAGS) -Wall -fPIC -g -Wno-unused-function #-pthread 
CXXLIBS    = $(ROOTLIBS)
CC = g++ #g++

all: $(DICTNAME).so $(TARGET) 

%: %.o
	$(CC) $(CXXFLAGS) $(CXXLIBS) $(DICTNAME).so $< -o $@

$(OBJS): $(SRCS)
	$(CC) $(CXXFLAGS) -c $^

$(SOBJS): $(SOURCES)
	$(CC) $(CXXFLAGS) -c $^

$(DICTNAME).cxx: 
	$(ROOTCINT) -f $@ -c -p $(HEADERS) LinkDef.h

$(DICTNAME).so: $(DICTNAME).cxx $(SOBJS)
	$(CC) -shared -o $@ $(ROOTLDFLAGS) $(CXXFLAGS) $(CXXLIBS) -I$(ROOTSYS)/include $^

clean:
	rm  -r  $(OBJS) $(SOBJS) $(DICTNAME).cxx 
