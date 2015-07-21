TARGET = ReconMain
CC = g++
ROOTCINT = $(ROOTSYS)/bin/rootcint
DICTNAME = ReconSelector_dict
SRCS = $(addsuffix .C, $(TARGET))
DIR = ./ReconSelector
SOURCES = $(DIR)/ReconSelector.C
HEADERS = $(DIR)/ReconSelector.h
OBJS = $(addsuffix .o, $(notdir $(basename $(SRCS))))
SOBJS = $(addsuffix .o, $(notdir $(basename $(SOURCES))))
ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags) 
ROOTLIBS   = $(shell $(ROOTSYS)/bin/root-config --libs) -lProof -lProofPlayer
ROOTLDFLAGS = $(shell $(ROOTSYS)/bin/root-config --ldflags)

all: $(DICTNAME).C
	$(CC) $(ROOTCFLAGS) $(ROOTLIBS) $(SRCS) $(SOURCES) $^ -o $(TARGET) 

$(DICTNAME).C: $(HEADERS)
	$(ROOTCINT) -f $@ -c $^

.PHONY: clean
	
clean: 
	rm -f $(DIR)/*.so