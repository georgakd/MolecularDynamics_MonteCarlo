# *************************************************************
# Makefile for helper scripts for MD/MC
# *************************************************************
CC = gcc

ALTEREXE = alter
ELONGEXE = elong
REMOVEEXE = remove
SIEXE = Si
EXPANDEXE = expand
PROJEXE = proj

#################### Options ###################################
CFLAGS = -O2 -g -Wall -c

#################### Exes ######################################

PROGS = $(ALTEREXE) $(ELONGEXE) $(REMOVEEXE) $(SIEXE) $(EXPANDEXE) $(PROJEXE)
 
################### Targets #####################################
all: $(PROGS)
$(ALTEREXE):
	$(CC) $(CFLAGS) alter.c
	$(CC) -o $(ALTEREXE) alter.o
$(ELONGEXE):
	$(CC) $(CFLAGS) elongation.c
	$(CC) -o $(ELONGEXE) elongation.o
$(REMOVEEXE):
	$(CC) $(CFLAGS) remove.c
	$(CC) -o $(REMOVEEXE) remove.o
$(SIEXE):
	$(CC) $(CFLAGS) Si.c
	$(CC) -lm -o $(SIEXE) Si.o
$(EXPANDEXE):
	$(CC) $(CFLAGS) xy_expand.c
	$(CC) -o $(EXPANDEXE) xy_expand.o
$(PROJEXE):
	$(CC) $(CFLAGS) z_projection.c
	$(CC) -o $(PROJEXE) z_projection.o		
	
clean: 
	rm -rf *.o
	rm -rf $(ALTEREXE)
