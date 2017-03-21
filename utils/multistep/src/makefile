COMPILER = ifort
LINK =     xilink

SRC_DLL =       TOOLS.F90\
                BDF.F90\
                SDBDF.F90\
                HSIZE.F90\
                INTER.F90\
                REJECT.F90\
                INTEG.F90\
#                LSODA.F90\
                IEDS.F90\
                INITSTEP.F90\

SRC_MOD  =      PARS.F90
SRC_MOD4 =      PRECISION4.F90 
SRC_MOD8 =      PRECISION8.F90 

ADD_LIBS = ..\..\..\bin\liblsoda_deb.lib

TAR_DLL4 = ..\..\..\bin\INT304_R4_deb.DLL
TAR_DLL8 = ..\..\..\bin\INT304_R8_deb.DLL

TAR_LIB4 = ..\..\..\bin\libINT304_R4_deb.LIB
TAR_LIB8 = ..\..\..\bin\libINT304_R8_deb.LIB


MOD = $(SRC_MOD:.F90=.obj)
MOD4 = $(SRC_MOD4:.F90=.obj)
MOD8 = $(SRC_MOD8:.F90=.obj)
OBJ = $(SRC_DLL:.F90=.obj)

#OPT = -c -fixed /iface:stdref -debug:full -check:all
#OPT = -c -fixed /iface:stdref /check:all /debug:full  /traceback /warn:unused /warn:argument_checking /warn:declarations /fpe:0 
OPT = -c -fixed  /traceback /heap-arrays:10240 #/iface:stdref
lflags = /nodefaultlib:libcmt
#/check:all /debug:full 

$(MOD): $(SRC_MOD) 
	$(COMPILER) $(OPT) $(SRC_MOD)

$(MOD4): $(SRC_MOD4) 
	$(COMPILER) $(OPT) $(SRC_MOD4) 

$(MOD8): $(SRC_MOD8) 
	$(COMPILER) $(OPT) $(SRC_MOD8) 

$(OBJ): $(SRC_DLL) 
	$(COMPILER) $(OPT) $(SRC_DLL)


#ALL: $(TAR_DLL4) $(TAR_DLL8) $(TAR_LIB4) $(TAR_LIB8)


DLL4: $(TAR_DLL4)

$(TAR_DLL4): $(MOD4) $(MOD) $(OBJ) 
	$(LINK) -dll $(lflags) $(ADD_LIBS) $(MOD) $(MOD4) $(OBJ) /out:$(TAR_DLL4) 

DLL8: $(TAR_DLL8)

$(TAR_DLL8): $(MOD8) $(MOD) $(OBJ) 
	$(LINK) -dll $(lflags) $(ADD_LIBS) $(MOD) $(MOD8) $(OBJ) /out:$(TAR_DLL8) 

LIB4: $(TAR_LIB4)

$(TAR_LIB4): $(MOD4) $(MOD) $(OBJ) 
	$(LINK) -lib $(lflags) $(ADD_LIBS) $(MOD) $(MOD4) $(OBJ) /out:$(TAR_LIB4) 

LIB8: $(TAR_LIB8)

$(TAR_LIB8): $(MOD8) $(MOD) $(OBJ) 
	$(LINK) -lib $(lflags) $(ADD_LIBS) $(MOD) $(MOD8) $(OBJ) /out:$(TAR_LIB8) 

CLEAN:
	del *.obj
#    del *.lib
    del *.mod
#    del *.exp
    del *.pdb 



