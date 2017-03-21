comp = ifort
link = xilink

bindir = ../../../bin

cflags = /c /I../../../include /debug:full /pdb /traceback #/iface:stdref
lflags = /nodefaultlib:libcmtd /dll

src=precision.f90 xedriver.f90 xesolver.f90

obj = $(src:.f90=.obj)

addlib = $(bindir)/libint304_r8_deb.lib

target = $(bindir)/xesolver_deb.dll

$(obj): $(src)
    $(comp) $(cflags) $(src)

$(target): $(obj)
    $(link) $(lflags) $(obj) $(addlib) /out:$(target)

all: $(target)

clean:
    del *.obj
    del *.mod
    del *.lib