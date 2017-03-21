comp=ifort
link=xilink

tar=../../../bin/liblsoda_deb.lib

src=opkda1.f opkda2.f opkdmain.f
obj=$(src:.f=.obj)

cflags = -c -fixed /iface:stdref
dflags = /debug:full /traceback
lflags = /lib /nodefaultlib:libcmt

all: $(tar) 

$(obj): $(src)
    $(comp) $(cflags) $(dflags) $(src)

$(tar): $(obj)
    $(link) $(lflags) $(obj) -out:$(tar)


clean:
    del *.mod *.obj *.lib