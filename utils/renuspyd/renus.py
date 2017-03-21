# -*- coding: utf-8 -*-
import os
from ctypes import *
from libhandle import libhandle
import numpy as np 
from numpy import array
import h5py
import copy

def frombuffer(loc=None, dim=None, tp=None):
    """
    """
    if (loc == 0 or dim == 0):
        print "WARNING: Can't get data from buffer"
        value = np.array([],dtype="double")
    else:
        buf = np.core.multiarray.int_asbuffer(loc, dim)
        value = np.frombuffer(buf, tp)
    return value


class XeODE(libhandle):
    def __init__(self):
        """
        """
        libhandle.__init__(self,dirdll="../../bin", dlls=["xesolver_deb"])
        mainlib = self.dlls.keys()[-1]
        self.odesolver = self.dlls[mainlib].xenon_odes

        p_c_double = POINTER(c_double)
        p_c_int = POINTER(c_int)

        CMPFUNC = CFUNCTYPE(None, p_c_double, p_c_double, p_c_double, 
            p_c_double, p_c_double, p_c_double, p_c_double,
            p_c_double, p_c_double, p_c_double, p_c_double,
            p_c_int, p_c_double, p_c_double)

        self.cmpXeODESolver = CMPFUNC(self.odesolver) #  Make sure you keep references to CFUNCTYPE() objects as long as they are used from C code !!

class RENUS(libhandle):
    def __init__(self):
        """
        """
        libhandle.__init__(self,dirdll="../../bin", dlls=["renus"])
        "xesolver_deb"
        mainlib = self.dlls.keys()[-1]
#       self.hello = self.dlls[mainlib].hello
        self.libload = self.dlls[mainlib].renusinput
        self.libinit = self.dlls[mainlib].RENUSINIT
        self.libstat = self.dlls[mainlib].renustate
        self.libqtrans = self.dlls[mainlib].QTRANS
        self.lib_set_dumpfunc = self.dlls[mainlib].set_dumpfunc
        self.lib_set_xeodefunc = self.dlls[mainlib].set_xeodefunc
        self.lib_join2data = self.dlls[mainlib].join2data
        self.libtran = self.dlls[mainlib].renustransient

        self.stage = 1

    def run(self, inpfile=None, dumpfile=None, dtime=[], mode="qtr"):
        """
        """
        self.load(filename=inpfile)
        self.libinit()
        self.join2var()        
        if dumpfile: self.dumph5(dumpfile)
        self.libstat()        

        if len(dtime)>0:
            dtime = array(dtime, dtype="double")
            len_dtime = len(dtime)
            lendt_ = c_int(len_dtime)
            if(mode == "qtr"): self.libqtrans(lendt_, dtime.ctypes.data)

        if(mode == "tr"): self.libtran()


    def load(self,filename):
        """Reading of input file
        """
        fname=filename.ljust(255)
        fnlen=c_int(len(fname))
        self.libload(c_char_p(fname))

    def set_dumpfunc(self,func):
        """
        """
        CMPFUNC = CFUNCTYPE(None)
        self.cmpfunc = CMPFUNC(func) #  Make sure you keep references to CFUNCTYPE() objects as long as they are used from C code !!
        self.lib_set_dumpfunc(self.cmpfunc)

    def set_xeodesolver(self,cmpfunc):
        """
        """
        self.lib_set_xeodefunc(cmpfunc)



    def join2var(self):
        """Join to variable from dll
        """
        lenptrs = 100
        # address of arrays 
        ct_i32ptrs = (c_long*lenptrs)(0)
        ct_i64ptrs = (c_long*lenptrs)(0)
        ct_f32ptrs = (c_long*lenptrs)(0)
        ct_f64ptrs = (c_long*lenptrs)(0)

        # lengths of arrays
        ct_i32lens = (c_long*lenptrs)(0)
        ct_i64lens = (c_long*lenptrs)(0)
        ct_f32lens = (c_long*lenptrs)(0)
        ct_f64lens = (c_long*lenptrs)(0)

        # get address and length from dll by the name of variable
        self.lib_join2data(byref(ct_i32ptrs), byref(ct_i32lens), byref(ct_i64ptrs), byref(ct_i64lens), 
                           byref(ct_f32ptrs), byref(ct_f32lens), byref(ct_f64ptrs), byref(ct_f64lens))

        # int32               
        keys = "ng nx ny nxy nz nchan nzth".split()
        data = frombuffer(loc=ct_i32ptrs[0], dim=4*ct_i32lens[0], tp="int32")
        self.sizes = dict(zip(keys,data))
#        for i,key in enumerate(keys):
#            self.sizes[key] = data[i]

        #float64        
        #self.time = frombuffer(loc=ct_f64ptrs[0], dim=8*ct_f64lens[0], tp="float64")       
        self.time = frombuffer(loc=ct_f64ptrs[1], dim=8*ct_f64lens[1], tp="float64")
        self.plev = frombuffer(loc=ct_f64ptrs[2], dim=8*ct_f64lens[2], tp="float64")
        self.flux = frombuffer(loc=ct_f64ptrs[3], dim=8*ct_f64lens[3], tp="float64")
        
        self.relp = frombuffer(loc=ct_f64ptrs[4], dim=8*ct_f64lens[4], tp="float64")
        self.absp = frombuffer(loc=ct_f64ptrs[5], dim=8*ct_f64lens[5], tp="float64")
        self.tcool= frombuffer(loc=ct_f64ptrs[6], dim=8*ct_f64lens[6], tp="float64")
        self.dcool= frombuffer(loc=ct_f64ptrs[7], dim=8*ct_f64lens[7], tp="float64")
        self.tfuel= frombuffer(loc=ct_f64ptrs[8], dim=8*ct_f64lens[8], tp="float64")
        self.boron= frombuffer(loc=ct_f64ptrs[9], dim=8*ct_f64lens[9], tp="float64")
        self.rnxe = frombuffer(loc=ct_f64ptrs[10], dim=8*ct_f64lens[10], tp="float64")
        self.rnio = frombuffer(loc=ct_f64ptrs[11], dim=8*ct_f64lens[11], tp="float64")
        self.rnsm = frombuffer(loc=ct_f64ptrs[12], dim=8*ct_f64lens[12], tp="float64")
        self.rnpm = frombuffer(loc=ct_f64ptrs[13], dim=8*ct_f64lens[13], tp="float64")
        self.keff = frombuffer(loc=ct_f64ptrs[14], dim=8*ct_f64lens[14], tp="float64")


    def dumph5(self, h5file):
        """Setup name of dump file and embed dump function into RENUS
        """
        self.f5 = h5py.File(h5file,"w")
        self.dumph5_data()
        self.set_dumpfunc(self.dumph5_dyn_stage)
        self.f5dyn = self.f5.create_group("dyn")
        self.count = 0        
        self.count_stage = 0

    def dumph5_data(self):
        """Dump basic data of model
        """
        f5dat = self.f5.create_group("data")
        for key in self.sizes:
            f5dat.create_dataset(key, data=[self.sizes[key]])

    def dumph5_dyn_stage(self):
        """
        """
        if self.count_stage == 0:
            self.count_stage = self.stage
            self.dumph5_dyn()
        self.count_stage = self.count_stage - 1        

    def dumph5_dyn(self):
        """Dump data from RENUS during transient calculations
        """
        scope = globals()
        for varname in self.sizes:
            scope[varname] = self.sizes[varname]
   
        f5gr = self.f5dyn.create_group("{0:06d}".format(self.count))
        f5gr.attrs["time,s"] = self.time[0]
        f5gr.attrs["plev,%"] = 100*self.plev[0]
        f5gr.attrs["boron,ppm"] = self.boron[0]
        f5gr.attrs["keff"] = self.keff[0]

        # dump neutron flux by groups
        f5flux = f5gr.create_group("flux")
        for k,flux in enumerate(self.flux.reshape((nz+1,nxy+1,ng)).transpose()):
            f5flux.create_dataset("phif_%i"%k, data=flux[1:].transpose()[1:])
        
        # convert's templates for RENUS's data // because indexation of most arrays in Renus begin from '0'
        conv0 = lambda x: x.reshape((nchan+1,nzth+1))[1:].transpose()[1:]
        conv1 = lambda x: x.reshape((nchan+2,nzth+1))[1:-1].transpose()[1:]
        conv2 = lambda x: x.reshape((nxy,nz)).transpose()

        #Reacor's power
        if(len(self.relp)>0): f5gr.create_dataset("relp", data=conv0(self.relp))
        if(len(self.absp)>0): f5gr.create_dataset("absp", data=conv0(self.absp))

        #Temperatures and densities
        if(len(self.tcool)>0): f5gr.create_dataset("tcool", data=conv1(self.tcool))
        if(len(self.dcool)>0): f5gr.create_dataset("dcool", data=conv1(self.dcool))

        #Poisons concentrations:
        if(len(self.rnxe)>0): f5gr.create_dataset("rnxe", data=conv2(self.rnxe))
        if(len(self.rnio)>0): f5gr.create_dataset("rnio", data=conv2(self.rnio))
        if(len(self.rnsm)>0): f5gr.create_dataset("rnsm", data=conv2(self.rnsm))
        if(len(self.rnpm)>0): f5gr.create_dataset("rnpm", data=conv2(self.rnpm))

        self.count = self.count + 1

    def close(self):
        """close h5 file if open and upload library RENUS
        """
        if hasattr(self,'f5'): self.f5.close()
        self.upload()
