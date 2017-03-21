from os import environ, mkdir, chdir
from os.path import join, abspath, exists
from shutil import copy, copytree, rmtree
from subprocess import call
import h5py
from numpy import array
import time



class sketch:
	"""SKETCH neutron-th code
	"""
	def __init__(self,task,curdir='.'):
		self.tooldir = join(environ["HOME"],"Codes/Sketch")
		self.curdir = abspath(curdir)
		self.outdir = None
		self.task = task
		self.outdir = abspath(task["outdir"])
		self.f5name = join(self.outdir,"sketch.h5")
		self.f5stat = join(self.outdir,"sketch_st.h5")
		self.dt = 0

	def run(self):
		"""Update input files for sketch input directories and run the task
		"""
		task = self.task
		if not exists(self.outdir): mkdir(self.outdir)

		if task.has_key("dat"): copy(task["dat"], join(self.tooldir,"input/sample.dat"))
		if task.has_key("parfh"): copy(task["parfh"], join(self.tooldir,"include/parameters.fh"))
		if task.has_key("history"): copy(task["history"], join(self.tooldir,"input/history.txt"))
		if task.has_key("xslib"): copy(task["xslib"], join(self.tooldir,"input/xslibrary.lib"))

		if task.has_key("st.ini"):
			if task.has_key("st.ini"): copy(task["st.ini"], join(self.tooldir,"input/sketch.ini"))		
			self.execute()
			self.saveh5("sketch_st.h5")

		if task.has_key("kin.ini"):
			if task.has_key("kin.ini"): copy(task["kin.ini"], join(self.tooldir,"input/sketch.ini"))		
			self.execute()
			self.saveh5("sketch.h5")


	def execute(self):
		"""Execute sketch.exe and return to current directory
		"""
		chdir(self.tooldir)
		t1 = time.clock()
		call("sketch.exe")
		self.dt = time.clock() - t1
		chdir(self.curdir)

	def saveh5(self,fname):
		"""Put results into output directory
		"""
		chdir(self.tooldir)
		call("python utils/grf2h5.py output/sketch.grf -o " + join(self.outdir,fname))
		chdir(self.curdir)

	def XeI135(self,f5name=None):
		"""Reading maximum xenon concentration
		"""
		if not f5name: f5name = self.f5name
		f5 = h5py.File(f5name)
		res = []
		for name, g in f5["dyn"].items()[:-1]:
			time = g.attrs["t"]
			xe135 = g['f/Xe135'].value.mean()
			i135  = g['f/i-135'].value.mean()
			res.append([time,xe135,i135])
		f5.close()
		return array(res)


	def offset(self,f5name=None):
		"""Reading offset
		"""
		convert = 24.
		if not f5name: f5name = self.f5name
		f5 = h5py.File(f5name)
		hz=f5['dsk/hz'].value		
		res = []		
		for name, g in f5["dyn"].items()[:-1]:		
			time = g.attrs["t"]
			power = (g['f/POWER, [Wt/cm^3]'].value).sum(1)
			power = (power * hz).reshape(2,len(power)/2)
			top, bottom = power.sum(1)
			offs = 100 * (bottom-top)/(top+bottom)
			res.append([time,offs])
		f5.close()
		return array(res)

	def CR_position(self,f5name=None):
		"""Reading of control rods position
		"""
		if not f5name: f5name = self.f5name
		f5 = h5py.File(f5name)
		res = []
		for name, g in f5["dyn"].items()[:-1]:
			pos = g["v/Control Rod Positions, [cm]"].value
			time = g.attrs["t"]
			res.append([time] + list(pos))
		f5.close()
		return array(res)

	def meshz(self,f5name=None):
		"""Open h5 arhive and get mesh by Oz
		"""
		if not f5name: f5name = self.f5stat
		f5 = h5py.File(f5name)
		addr_hz='dsk/hz'
		hz=f5[addr_hz].value
		f5.close()
		return hz

	def boron_concentration(self,f5name=None):
		"""Open h5 arhive and read boron concentration
		"""
		if not f5name: f5name = self.f5name
		f5 = h5py.File(f5name)
		boron = []
		for name, g in f5["dyn"].items()[:-1]:
			addr_boron='f/BORON CONCENTRATION, [ppm]'
			time = g.attrs["t"]
			boron.append([time,g[addr_boron].value.max()])
		f5.close()
		return array(boron)

	def powert(self,f5name=None):
		"""Open h5 arhive and read boron concentration
		"""
		if not f5name: f5name = self.f5name
		f5 = h5py.File(f5name)
		power = []
		times = []
		for name, g in f5["dyn"].items()[:-1]:
			pars = g["scal_tr"].value
			time = g.attrs["t"]
			power.append([pars[0]])
			times.append([time])
		f5.close()
		return array([times, 100 * array(power) / power[0]]).transpose()

	def power(self,f5name=None):
		"""Open h5 arhive and read power distribution by one channel
		"""
		if not f5name: f5name = self.f5stat
		f5 = h5py.File(f5name)
		addr_pow='dyn/000000/f/POWER, [Wt/cm^3]'
		power = f5[addr_pow].value
		f5.close()
		return power.transpose()[0]

	def boron(self,f5name=None):
		"""Open h5 arhive and read boron concentration
		"""
		if not f5name: f5name = self.f5stat
		f5 = h5py.File(f5name)
		addr_boron='dyn/000000/f/BORON CONCENTRATION, [ppm]'
		boron=f5[addr_boron].value.max()
		f5.close()
		return boron




