import re
from os import environ, mkdir, chdir
from os.path import join, abspath, exists, isfile
from shutil import copy, copytree, rmtree
from subprocess import call
from numpy import array, loadtxt
import time

class renus:
	"""RENUS code
	"""
	def __init__(self, task, curdir='.',mode="release"):
		self.mode = mode
		self.curdir = abspath(curdir)
		self.tooldir = join(environ["HOME"],"Codes/Renus")
		self.outdir = task["outdir"]
		self.inpdir = task["inpdir"]
		self.dt = 0

	def run(self):
		if exists(join(self.tooldir,"work")): rmtree(join(self.tooldir,"work"))
		copytree(self.inpdir, join(self.tooldir,"work"))

		self.execute()

		if not exists(self.outdir): mkdir(self.outdir)

		f1 = join(self.tooldir,"work/out/sample.out")
		f2 = join(self.tooldir,"work/out/sample.plt")
		if isfile(f1):
			copy(f1, join(self.outdir,"sample.out"))
		if isfile(f2):
			copy(f2, join(self.outdir,"sample.plt"))

	def execute(self):
		chdir(join(self.tooldir, "work"))
		t1 = time.clock()
		call("../bin/renus.exe sample.inp") #%self.mode)
		self.dt = time.clock() - t1
		chdir(self.curdir)

	def offset(self):
		convert = 1./3600
		out = loadtxt(join(self.outdir,"sample.plt"), skiprows=1).transpose()
		time, value = out[0], 100*out[-1]
		return array([time,value]).transpose()

	def powert(self):
		out = loadtxt(join(self.outdir,"sample.plt"), skiprows=1).transpose()
		time, value = out[0], out[1]
		return array([time,value]).transpose()

	def tcoolt(self):
		out = loadtxt(join(self.outdir,"sample.plt"), skiprows=1).transpose()
		time, value = out[0], out[4]
		return array([time,value]).transpose()


	def power(self):
		"""Open output file of RENUS and read power distribution by one channel
		"""
		power = []
		f = open(join(self.outdir,"sample.out"),"r")
		for line in f:
			if re.search(r"Core Average Axial Power Distribution", line):
				for line in f:
					if line == "\n":
						break
					else:
						power.append(line.split()[1])
				break
		power = map(float,power)
		f.close()
		return array(power)

	def meshz(self):
		"""Open output file of RENUS and get mesh by Oz
		"""
		f = open(join(self.outdir,"sample.out"),"r")
		hz = []
		for line in f:
			if re.search(r"grid_z", line):
				data = line.split()[1:]
				for d in data:
					a = d.split('*')
					if len(a) == 2:
						hz.extend(int(a[0])*[a[1]])
					else:
						hz.append(a[0])
				break
		f.close()
		hz = map(float,hz)
		return hz

	def boron(self):
		"""Open output file of RENUS and read boron concentration
		"""
		f = open(join(self.outdir,"sample.out"),"r")
		for line in f:
			if re.search(r"BORON  ", line):
				boron = float(line.split(":")[-1])
		return boron