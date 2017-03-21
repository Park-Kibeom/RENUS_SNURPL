import h5py

class DumpData:
	def __init__(self,model)


f5 = h5py.File("renus.h5","w")
f5.create_group("dsk")
f5dyn = f5.create_group("dyn")
bank = []
time = []
plev = []
num = 0
def renus_dumph5():
	ng, nx, ny, nxy, nz, nchan, nzth = model.sizes[:7]
	f5gr = f5dyn.create_group("{0:06d}".format(len(time)))
	f5gr.attrs["time,s"] = model.time[0]
	f5gr.attrs["plev,rel"] = model.plev[0]
	f5gr.attrs["boron,ppm"] = model.boron[0]
	# dump neutron flux by groups
	f5flux = f5gr.create_group("flux")
	for k,flux in enumerate(model.flux.reshape((nz+1,nxy+1,ng)).transpose()):
		f5flux.create_dataset("phif_%i"%k, data=flux[1:].transpose()[1:])
	
	# convert's templates for RENUS's data // because indexation of most arrays in Renus begin from '0'
	conv0 = lambda x: x.reshape((nchan+1,nzth+1))[1:].transpose()[1:]
	conv1 = lambda x: x.reshape((nchan+2,nzth+1))[1:-1].transpose()[1:]
	conv2 = lambda x: x.reshape((nxy,nz)).transpose()

	#Reacor's power
	f5gr.create_dataset("relp", data=conv0(model.relp))
	f5gr.create_dataset("absp", data=conv0(model.absp))

	#Temperatures and densities
	f5gr.create_dataset("tcool", data=conv1(model.tcool))
	f5gr.create_dataset("dcool", data=conv1(model.dcool))

	#Poisons concentrations:
	f5gr.create_dataset("rnxe", data=conv2(model.rnxe))
	f5gr.create_dataset("rnio", data=conv2(model.rnio))
#	f5gr.create_dataset("rnsm", data=conv2(model.rnsm))
#	f5gr.create_dataset("rnpm", data=conv2(model.rnpm))


	time.append(model.time[0])
	plev.append(model.plev[0])
	res = model.flux.reshape((19,2,2)).transpose()[0][1][1:]
	bank.append(copy.copy(res)/np.linalg.norm(res))



f5.close()