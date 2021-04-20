#!python3

def make_dir(dir_name):
	"""
	A function that creates a directory at the path dir_name, if it doesn't exist.
		dir_name : path of the desired directory(ies)
	"""
	path = dir_name
	try:
		os.makedirs(path)
	except FileExistsError:
		print(f"Directory {path} already exists")
	except OSError:
		print (f"Creation of the directory %s failed" % path)
	else:
		print ("Successfully created the directory %s " % path)

def npy2csv(data_loc):
	"""
		A helper function that converts .npy to .csv 
			data_loc : the path pointing to .npy file
	"""
	import numpy as np
	# loads the .npy file and saves the 2D array into a .csv file
	np.savetxt(f"{data_loc[:-4]}.csv", np.load(data_loc), delimiter=",")
	# since this function is source-agnostic, there's no real way to tell 
	# what the headers (or labels for the features) of the datset is
