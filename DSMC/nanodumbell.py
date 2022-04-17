#!python3

class gasMolecule:
	def __init__(self, pos, vel):
		self.pos = pos
		self.vel = vel


# define a particle geometry (a nanodumbell with radius a)

nanoparticle_loc = [0.,0.,0.]
nanoparticle_rad = 1.0 # all rescaled to the size of the nanoparticle

# define all relevant physical constants

# define a "radius" of attack

attack_sphere = 10.0 * nanoparticle_rad

# generate random numbers from correct CDF for incoming velocities

# check for collisions, find collision point

# keep checking for collisions for each particle until you cannot anymore

# find the exit momentum

# calculate "damping force"