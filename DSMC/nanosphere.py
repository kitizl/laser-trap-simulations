#!python3

class gasMolecule:
	def __init__(self, pos, vel):
		self.pos = pos
		self.vel = vel

def generate_velocity():
	return -1

def generate_position():
	return -1

def generate_exit_velocity():
	return -1

def collision(molecule, radius_of_nanosphere):
	
	if isCollide:
		x_col, y_col, z_col = [0,0,0]
		return x_col, y_col, z_col # or whatever floats their boat
	else:
		return "Error"

def eigenvals(matrix):
	# returns the eigenvalues of a square matrix
	vals = (0.0,0.0,0.0)
	return vals

# define a particle geometry (a nanosphere with radius a)

nanoparticle_loc = [0.,0.,0.]
nanoparticle_rad = 1.0 # all rescaled to the size of the nanoparticle
nanoparticle_motion_x = 0.5 # motion of the particle in x-axis
nanoparticle_motion_y = 0.5 # motion of the particle in y-axis
nanoparticle_motion_z = 0.5 # motion of the particle in z-axis
# define all relevant physical constants

factors = 10.0 # insert the clusterfuck of dimensionalizing constants

# define a "radius" of attack

attack_sphere = 10.0 * nanoparticle_rad

# instantiate array of gasmolecules with randomly generated u,v,w

num_of_molecules = 100 # to start with, answers should converge for large values of 100

molecules = [ gasMolecule(generate_position(), generate_velocity()) for i in range(num_of_molecules)]

# check for collisions, find collision point

number_of_collisions = 0

force_x = [0.0,0.0,0.0]

for m in molecules:
	incoming_momentum = m.momentum
	col_coords = collision(m , nanoparticle_rad)
	if isinstance(col_coords, tuple):
		# only if a collision is registered
		number_of_collisions += 1
		# diffuse and reflection
		new_vels = generate_exit_velocity()
		m.vel = new_vels
		# nanosphere, so no need to check for new collisions
		outgoing_momentum = m.momentum
		# change in thing
		d_mom = outgoing_momentum - incoming_momentum

		force_x += factors*d_mom

force_y = [0.0,0.0,0.0]

for m in molecules:
	incoming_momentum = m.momentum
	col_coords = collision(m , nanoparticle_rad)
	if isinstance(col_coords, tuple):
		# only if a collision is registered
		number_of_collisions += 1
		# diffuse and reflection
		new_vels = generate_exit_velocity()
		m.vel = new_vels
		# nanosphere, so no need to check for new collisions
		outgoing_momentum = m.momentum
		# change in thing
		d_mom = outgoing_momentum - incoming_momentum

		force_y += factors*d_mom

force_y = [0.0,0.0,0.0]

for m in molecules:
	incoming_momentum = m.momentum
	col_coords = collision(m , nanoparticle_rad)
	if isinstance(col_coords, tuple):
		# only if a collision is registered
		number_of_collisions += 1
		# diffuse and reflection
		new_vels = generate_exit_velocity()
		m.vel = new_vels
		# nanosphere, so no need to check for new collisions
		outgoing_momentum = m.momentum
		# change in thing
		d_mom = outgoing_momentum - incoming_momentum

		force_y += factors*d_mom

force_tensor = [force_x, force_y, force_z]

gamma_x, gamma_y, gamma_z = eigenvals(force_tensor)

# calculate "damping force"

print(f"Gamma X : {gamma_x/1000} kHz")
print(f"Gamma Y : {gamma_y/1000} kHz")
print(f"Gamma Z : {gamma_z/1000} kHz")

print("\n\nWe expect these values to be as close to each other as possible")