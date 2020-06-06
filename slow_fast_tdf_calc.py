# Code from the Medium article: Why is Numpy awesome? 
# Simple code for calculating the total distribution function of an arbitary cubic box of atoms.
# This is just a interatomic histogram with no normalisation or density accounted for. 
# Comment/Uncomment slow_method/fast_method below to compare the calculation method 

# Import the relevant modules
import numpy as np
import matplotlib.pyplot as plt
import math
import time

# This function is used to generate a random x,y,z vector
# This is useful for adding some disorder into your box of atoms 
def random_three_vector(npts, radius):
    """
    Taken from: https://gist.github.com/andrewbolster/10274979
    """
    np.random.seed(1)
    phi = np.random.uniform(0,np.pi*2, npts)
    costheta = np.random.uniform(-1,1, npts)

    theta = np.arccos( costheta )
    x = radius*np.sin( theta) * np.cos( phi )
    y = radius*np.sin( theta) * np.sin( phi )
    z = radius*np.cos( theta )

    return np.column_stack((x,y,z))

# ---------------------------------------------------------------------

# Start the clock to check how long the calculation takes
start_code = time.clock()

# Set the fast/slow method for calculting the histogram
# fast_method=True; slow_method=False; plot_title='fast method'
slow_method=True; fast_method=False; plot_title='slow method'

# ---------------------------------------------------------------------
# Input parameters

# Set the size of the cubic box in Angstroms
a_lat = 10.0
box_cell = [a_lat,a_lat,a_lat]

# Set the distance between each atom in a single direction
# Atoms are arranged in a simple 3D grid 
one_d_sep = 1 # Angstroms

# parameter for setting a bit of disorder to the atoms
rand_radius = 0.10

# Set the delta_r value for the PDF
del_r = 0.01
# Set the minimum r-value
min_r = 0.00

# ---------------------------------------------------------------------
# Generating the box of atoms

# Calculate the number of atoms in the box
n_atoms_dir = a_lat/one_d_sep # number of atoms in one direction
n_atoms_tot = int(n_atoms_dir**3)

print('total atoms:', n_atoms_tot)

# Generate the array for the x,y,z coordinates of the atoms in Angstroms
xyz_Ang = np.empty((int(n_atoms_tot),3))
# Assign values to the array
rc = 0
for i in range(int(n_atoms_dir)):
    x = i*one_d_sep
    for j in range(int(n_atoms_dir)):
        y = j*one_d_sep
        for k in range(int(n_atoms_dir)):
            z = k*one_d_sep

            xyz_Ang[rc,0] = x
            xyz_Ang[rc,1] = y
            xyz_Ang[rc,2] = z

            rc += 1

# Add some disorder to the atoms
rand_arr = random_three_vector(n_atoms_tot, rand_radius)
for c in range(3):
    xyz_Ang[:,c] = xyz_Ang[:,c] + rand_arr[:,c]

# See how long that takes
print('time to make box of atoms:', time.clock() - start_code)

# ---------------------------------------------------------------------
# Calculating the histogram

# determine maximum r-value for the rdf
max_r = np.min(box_cell)/2

# max value for the histogram    
max_r_hist = math.ceil(max_r)

# Set the histogram bins
# The following edges mean that the middle of the bins are at del_r values
# remember: range goes from start to end-1
min_bin_edge = min_r + del_r/2
max_bin_edge = max_r_hist + del_r/2 + del_r
npts = int((max_bin_edge - min_bin_edge)/del_r)
bin_edges = np.linspace(min_bin_edge, max_bin_edge, npts, endpoint=False)

# Set the x-values for the plotting the histogram later
r_rdf = bin_edges[0:len(bin_edges)-1] + del_r/2

# Loop through each atom in the box
for i in range(n_atoms_tot):

    # Select out the reference atom
    xi, yi, zi = xyz_Ang[i,0], xyz_Ang[i,1], xyz_Ang[i,2]

    # Fast method for claculating the histogram
    # Makes use of np.historgram
    if fast_method: 

        # Calculate the distance relative to the reference atom i for all atoms
        # Don't need to worry about distance between the atom and itself as the histogram starts from del_r 
        dist_vals = np.sqrt((xyz_Ang[:,0]-xi)**2 + (xyz_Ang[:,1]-yi)**2 + (xyz_Ang[:,2]-zi)**2)

        # Calculate the histogram and add histograms together
        dist_hist_i, _ = np.histogram(dist_vals, bins=bin_edges)
        if i==0: dist_hist = dist_hist_i # Initiate the histogram array
        else: dist_hist += dist_hist_i # Add the extra counts
    
    # Slow method for calculating the histogram
    if slow_method:

        # Create the histogram array
        if i == 0:
            dist_hist = np.empty((len(r_rdf)))
            dist_hist[:] = 0
        
        # Loop through all other atoms
        for j in range(n_atoms_tot):

            # Select out jth atom
            xj, yj, zj = xyz_Ang[j,0], xyz_Ang[j,1], xyz_Ang[j,2]

            # Calculate the distance relative to the reference atom i
            dist_val = np.sqrt((xj-xi)**2 + (yj-yi)**2 + (zj-zi)**2)

            # Bin the distance value
            index = math.floor(dist_val/del_r)
            if 0 < index < len(r_rdf): dist_hist[index] += 1
        

# See how long that takes
print('time to calc histogram:', time.clock() - start_code)

# Plot the histogram
plt.plot(r_rdf, dist_hist)

# Figure aethetics
plt.title(plot_title)
plt.ylabel('no. pairs')
Angst = r'$\AA$'
plt.xlabel('r ('+Angst+')')
plt.tight_layout()

plt.show()
