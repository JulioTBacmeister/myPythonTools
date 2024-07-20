from mpi4py import MPI
import numpy as np

def create_halo_region(sub_array, rank, size, halo_size):
    # Add halo regions to the sub-array
    top_halo = np.zeros((halo_size, sub_array.shape[1]))
    bottom_halo = np.zeros((halo_size, sub_array.shape[1]))

    # Send and receive halo regions
    if rank > 0:
        comm.Sendrecv(sub_array[:halo_size, :], dest=rank-1, recvbuf=top_halo)
    if rank < size-1:
        comm.Sendrecv(sub_array[-halo_size:, :], dest=rank+1, recvbuf=bottom_halo)
    
    # Concatenate halos to sub-array
    if rank > 0:
        sub_array = np.vstack((top_halo, sub_array))
    if rank < size-1:
        sub_array = np.vstack((sub_array, bottom_halo))

    return sub_array

def compute_differences(sub_array, halo_size):
    # Compute differences in the interior of the sub-array
    differences = np.diff(sub_array, axis=0)
    return differences[halo_size:-halo_size, :]  # Remove halo regions from result

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Assume a 500x500 array
global_array_shape = (500, 500)
local_row_count = global_array_shape[0] // size
halo_size = 1  # One row of halo on each side

# Initialize the global array on rank 0
if rank == 0:
    global_array = np.random.rand(*global_array_shape)
else:
    global_array = None

# Scatter the global array to all ranks
local_array = np.zeros((local_row_count, global_array_shape[1]))
comm.Scatter(global_array, local_array, root=0)

# Create halo regions
local_array_with_halo = create_halo_region(local_array, rank, size, halo_size)

# Compute differences in the local sub-array
local_differences = compute_differences(local_array_with_halo, halo_size)

# Gather the results back to the root process
if rank == 0:
    global_differences = np.zeros((global_array_shape[0]-1, global_array_shape[1]))
else:
    global_differences = None

comm.Gather(local_differences, global_differences, root=0)

# Print the result on the root process
if rank == 0:
    print(global_differences)
