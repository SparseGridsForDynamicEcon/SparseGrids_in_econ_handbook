import math
import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def pad(arrayToPad):
    """Takes an input array and pads it across the first axis such that it is 
    divisible by the number of points per worker.

    Keyword arguments:
    arrayToPad -- numpy array of inputs to pad
    """
    pointsPerWorker = math.ceil(arrayToPad.shape[0] / comm.Get_size())

    return np.pad(
        arrayToPad,
        (
            (
                0,
                pointsPerWorker * comm.Get_size() - arrayToPad.shape[0],
            ),
            (0, 0),
        ),
        mode="constant",
    )

def scatter(arrayToScatter):
    """Takes an input array and scatters it across workers.

    Keyword arguments:
    arrayToScatter -- numpy array of inputs to scatter
    """
    pointsPerWorker = math.ceil(arrayToScatter.shape[0] / comm.Get_size())
    sendbuf = pad(arrayToScatter)
    recvbuf = np.empty((pointsPerWorker, arrayToScatter.shape[1]))
    comm.Scatter(sendbuf, recvbuf, root=0)

    return recvbuf