# Copyright (c) 2013-2022, The MercuryDPM Developers Team. All rights reserved.
# For the list of developers, see <http://www.MercuryDPM.org/Team>.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#   * Neither the name MercuryDPM nor the
#     names of its contributors may be used to endorse or promote products
#     derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE
# ------------This function loads defaults for MClump tool-------------------------

import numpy as np
from stl import mesh

def baseline():

    OPT = {
        'mode' : 4,                # 1 - start with the list of pebbles, compute inertia by summation over pebbles
                                   # 2 - start with the list of pebbles, compute inertia by voxelization
                                   # 3 - start with triangulated surface (stl) file, generate pebbles
                                   # 4 - use both stl and the list of pebbles generated by external library
                                   # 5 - generate stl sequence for blender
        # Input
        'clumpFileName': 'rattleback.txt',         # list of pebbles input file name (mode 1,2,4)
        'stlFileName' : 'rattleback.stl',          # clump surface (mode 3, 4)
        'clumpInputDir' : './input/clump/',     # list of pebbles input dir (mode 1,2,4)
        'stlInputDir' : './input/stl/',         # clump surface input dir (mode 3, 4)
        'clumpSeqDir' : './blender/clump_seq/', # where to search for a clump sequence

        # Options
        'voxNum'  : 100,                        # Voxel grid definition (mode 2)
        'verbose' : True,                       # Detailed messages
        'useColors': True,                      # Use colors in the output
        'useNumba': True,                       # Use just in time pre-compilation with numba
                                                # (available only if numba and llvmlite are installed)
        'rotateToPD' : True,                    # Aligning principal directions with global cartesian axes
        'stlThetaRes': 32,                      # Theta in (0, pi) resolution - stl pebble model
        'stlPhiRes': 32,                        # Phi in (0, 2pi) resolution - stl pebble model

        # Output
        'clumpName':            'clump04',                             # Used in output dirs
        'paraviewOutput':       True,                            # Save paraview files for clump visualization
        'VoxelStlOutput':       True,                            # Save voxels as stl
        'PebbleStlOutput':      True,                            # Save pebbles as stl
        'animate':              True,  # make animation of the rotating clump, otherwise - static
        'numTimesteps':         100,  # Number of timesteps in the animation - 'animate' = True
        'clumpOutputDir' :      'clump',                         # Directory with the output list of pebbles
        'inertiaOutputDir' :    'inertia',                     # Mass, toi, principal directions
        'paraviewOutputDir' :   'vtu',                        # vtu, pebbles + principal axes
        'voxelStlOutputDir':    'voxelgrid',  # stl (voxels)
        'stlOutputDir':         'stl',                        # stl (voxels)
        'vtuFileName':          'body',                              # vtu output file name
        'voxStlFileName':       'voxels.stl',                     # voxels stl output file name
        'outStlFileName':       'model.stl',  # voxels stl output file name
        'outClumpFileName':     'clump.txt',                    # output clump file name
        'inertiaFileName':      'toi.txt',                       # tensor of inertia output file name
        'pdFileName':           'pd.txt',                             # paraview output file name
        'massFileName':         'mass.txt',                         # paraview output file name
        'stlSeqDir':            './blender/stl_seq',
        # Clump generation options
        'meshToClumpAlg':   1,                              # Clump reneration algorithm (1 - regular placement)
        'regPlacementDefinition': 10                        # Number of pebbles along the side of the bounding box (regular placement algorithm)

    }

    DATA = {
        # Model data
        'density': 1,       # Default density of the clump material
        'stlMesh': None,    # stl mesh
        'toi': None,        # Tensor of inertia
        'pd': None,         # Principal directions TODO: PD should form the right basis - otherwise - switch PDs
        'mass': 0,          # Mass of the clump
        'pebbles': None,    # List of pebbles
        'voxelGrid': None,   # Voxel grid resolving clump (see OPT['voxNum'] for voxel grid resolution)
        'clumpSequence': None
    }

    return OPT, DATA
