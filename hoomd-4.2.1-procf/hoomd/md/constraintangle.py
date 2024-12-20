# Copyright (c) 2009-2022 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

r"""Angle forces.

Angle force classes apply a force and virial on every particle in the simulation
state commensurate with the potential energy:

.. math::

    U_\mathrm{angle} = \sum_{(i,j,k) \in \mathrm{angles}} U_{ijk}(\theta)

Each angle is defined by an ordered triplet of particle tags in the
`hoomd.State` member ``angle_group``. HOOMD-blue does not construct angle
groups, users must explicitly define angles in the initial condition.

.. image:: md-angle.svg
    :alt: Definition of the angle bond between particles i, j, and k.

In the angle group (i,j,k), :math:`\theta` is the angle between the vectors
:math:`\vec{r}_{ij}` and :math:`\vec{r}_{kj}`.

.. rubric Per-particle energies and virials

Angle force classes assign 1/3 of the potential energy to each of the particles
in the angle group:

.. math::

    U_l = \frac{1}{3} \sum_{(i,j,k) \in \mathrm{angles}}
    U_{ijk}(\theta) [l=i \lor l=j \lor l=k]

and similarly for virials.
"""

from hoomd.md import _md
from hoomd.md.force import Force
from hoomd.data.parameterdicts import ParameterDict, TypeParameterDict
from hoomd.data.typeparam import TypeParameter
import hoomd
import numpy

class constangle(Force):
    """Base class for angle forces.

    Note:
        :py:class:`ConstAngle` is the base class for all angle forces. Users should
        not instantiate this class directly.
    """

    def __init__(self, nlist):
        super().__init__()
        self._nlist = nlist  # Directly store the NeighborList object
        params = TypeParameter(
            'params', 'angle_types',
            TypeParameterDict(t0=float, k=float, len_keys=1)
        )
        self._add_typeparam(params)

    def _attach_hook(self):
        # Ensure the neighbor list is part of the simulation
        if not self._nlist._attached:
            self._nlist._attach(self._simulation)

        # Select the appropriate compute class based on the device
        cls = _md.ConstraintAngleForceCompute

        # Create the C++ object
        self._cpp_obj = cls(self._simulation.state._cpp_sys_def,
                            self._nlist._cpp_obj)
        super()._attach_hook()

    @property
    def nlist(self):
        """Neighbor list used to compute the real-space term."""
        return self._nlist

    @nlist.setter
    def nlist(self, value):
        if self._attached:
            raise RuntimeError("nlist cannot be set after scheduling.")
        else:
            self._nlist = value  # Directly assign the NeighborList object

    @property
    def _children(self):
        """Return the child objects (for the simulation tree)."""
        return [self.nlist]

