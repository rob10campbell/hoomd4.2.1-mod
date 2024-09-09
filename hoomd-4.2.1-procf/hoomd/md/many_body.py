# Copyright (c) 2009-2023 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

########## Modified by Rheoinformatic //~ [RHEOINF] ##########

r"""Implement many body potentials.

Triplet force classes apply a force and virial on every particle in the
simulation state commensurate with the potential energy:

.. math::

    U_\mathrm{many-body} = \frac{1}{2} \sum_{i=0}^\mathrm{N_particles-1}
                           \sum_{j \ne i}
                           \sum_{j \ne k} U(\vec{r}_{ij}, \vec{r}_{ik})

where :math:`\vec{r}_{ij} = \mathrm{minimum\_image}(\vec{r}_j - \vec{r}_i)`.
`Triplet` applies a short range cutoff for performance and assumes that both
:math:`U(\vec{r}_{ij}, \vec{r}_{ik})` and its derivatives are 0 when
:math:`r_{ij} > r_\mathrm{cut}` or :math:`r_{ik} > r_\mathrm{cut}`.

Specifically, the force :math:`\vec{F}` applied to each particle :math:`i` is:

.. math::
    \vec{F_i} =
    \begin{cases}
     -\nabla V(\vec r_{ij}, \vec r_{ik})
            & r_{ij} < r_{\mathrm{cut}}
            \land r_{ik} < r_{\mathrm{cut}} \\
        0 & \mathrm{otherwise}
    \end{cases}

The per particle energy terms are:

.. math::

    U_i = \frac{1}{2} \sum_{j \ne i}
          \sum_{j \ne k} U(\vec{r}_{ij}, \vec{r}_{ik})
          [r_{ij} < r_{\mathrm{cut}} \land r_{ik} < r_{\mathrm{cut}}]
"""

import copy
import warnings

import hoomd
from hoomd.data.parameterdicts import ParameterDict, TypeParameterDict
from hoomd.data.typeconverter import positive_real
from hoomd.data.typeparam import TypeParameter
from hoomd.md import _md
from hoomd.md.force import Force


class Triplet(Force):
    """Base class triplet force.

    `Triplet` is the base class for many-body triplet forces.

    Warning:
        This class should not be instantiated by users. The class can be used
        for `isinstance` or `issubclass` checks.

    .. py:attribute:: r_cut

        *r_cut* :math:`[\\mathrm{length}]`, *optional*: defaults to the value
        ``default_r_cut`` specified on construction.

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `float`])

    .. py:attribute:: nlist

        Neighbor list used to compute the triplet potential.

        Type: `hoomd.md.nlist.NeighborList`

    Warning:
        Currently HOOMD-blue does not support reverse force communication
        between MPI domains on the GPU. Since reverse force communication is
        required for the calculation of three-body forces, attempting to use
        this potential on the GPU with MPI will result in an error.
    """

    def __init__(self, nlist, default_r_cut=None):
        super().__init__()
        r_cut_param = TypeParameter(
            'r_cut', 'particle_types',
            TypeParameterDict(positive_real, len_keys=2))
        if default_r_cut is not None:
            r_cut_param.default = default_r_cut
        self._add_typeparam(r_cut_param)
        self._param_dict.update(
            ParameterDict(nlist=hoomd.md.nlist.NeighborList))
        self.nlist = nlist

    def _setattr_param(self, attr, value):
        if attr == "nlist":
            self._nlist_setter(value)
            return
        super()._setattr_param(attr, value)

    def _nlist_setter(self, new_nlist):
        if self._attached:
            raise RuntimeError("nlist cannot be set after scheduling.")
        self._param_dict._dict["nlist"] = new_nlist

    def _attach_hook(self):
        if self.nlist._attached and self.nlist._simulation != self._simulation:
            warnings.warn(
                f"{self} object is creating a new equivalent neighbor list."
                f" This is happending since the force is moving to a new "
                f"simulation. Explicitly set the nlist to hide this warning.",
                RuntimeWarning)
            self.nlist = copy.deepcopy(self.nlist)
        self.nlist._attach(self._simulation)
        self.nlist._cpp_obj.setStorageMode(_md.NeighborList.storageMode.full)
        if isinstance(self._simulation.device, hoomd.device.CPU):
            cls = getattr(_md, self._cpp_class_name)
        else:
            cls = getattr(_md, self._cpp_class_name + "GPU")

        self._cpp_obj = cls(self._simulation.state._cpp_sys_def,
                            self.nlist._cpp_obj)


class Tersoff(Triplet):
    r"""Tersoff force.

    Args:
        nlist (hoomd.md.nlist.NeighborList): Neighbor list
        default_r_cut (float): Default cutoff radius :math:`[\mathrm{length}]`.

    The Tersoff potential is a bond-order potential based on the Morse potential
    that accounts for the weakening of individual bonds with increasing
    coordination number. It does this by computing a modifier to the attractive
    term of the potential. The modifier contains the effects of third-bodies on
    the bond energies. The potential also includes a smoothing function around
    the cutoff. The implemented smoothing function is exponential in nature as
    opposed to the sinusoid used by `J. Tersoff 1988`_.

    .. _J. Tersoff 1988:
      https://journals.aps.org/prb/abstract/10.1103/PhysRevB.38.9902

    `Tersoff` computes the Tersoff three-body force on every particle in the
    simulation state. Despite the fact that the Tersoff potential accounts for
    the effects of third bodies, the species of the third body is irrelevant. It
    can thus use type-pair parameters similar to those of the pair potentials:

    .. math::

        U_{ij}(r) = \frac{1}{2} f_C(r_{ij})
            \left[f_R(r_{ij}) + b_{ij}f_A(r_{ij})\right]

    where

    .. math::

        f_R(r) = A_1e^{\lambda_1(r_D-r)}

        f_A(r) = A_2e^{\lambda_2(r_D-r)}

    .. math::
        f_C(r) =
        \begin{cases}
            1 & r < r_{\mathrm{cut}} - r_{CT} \\
            \exp \left[-\alpha\frac{x(r)^3}{x(r)^3 - 1} \right]
                & r_{\mathrm{cut}} - r_{CT} < r < r_{\mathrm{cut}} \\
            0 & r > r_{\mathrm{cut}}
        \end{cases}

    .. math::

        b_{ij} = (1 + \gamma^n\chi_{ij}^n)^{\frac{-1}{2n}}

    In the definition of :math:`f_C(r)`, there is a quantity :math:`x(r)`, which
    is defined as

    .. math::

        x(r) = \frac{r - (r_{\mathrm{cut}} - r_{CT})}{r_{CT}}

    which ensures continuity between the different regions of the potential. In
    the definition of :math:`b_{ij}`, there is a quantity :math:`\chi_{ij}`
    which is defined as

    .. math::

        \chi_{ij} = \sum_{k \neq i,j} f_C(r_{ik})
                     \cdot e^{\lambda_3^3 |r_{ij} - r_{ik}|^3}
        \cdot g(\theta_{ijk})

        g(\theta_{ijk}) = 1 + \frac{c^2}{d^2}
                           - \frac{c^2}{d^2 + |m - \cos(\theta_{ijk})|^2}

    .. py:attribute:: params

        The Tersoff potential parameters. The dictionary has the following
        keys:

        * ``magnitudes`` (tuple[`float`, `float`]) - :math:`(A_1, A_2)` -
          Magnitudes of the repulsive and attractive
          terms (*default*: (1.0, 1.0)) :math:`[\mathrm{energy}]`
        * ``exp_factors`` (tuple[`float`, `float`]) -
          :math:`(\lambda_1, \lambda_2)` - exponential factors of the
          repulsive and attractive
          terms (*default*: 2.0) :math:`[\mathrm{length}^{-1}]`
        * ``lambda3`` (`float`) - :math:`\lambda_3` - exponential factor in
          :math:`\chi_{ij}` (*default*: 0.0) :math:`[\mathrm{length}^{-1}]`
        * ``dimer_r`` (`float`) - :math:`r_D` - length shift in attractive
          and repulsive terms (*default*: 1.5) :math:`[\mathrm{length}]`
        * ``cutoff_thickness`` (`float`) - :math:`r_{CT}` - distance which
          defines the different regions of the
          potential (*default*: 0.2) :math:`[\mathrm{length}]`
        * ``alpha`` (`float`) - :math:`\alpha` - decay rate of the cutoff
          term
          :math:`f_C(r)` (*default*: 3.0) :math:`[\mathrm{dimensionless}]`
        * ``n`` (`float`) - :math:`n` - power in
          :math:`b_{ij}` (*default*: 0.0) :math:`[\mathrm{dimensionless}]`
        * ``gamma`` (`float`) - :math:`\gamma` - coefficient in
          :math:`b_{ij}` (*default*: 0.0) :math:`[\mathrm{dimensionless}]`
        * ``c`` (`float`) - :math:`c` - coefficient in
          :math:`g(\theta)` (*default*: 0.0) :math:`[\mathrm{dimensionless}]`
        * ``d`` (`float`) - :math:`d` - coefficient in
          :math:`g(\theta)` (*default*: 1.0) :math:`[\mathrm{dimensionless}]`
        * ``m`` (`float`) - :math:`m` - coefficient in
          :math:`g(\theta)` (*default*: 0.0) :math:`[\mathrm{dimensionless}]`

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = md.nlist.Cell()
        tersoff = md.many_body.Tersoff(default_r_cut=1.3, nlist=nl)
        tersoff.params[('A', 'B')] = dict(magnitudes=(2.0, 1.0), lambda3=5.0)
    """
    _cpp_class_name = "PotentialTersoff"

    def __init__(self, nlist, default_r_cut=None):
        super().__init__(nlist, default_r_cut)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(cutoff_thickness=0.2,
                              magnitudes=(1.0, 1.0),
                              exp_factors=(2.0, 1.0),
                              lambda3=0.0,
                              dimer_r=1.5,
                              n=0.0,
                              gamma=0.0,
                              c=0.0,
                              d=1.0,
                              m=0.0,
                              alpha=3.0,
                              len_keys=2))
        self._add_typeparam(params)


class RevCross(Triplet):
    r"""Reversible crosslinker three-body force.

    Args:
        nlist (hoomd.md.nlist.NeighborList): Neighbor list
        default_r_cut (float): Default cutoff radius :math:`[\mathrm{length}]`.

    `RevCross` computes the revcross three-body force on every particle in the
    simulation state. Despite the fact that the revcross potential accounts for
    the effects of third bodies, it is actually just a combination of two body
    potential terms. It can thus use type-pair parameters similar to those of
    the pair potentials.

    The RevCross potential has been described in detail in
    `S. Ciarella and W.G. Ellenbroek 2019 <https://arxiv.org/abs/1912.08569>`_.
    It is based on a generalized Lennard-Jones pairwise attraction to form bonds
    between interacting particles:

    .. math::
        U_{ij}(r) =
        \begin{cases}
        4 \varepsilon \left[
            \left( \dfrac{ \sigma}{r_{ij}} \right)^{2n}
            - \left( \dfrac{ \sigma}{r_{ij}} \right)^{n}
        \right] & r < r_\mathrm{cut} \\
        0 & r \ge r_\mathrm{cut}
        \end{cases}

    Then an additional three-body repulsion is evaluated to compensate the bond
    energies imposing single bond per particle condition:

    .. math::

        v^{\left( 3b \right)}_{ijk} = \lambda \epsilon
        \hat{v}^{ \left( 2b \right)}_{ij}
        \left(\vec{r}_{ij}\right)
        \cdot \hat{v}^{ \left( 2b \right)}_{ik}
        \left(\vec{r}_{ik}\right)~,

    where the two body potential is rewritten as:

    .. math::

        \hat{v}^{ \left( 2b \right)}_{ij}\left(\vec{r}_{ij}\right) =
        \begin{cases}
        1 & r \le r_{min} \\
        - \dfrac{v_{ij}\left(\vec{r}_{ij}\right)}{\epsilon}
            & r > r_{min} \\
        \end{cases}

    .. attention::

        The RevCross potential models an asymmetric interaction between two
        different chemical moieties that can form a reversible bond. This
        requires the definition of (at least) two different types of particles.
        A reversible bond is only possible between two different species,
        otherwise :math:`v^{\left( 3b \right)}_{ijk}`, would prevent any bond.
        In our example we then set the interactions for types A and B with
        ``rev_c.params[[('A','B'),('A','B')]]`` to
        ``{"sigma": 0.0, "n": 0, "epsilon": 0, "lambda3": 0}``
        and the only non-zero energy only between the different types with
        setting ``rev_c.params[('A','B')]`` to
        ``{"sigma":1, "n": 100, "epsilon": 100, "lambda3": 1}``.
        Notice that the number of the minority species corresponds to the
        maximum number of bonds.


    This three-body term also tunes the energy required for a bond swap through
    the coefficient:- :math:`\lambda` - *lambda3* (unitless)
    in `S. Ciarella and W.G. Ellenbroek 2019
    <https://arxiv.org/abs/1912.08569>`__
    is explained that setting :math:`\lambda=1` corresponds to no energy
    requirement to initiate bond swap, while this energy barrier scales roughly
    as :math:`\beta \Delta E_\text{sw} =\beta \varepsilon(\lambda-1)`.

    Note:
        Choosing :math:`\lambda=1` pushes the system to cluster because the
        three-body term is not enough to compensate the energy of multiple
        bonds, so it may cause nonphysical situations.

    .. py:attribute:: params

        The revcross potential parameters. The dictionary has the following
        keys:

        * ``epsilon`` (`float`, **required**) - :math:`\varepsilon`
          :math:`[\mathrm{energy}]`

        * ``sigma`` (`float`, **required**) - :math:`\sigma`
          :math:`[\mathrm{length}]`

        * ``n`` (`float`, **required**) - :math:`n`
          :math:`[\mathrm{dimensionless}]`

        * ``lambda3`` (`float`, **required**) - :math:`\lambda_3`
          :math:`[\mathrm{dimensionless}]`

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = md.nlist.Cell()
        bond_swap = md.many_body.RevCross(default_r_cut=1.3,nlist=nl)
        bond_swap.params[('A', 'A'), ('B', 'B')] = {
            "sigma":0,"n": 0, "epsilon": 0, "lambda3": 0}
        # a bond can be made only between A-B and not A-A or B-B
        bond_swap.params[('A','B')] = {
            "sigma": 1, "n": 100, "epsilon": 10, "lambda3": 1}
    """
    _cpp_class_name = "PotentialRevCross"

    def __init__(self, nlist, default_r_cut=None):
        super().__init__(nlist, default_r_cut)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(sigma=2.0,
                              n=1.0,
                              epsilon=1.0,
                              lambda3=1.0,
                              len_keys=2))
        self._add_typeparam(params)


class SquareDensity(Triplet):
    r"""Soft force for simulating a van der Waals liquid.

    Args:
        nlist (hoomd.md.nlist.NeighborList): Neighbor list
        default_r_cut (float): Default cutoff radius :math:`[\mathrm{length}]`.

    `SquareDensity` that the square density three-body force should on every
    particle in the simulation state.

    The self energy per particle takes the form:

    .. math::

        \Psi^{ex} = B (\rho - A)^2

    which gives a pair-wise additive, three-body force:

    .. math::

        \vec{f}_{ij} = \left( B (n_i - A) + B (n_j - A) \right) w'_{ij}
        \vec{e}_{ij}

    Here, :math:`w_{ij}` is a quadratic, normalized weighting function,

    .. math::

        w(x) = \frac{15}{2 \pi r_{c,\mathrm{weight}}^3}
               (1-r/r_{c,\mathrm{weight}})^2

    The local density at the location of particle *i* is defined as

    .. math::

        n_i = \sum\limits_{j\neq i} w_{ij}
              \left(\big| \vec r_i - \vec r_j \big|\right)

    .. py:attribute:: params

        The SquareDensity potential parameters. The dictionary has the
        following keys:

        * ``A`` (`float`, **required**) - :math:`A` - mean density
          (*default*:0) :math:`[\mathrm{length}^{-2}]` in 2D and
          :math:`[\mathrm{length}^{-3}]` in 3D
        * ``B`` (`float`, **required**) - :math:`B` - coefficient of the
          harmonic density term
          :math:`[\mathrm{energy} \cdot \mathrm{length}^4]` in 2D and
          :math:`[\mathrm{energy} \cdot \mathrm{length}^6]` in 3D

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = nlist.Cell()
        sqd = md.many_body.SquareDensity(nl, default_r_cut=3.0)
        sqd.params[('A', 'B')] = dict(A=1.0, B=2.0)
        sqd.params[('B', 'B')] = dict(A=2.0, B=2.0, default_r_on=1.0)

    For further details regarding this multibody potential, see

    [1] P. B. Warren, "Vapor-liquid coexistence in many-body dissipative
    particle dynamics" Phys. Rev. E. Stat. Nonlin. Soft Matter Phys., vol. 68,
    no. 6 Pt 2, p. 066702, 2003.
    """
    _cpp_class_name = "PotentialSquareDensity"

    def __init__(self, nlist, default_r_cut=None):
        super().__init__(nlist, default_r_cut)
        params = TypeParameter('params', 'particle_types',
                               TypeParameterDict(A=0.0, B=float, len_keys=2))
        self._add_typeparam(params)


##~ add AngularRepulsion [RHEOINF]
class AngularRepulsion(Triplet):
    r"""Angular repulsion three-body force.

    Args:
        nlist (hoomd.md.nlist.NeighborList): Neighbor list
        default_r_cut (float): Default cutoff radius :math:`[\mathrm{length}]`.

    `AngularRepulsion` computes the angle-based three-body force on every particle in the
    simulation state. Despite the fact that this potential accounts for
    the effects of third bodies, it is actually just a combination of two body
    potential terms. It can thus use type-pair parameters similar to those of
    the pair potentials.

    The AngularRepulsion potential is based on the Stillinger-Weber (SW) potential, and
    was  developed by Minaspi Bantawa, Wayan Fontaine-Seiler, Peter Olmsted, and 
    Emanuela Del Gado for coarse-grained gel structure analysis. It is described 
    in M. Bantawa et.al. 2021 <https://doi.org/10.1088/1361-648X/ac14f6>

    It has been modified to include a 2-body Morse potential:

    .. math::
        U_{ij},(r_{ij}) =
        \begin{cases}
        D_0 \left[ \exp \left(-2\alpha\left(
            r_ij - \left( a_i + a_j \right) - r_0 \right)\right) 
            - 2 \exp \left(-\alpha\left(r_ij -
            \left( a_i + a_j \right) - r_0 \right)
            \right) \right]
        \end{cases}

    Then an additional three-body angular repulsion is evaluated:

    .. math::

        U_{ij},(\vec{r}_{ij},\vec{r}_{ik}) =
        \begin{cases}
        B \Lambda \left( r \right) \Lambda \left( r' \right) \exp \left[
            - \left( \frac{{r}_{ij} \dot \vec{r}_{ik}}{r_{ij}r_{ik}}
            - \cos \bar{\theta} \right)^2 / w^2 \right]
        \end{cases}


    .. py:attribute:: params

        The potential parameters. The dictionary has the following
        keys:

        * ``D0`` (`float`, **required**) - depth of the potential at its
          minimum :math:`D_0` :math:`[\mathrm{energy}]`
        * ``alpha`` (`float`, **required**) - the width of the potential well
          :math:`\alpha` :math:`[\mathrm{length}^{-1}]`
        * ``r0`` (`float`, **required**) - position of the minimum
          :math:`r_0` :math:`[\mathrm{length}]`
        * ``scaled_D0`` (bool) - on/off class attribute for scaling D0 by particle size (D0*((radius_i_radius_j)/2); defaults to False [RHEOINF]
          :math: `scaled_D0` :math: `[true/false]` [RHEOINF]

        * ``theta_bar`` (`float`, **required**) - size of the reference angle limit (in radians)
          :math:`theta_bar` :math:`[\mathrm{radians}]`
        * ``B`` (`float`, **required**) - magnitude of angular repulsion 
          :math:`B` :math:`[\mathrm{energy}]`
        * ``w`` (`float`, **required**) - width of the angular repulsion potential (about peak B at position theta_bar)
          :math:`w` :math:`[\mathrm{length?}]`

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = md.nlist.Cell()
        angular = md.many_body.AngularRepulsion(default_r_cut=3.0,nlist=nl)
        angular.params[('A', 'A'), ('B', 'B')] = {
            "D0": 1.0,"alpha": 3.0, "r0": 1.0, "theta_bar":1.13446, "B":3, "w":0.30, "scaled_D0":false}
        angular.params[('A','B')] = {
            "D0": 1.0,"alpha": 3.0, "r0": 1.0, "theta_bar":1.13446, "B":1, "w":5, "scaled_D0":false}

    """
    _cpp_class_name = "PotentialAngularRepulsion"

    def __init__(self, nlist, default_r_cut=None):
        super().__init__(nlist, default_r_cut)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(D0=float,
                              alpha=float,
                              r0=float,
                              theta_bar=float,
                              B=float,
                              w=float,
                              scaled_D0=bool,
                              len_keys=2))
        self._add_typeparam(params)
##~
