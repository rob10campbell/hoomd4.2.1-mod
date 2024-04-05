# Copyright (c) 2009-2022 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Implement BoxResize."""

from hoomd.operation import Updater
from hoomd.data.parameterdicts import ParameterDict
from hoomd.variant import Variant, Constant
from hoomd import _hoomd
from hoomd.trigger import Periodic


class BoxShear(Updater):

    def __init__(self, trigger, erate, deltaT, flip):
        params = ParameterDict(erate=Variant,
                               deltaT = float,
                               flip=bool)
        params['erate'] = erate
        params['trigger'] = trigger
        params['deltaT'] = deltaT
        params['flip'] = flip
        self._param_dict.update(params)
        super().__init__(trigger)

    def _attach(self):
        self._cpp_obj = _hoomd.BoxShearUpdater(
            self._simulation.state._cpp_sys_def, self.trigger, self.erate, self.deltaT, self.flip)
        super()._attach()
