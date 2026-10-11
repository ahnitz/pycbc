# Copyright (C) 2026 PyCBC Team
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.

"""Legacy alias module for LineContinuumPSDModel."""

import warnings
from .line_continuum import LineContinuumPSDModel


class BayesLinePSDModel(LineContinuumPSDModel):
    """Deprecated alias for LineContinuumPSDModel.
    Please use pycbc.psd.estimators.LineContinuumPSDModel instead.
    """
    def __init__(self, *args, **kwargs):
        warnings.warn(
            "BayesLinePSDModel is deprecated and will be removed in a future release. "
            "Use LineContinuumPSDModel instead.",
            FutureWarning,
            stacklevel=2,
        )
        super().__init__(*args, **kwargs)
