#
# LSST Data Management System
# Copyright 2016 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import astropy.io.fits as pf
import sys

# '15A field': ('14A field', visit chosen as template)
FIELD_PAIRS = {'Blind15A_26': ('Blind14A_04', 289444),
               'Blind15A_42': ('Blind14A_09', 289449),
               'Blind15A_40': ('Blind14A_10', 289450)}

templateGroups = dict([(x[1], (list(), x[0])) for x in FIELD_PAIRS.values()])


def main():
    """Process the given HiTS2015 files on the command line and
    collect their visit numbers with their corresponding 2014
    template visit and field name. Supported fields are:
    ``Blind15A_26``, ``Blind15A_42``, ``Blind15A_40``

    Output to be inserted as ``--id visit=``, ``--templateId visit=``
    and ``--templateId object=`` fields for ``ap_pipe.py``.
    """
    for F in sys.argv[1:]:
        with pf.open(F) as hdul:

            x = hdul[0].header['OBJECT']
            if x not in FIELD_PAIRS.keys():
                raise ValueError("Unknown field")

            v15 = hdul[0].header['EXPNUM']
            templateGroups[FIELD_PAIRS[x][1]][0].append(str(v15))

    for x in templateGroups.items():
        print("^".join(x[1][0]), "   ", x[0], "  ", x[1][1])


if __name__ == "__main__":
    main()
