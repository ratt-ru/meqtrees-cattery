# -*- coding: utf-8 -*-
#
#% $Id$
#
#
# Copyright (C) 2002-2007
# The MeqTree Foundation &
# ASTRON (Netherlands Foundation for Research in Astronomy)
# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>,
# or write to the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
"""<P>This module implements one solvable instance of the WSRTCos3Beam class from wsrt_cos3_beam. It is preferable to
use the latter module directly -- this module is only retained for backwards compatibility with old scripts.</P>
<P align="right">Author: O. Smirnov &lt;<tt>smirnov@astron.nl</tt>&gt;</P>""";

__default_label__ = "E";
__default_name__  = "WSRT cos^3 beam";

from wsrt_cos3_beam import WSRTCos3Beam

_beam = WSRTCos3Beam("E",solvable=True);

compile_options = _beam.compile_options;
compute_jones = _beam.compute_jones;
compute_jones_tensor = _beam.compute_jones_tensor;
