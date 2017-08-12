#
#
#    This file is part of APSYNSIM: A real-time Aperture Synthesis Simulator
#
#    Copyright (C) 2014  Ivan Marti-Vidal (Nordic ARC Node, OSO, Sweden)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

# TODO: atm produces a number of deprecation issues around non-int indices
# from __future__ import (absolute_import, division, print_function)
from __future__ import (absolute_import, print_function)

try:
    import Tkinter as Tk
except:
    import tkinter as Tk

import interferometer

__version__ = '2.1-a'

if __name__ == "__main__":

    root = Tk.Tk()
    TITLE = ('Aperture Synthesis Simulator - ESO Supernova fork '
             '(based on software by I. Marti-Vidal, '
             'Onsala Space Observatory) - version  %s' % __version__)
    root.wm_title(TITLE)

    myint = interferometer.Interferometer(tkroot=root)
    Tk.mainloop()
