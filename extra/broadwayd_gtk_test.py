# Dylan Kenneth Eliot

"""
use with broadwayd, and all the gtk3 apps you build within python are browser UI adapted.

Once it is loaded separately, it can then be iframe loaded. 

"""


import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk

win = Gtk.Window(title="Hello World")
win.connect("destroy", Gtk.main_quit)
win.show_all()
Gtk.main()
