import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
import gspread
import datetime

class LoginWindow(Gtk.Window):
    def __init__(self):
        Gtk.Window.__init__(self, title="Login")

        self.set_border_width(10)

        # Create a vertical box to hold the widgets
        vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
        self.add(vbox)

        # Create a label for the username
        username_label = Gtk.Label("Username")
        vbox.pack_start(username_label, True, True, 0)

        # Create an entry for the username
        self.username_entry = Gtk.Entry()
        vbox.pack_start(self.username_entry, True, True, 0)

        # Create a label for the password
        password_label = Gtk.Label("Password")
        vbox.pack_start(password_label, True, True, 0)

        # Create an entry for the password
        self.password_entry = Gtk.Entry()
        self.password_entry.set_visibility(False)
        vbox.pack_start(self.password_entry, True, True, 0)

        # Create a login button
        login_button = Gtk.Button(label="Login")
        login_button.connect("clicked", self.on_login_clicked)
        vbox.pack_start(login_button, True, True, 0)

        # Create a logout button
        logout_button = Gtk.Button(label="Logout")
        logout_button.connect("clicked", self.on_logout_clicked)
        vbox.pack_start(logout_button, True, True, 0)
        logout_button.set_sensitive(False)
        self.logout_button = logout_button

        self.logged_in = False

    def on_login_clicked(self, button):
        # Get the username and password
        username = self.username_entry.get_text()
        password = self.password_entry.get_text()

        # Verify the username and password
        if username == "admin" and password == "password":
            print("Login Successful at ", self.get_timestamp())
            self.logged_in = True
            self.logout_button.set_sensitive(True)

            # Log the login time using gspread
            gc = gspread.service_account("service_account.json")
            sh = gc.open("Login Logs").sheet1
            sh.append_row([username, self.get_timestamp(), "Login"])
        else:
            print("Login Failed at ", self.get_timestamp())

    def on_logout_clicked(self, button):
        if self.logged_in:
            print("Logout Successful at ", self.get_timestamp())
            self.logged_in = False
            self.logout_button.set_sensitive(False)

            # Log the logout time using gspread
            # Log the logout time using gspread
            gc = gspread.service_account("service_account.json")
            sh = gc.open("Login Logs").sheet1
            sh.append_row([self.username_entry.get_text(), self.get_timestamp(), "Logout"])
        else:
            print("Logout Failed at ", self.get_timestamp())

    def get_timestamp(self):
        return str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

win = LoginWindow()
win.connect("delete-event", Gtk.main_quit)
win.show_all()
Gtk.main()

