import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk
from PIL import Image

class AreaMapWindow(Gtk.Window):

    def __init__(self):
        Gtk.Window.__init__(self, title="Image Map Creator")

        # Set up the GUI
        vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
        self.add(vbox)

        hbox = Gtk.Box(spacing=6)
        vbox.pack_start(hbox, True, True, 0)

        self.image = Gtk.Image()
        hbox.pack_start(self.image, True, True, 0)

        self.area_map = Gtk.DrawingArea()
        hbox.pack_start(self.area_map, True, True, 0)

        self.area_map.connect('draw', self.draw_areas)

        self.coord_list = Gtk.TextView()
        self.coord_list.set_editable(False)
        self.coord_list.set_cursor_visible(False)
        self.coord_list_buffer = self.coord_list.get_buffer()

        vbox.pack_start(self.coord_list, True, True, 0)

        # Load the image and set up the areas
        self.img = Image.open('image.jpg')
        self.width, self.height = self.img.size

        self.polygons = []
        self.current_polygon = []

        self.image.set_from_file('image.jpg')

        # Set up the event handlers
        self.area_map.connect('button-press-event', self.on_area_map_button_press)
        self.area_map.connect('motion-notify-event', self.on_area_map_motion_notify)
        self.connect('key-press-event', self.on_key_press)

    def draw_areas(self, area, cr):
        # Draw the polygons on the area map
        for poly in self.polygons:
            cr.set_source_rgba(1, 0, 0, 0.5)
            cr.move_to(*poly[0])
            for point in poly[1:]:
                cr.line_to(*point)
            cr.fill()

    def update_area_map(self):
        # Redraw the area map
        self.area_map.queue_draw()

        # Update the coordinate list
        self.coord_list_buffer.set_text('')
        for i, poly in enumerate(self.polygons):
            self.coord_list_buffer.insert_at_cursor('Polygon {}: {}\n'.format(i+1, poly))

    def on_area_map_button_press(self, widget, event):
        # Start a new polygon on left-click
        if event.button == 1:
            self.current_polygon = [self.convert_coordinates((event.x, event.y))]
            self.polygons.append(self.current_polygon)
            self.update_area_map()

    def on_area_map_motion_notify(self, widget, event):
        # Add points to the current polygon on mouse motion
        if event.state & Gdk.ModifierType.BUTTON1_MASK:
            self.current_polygon.append(self.convert_coordinates((event.x, event.y)))
            self.update_area_map()

    def on_key_press(self, widget, event):
        # Remove the last point from the current polygon on backspace
        if event.keyval == Gdk.KEY_BackSpace and self.current_polygon:
            self.current_polygon.pop()
            self.update_area_map()

    def convert_coordinates(self, coord):
        # Convert screen coordinates to image coordinates
        return (int(coord[0] / self.area_map.get_allocated_width() * self.width),
                int(coord[1] / self.area_map.get_allocated_height() * self.height))

    def create_image_map(self, filename):
        # Create the HTML for the image map
        html = '<img src="image.jpg" usemap="#image-map" id="image">\n'
        html += '<map name="image-map">\n'
        for i, poly in enumerate(self.polygons):
            coords = ','.join(['{},{}'.format(*p) for p in poly])
            html += '\t<area shape="poly" coords="{}" href="#" alt="Polygon {}">\n'.format(coords, i+1)
        html += '</map>'

        # Save the HTML to a file
        with open(filename, 'w') as f:
            f.write(html)

        print('Image map created: {}'.format(filename))

win = AreaMapWindow()
win.connect('destroy', Gtk.main_quit)
win.show_all()
Gtk.main()
