# Dylan Kenneth Eliot
"""
Fixed bouncing-ball demo (ctypes + X11).

Fixes:
 - creates its own GC
 - allocates "white" color from the default colormap and sets GC foreground
 - draws using that GC so the ball is visible

the more you define, the more control ya have.

Good luck.
"""
import ctypes
import ctypes.util
import time
import sys
from ctypes import (
    Structure, Union, POINTER, byref,
    c_int, c_uint, c_ulong, c_long, c_char_p, c_ushort, c_char, c_void_p
)

# --- load libX11 ---
libname = ctypes.util.find_library("X11")
if not libname:
    print("libX11 not found. Install X11 (e.g. libx11-dev).")
    sys.exit(1)
X11 = ctypes.CDLL(libname)

# --- basic types & XEvent with padding (64-bit safe) ---
class Display(Structure):
    pass
DisplayPtr = POINTER(Display)

class XKeyEvent(Structure):
    _fields_ = [
        ("type", c_int),
        ("serial", c_ulong),
        ("send_event", c_int),
        ("display", DisplayPtr),
        ("window", c_ulong),
        ("root", c_ulong),
        ("subwindow", c_ulong),
        ("time", c_ulong),
        ("x", c_int),
        ("y", c_int),
        ("x_root", c_int),
        ("y_root", c_int),
        ("state", c_uint),
        ("keycode", c_uint),
        ("same_screen", c_int),
    ]

class XEvent(Union):
    _fields_ = [
        ("type", c_int),
        ("xkey", XKeyEvent),
        ("pad", c_long * 24),
    ]

# --- XColor (for XAllocNamedColor) ---
class XColor(Structure):
    _fields_ = [
        ("pixel", c_ulong),
        ("red", c_ushort),
        ("green", c_ushort),
        ("blue", c_ushort),
        ("flags", c_char),
        ("pad", c_char),
    ]

# --- prototypes (important!) ---
X11.XOpenDisplay.restype = DisplayPtr
X11.XOpenDisplay.argtypes = [c_char_p]

X11.XDefaultScreen.restype = c_int
X11.XDefaultScreen.argtypes = [DisplayPtr]

X11.XRootWindow.restype = c_ulong
X11.XRootWindow.argtypes = [DisplayPtr, c_int]

X11.XCreateSimpleWindow.restype = c_ulong
X11.XCreateSimpleWindow.argtypes = [
    DisplayPtr, c_ulong, c_int, c_int, c_uint, c_uint, c_uint, c_ulong, c_ulong
]

X11.XSelectInput.restype = c_int
X11.XSelectInput.argtypes = [DisplayPtr, c_ulong, c_long]

X11.XMapWindow.restype = c_int
X11.XMapWindow.argtypes = [DisplayPtr, c_ulong]

X11.XDefaultGC.restype = c_ulong
X11.XDefaultGC.argtypes = [DisplayPtr, c_int]

X11.XPending.restype = c_int
X11.XPending.argtypes = [DisplayPtr]

X11.XNextEvent.restype = None
X11.XNextEvent.argtypes = [DisplayPtr, POINTER(XEvent)]

X11.XClearWindow.restype = None
X11.XClearWindow.argtypes = [DisplayPtr, c_ulong]

X11.XFillArc.restype = None
X11.XFillArc.argtypes = [DisplayPtr, c_ulong, c_ulong, c_int, c_int, c_uint, c_uint, c_int, c_int]

X11.XFlush.restype = None
X11.XFlush.argtypes = [DisplayPtr]

X11.XCloseDisplay.restype = None
X11.XCloseDisplay.argtypes = [DisplayPtr]

# Color / GC prototypes
X11.XDefaultColormap.restype = c_ulong
X11.XDefaultColormap.argtypes = [DisplayPtr, c_int]

X11.XAllocNamedColor.restype = c_int
X11.XAllocNamedColor.argtypes = [DisplayPtr, c_ulong, c_char_p, POINTER(XColor), POINTER(XColor)]

X11.XCreateGC.restype = c_ulong
X11.XCreateGC.argtypes = [DisplayPtr, c_ulong, c_ulong, c_void_p]

X11.XSetForeground.restype = None
X11.XSetForeground.argtypes = [DisplayPtr, c_ulong, c_ulong]

X11.XFreeGC.restype = None
X11.XFreeGC.argtypes = [DisplayPtr, c_ulong]

# --- X event/type constants used here ---
Expose = 12
KeyPress = 2
ButtonPress = 4
PointerMotion = 6

ExposureMask = (1 << 15)
KeyPressMask = (1 << 0)
ButtonPressMask = (1 << 2)
PointerMotionMask = (1 << 6)

# --- demo parameters ---
WIN_X, WIN_Y = 100, 100
WIN_W, WIN_H = 640, 480
BALL_RADIUS = 20
FPS = 60.0
FRAME_TIME = 1.0 / FPS

def main():
    dpy = X11.XOpenDisplay(None)
    if not dpy:
        print("Unable to open X display (is $DISPLAY set?)")
        return

    scr = X11.XDefaultScreen(dpy)
    root = X11.XRootWindow(dpy, scr)

    # Create a window (parent is root_window)
    win = X11.XCreateSimpleWindow(
        dpy, root,
        WIN_X, WIN_Y,      # x, y
        WIN_W, WIN_H,      # width, height
        1,                 # border width
        0, 0               # border and background (use defaults)
    )

    # Select events we care about
    masks = ExposureMask | KeyPressMask | ButtonPressMask | PointerMotionMask
    X11.XSelectInput(dpy, win, masks)
    X11.XMapWindow(dpy, win)

    # Create a dedicated GC for drawing
    # Get default colormap and allocate the "white" color
    colormap = X11.XDefaultColormap(dpy, scr)
    screen_color = XColor()
    exact_color = XColor()
    ok = X11.XAllocNamedColor(dpy, colormap, b"white", byref(screen_color), byref(exact_color))
    if not ok:
        print("XAllocNamedColor failed; drawing may be invisible.")
        white_pixel = 0  # fallback; might be black on many setups
    else:
        white_pixel = screen_color.pixel
    # create a GC attached to our window drawable
    gc = X11.XCreateGC(dpy, win, 0, None)
    # set the foreground of this GC to the allocated white pixel
    X11.XSetForeground(dpy, gc, white_pixel)
    print(f"Allocated white pixel = {white_pixel}, gc = {gc}")

    # ball state
    x = WIN_W // 2
    y = WIN_H // 2
    dx = 200  # pixels per second
    dy = 150  # pixels per second

    ev = XEvent()
    running = True
    last_time = time.perf_counter()

    print("Bouncing ball running. Click -> move ball. Press any key -> exit.")

    while running:
        # handle events (process all pending first)
        while X11.XPending(dpy):
            X11.XNextEvent(dpy, byref(ev))
            ev_type = ev.type
            if ev_type == Expose:
                # explicit redraw on expose if desired
                pass
            elif ev_type == KeyPress:
                running = False
                break
            elif ev_type == ButtonPress:
                bx = ev.xkey.x
                by = ev.xkey.y
                x = max(BALL_RADIUS, min(WIN_W - BALL_RADIUS, bx))
                y = max(BALL_RADIUS, min(WIN_H - BALL_RADIUS, by))

        # timing
        now = time.perf_counter()
        dt = now - last_time
        if dt <= 0:
            dt = FRAME_TIME
        last_time = now

        # update physics
        x += int(dx * dt)
        y += int(dy * dt)

        # bounce on edges
        if x - BALL_RADIUS <= 0:
            x = BALL_RADIUS
            dx = -dx
        if x + BALL_RADIUS >= WIN_W:
            x = WIN_W - BALL_RADIUS
            dx = -dx
        if y - BALL_RADIUS <= 0:
            y = BALL_RADIUS
            dy = -dy
        if y + BALL_RADIUS >= WIN_H:
            y = WIN_H - BALL_RADIUS
            dy = -dy

        # clear & draw using our GC (which has the white foreground)
        X11.XClearWindow(dpy, win)
        X11.XFillArc(dpy, win, gc,
                     x - BALL_RADIUS, y - BALL_RADIUS,
                     BALL_RADIUS * 2, BALL_RADIUS * 2,
                     0, 360 * 64)
        X11.XFlush(dpy)

        # frame pacing
        elapsed = time.perf_counter() - now
        to_sleep = FRAME_TIME - elapsed
        if to_sleep > 0:
            time.sleep(to_sleep)

    # cleanup
    X11.XFreeGC(dpy, gc)
    X11.XCloseDisplay(dpy)
    print("Exited cleanly.")

if __name__ == "__main__":
    main()
