import ctypes
from ctypes import wintypes

# ======================================
# Correct Win64 Types
# ======================================

PTR_SIZE = ctypes.sizeof(ctypes.c_void_p)

if PTR_SIZE == 8:
    WPARAM  = ctypes.c_uint64
    LPARAM  = ctypes.c_int64
    LRESULT = ctypes.c_int64
else:
    WPARAM  = ctypes.c_uint32
    LPARAM  = ctypes.c_int32
    LRESULT = ctypes.c_int32

HANDLE = ctypes.c_void_p
HWND = HANDLE
HDC = HANDLE
HGLRC = HANDLE
HINSTANCE = HANDLE
LPCWSTR = ctypes.c_wchar_p
UINT = ctypes.c_uint
DWORD = ctypes.c_uint32
WORD = ctypes.c_uint16

# ======================================
# Load DLLs
# ======================================

user32 = ctypes.WinDLL("user32", use_last_error=True)
gdi32 = ctypes.WinDLL("gdi32", use_last_error=True)
opengl32 = ctypes.WinDLL("opengl32", use_last_error=True)
kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)

# Define DefWindowProcW correctly
user32.DefWindowProcW.argtypes = [HWND, UINT, WPARAM, LPARAM]
user32.DefWindowProcW.restype = LRESULT

# ======================================
# Constants
# ======================================

CS_OWNDC = 0x0020
WM_DESTROY = 0x0002
WM_CLOSE = 0x0010
WM_QUIT = 0x0012
WS_OVERLAPPEDWINDOW = 0x00CF0000
CW_USEDEFAULT = 0x80000000
SW_SHOW = 5

PFD_TYPE_RGBA = 0
PFD_MAIN_PLANE = 0
PFD_DOUBLEBUFFER = 0x00000001
PFD_DRAW_TO_WINDOW = 0x00000004
PFD_SUPPORT_OPENGL = 0x00000020

GL_COLOR_BUFFER_BIT = 0x00004000
GL_DEPTH_BUFFER_BIT = 0x00000100
GL_QUADS = 0x0007
GL_DEPTH_TEST = 0x0B71
GL_PROJECTION = 0x1701
GL_MODELVIEW = 0x1700

# ======================================
# OpenGL Function Signatures (CRITICAL)
# ======================================

opengl32.glClearColor.argtypes = [ctypes.c_float]*4
opengl32.glClearColor.restype = None

opengl32.glClear.argtypes = [ctypes.c_uint]
opengl32.glClear.restype = None

opengl32.glEnable.argtypes = [ctypes.c_uint]
opengl32.glEnable.restype = None

opengl32.glMatrixMode.argtypes = [ctypes.c_uint]
opengl32.glMatrixMode.restype = None

opengl32.glLoadIdentity.argtypes = []
opengl32.glLoadIdentity.restype = None

opengl32.glFrustum.argtypes = [ctypes.c_double]*6
opengl32.glFrustum.restype = None

opengl32.glTranslatef.argtypes = [ctypes.c_float]*3
opengl32.glTranslatef.restype = None

opengl32.glRotatef.argtypes = [ctypes.c_float]*4
opengl32.glRotatef.restype = None

opengl32.glBegin.argtypes = [ctypes.c_uint]
opengl32.glBegin.restype = None

opengl32.glVertex3f.argtypes = [ctypes.c_float]*3
opengl32.glVertex3f.restype = None

opengl32.glEnd.argtypes = []
opengl32.glEnd.restype = None

# ======================================
# Window Procedure
# ======================================

@ctypes.WINFUNCTYPE(LRESULT, HWND, UINT, WPARAM, LPARAM)
def WndProc(hwnd, msg, wparam, lparam):
    if msg == WM_CLOSE:
        user32.DestroyWindow(hwnd)
        return 0
    elif msg == WM_DESTROY:
        user32.PostQuitMessage(0)
        return 0
    return user32.DefWindowProcW(hwnd, msg, wparam, lparam)

# ======================================
# Create Window
# ======================================

def create_window():
    hInstance = kernel32.GetModuleHandleW(None)
    className = "OpenGLWindow"

    class WNDCLASS(ctypes.Structure):
        _fields_ = [
            ("style", UINT),
            ("lpfnWndProc", ctypes.c_void_p),
            ("cbClsExtra", ctypes.c_int),
            ("cbWndExtra", ctypes.c_int),
            ("hInstance", HINSTANCE),
            ("hIcon", HANDLE),
            ("hCursor", HANDLE),
            ("hbrBackground", HANDLE),
            ("lpszMenuName", LPCWSTR),
            ("lpszClassName", LPCWSTR),
        ]

    wndclass = WNDCLASS()
    wndclass.style = CS_OWNDC
    wndclass.lpfnWndProc = ctypes.cast(WndProc, ctypes.c_void_p)
    wndclass.hInstance = hInstance
    wndclass.lpszClassName = className

    user32.RegisterClassW(ctypes.byref(wndclass))

    hwnd = user32.CreateWindowExW(
        0, className, "Wine64 Cube",
        WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, CW_USEDEFAULT,
        800, 600,
        None, None,
        hInstance, None
    )

    user32.ShowWindow(hwnd, SW_SHOW)
    return hwnd

# ======================================
# Setup OpenGL
# ======================================

def setup_opengl(hwnd):
    hdc = user32.GetDC(hwnd)

    class PIXELFORMATDESCRIPTOR(ctypes.Structure):
        _fields_ = [
            ("nSize", WORD),
            ("nVersion", WORD),
            ("dwFlags", DWORD),
            ("iPixelType", ctypes.c_byte),
            ("cColorBits", ctypes.c_byte),
            ("cRedBits", ctypes.c_byte),
            ("cRedShift", ctypes.c_byte),
            ("cGreenBits", ctypes.c_byte),
            ("cGreenShift", ctypes.c_byte),
            ("cBlueBits", ctypes.c_byte),
            ("cBlueShift", ctypes.c_byte),
            ("cAlphaBits", ctypes.c_byte),
            ("cAlphaShift", ctypes.c_byte),
            ("cAccumBits", ctypes.c_byte),
            ("cAccumRedBits", ctypes.c_byte),
            ("cAccumGreenBits", ctypes.c_byte),
            ("cAccumBlueBits", ctypes.c_byte),
            ("cAccumAlphaBits", ctypes.c_byte),
            ("cDepthBits", ctypes.c_byte),
            ("cStencilBits", ctypes.c_byte),
            ("cAuxBuffers", ctypes.c_byte),
            ("iLayerType", ctypes.c_byte),
            ("bReserved", ctypes.c_byte),
            ("dwLayerMask", DWORD),
            ("dwVisibleMask", DWORD),
            ("dwDamageMask", DWORD),
        ]

    pfd = PIXELFORMATDESCRIPTOR()
    pfd.nSize = ctypes.sizeof(PIXELFORMATDESCRIPTOR)
    pfd.nVersion = 1
    pfd.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER
    pfd.iPixelType = PFD_TYPE_RGBA
    pfd.cColorBits = 24
    pfd.cDepthBits = 24
    pfd.iLayerType = PFD_MAIN_PLANE

    pf = gdi32.ChoosePixelFormat(hdc, ctypes.byref(pfd))
    gdi32.SetPixelFormat(hdc, pf, ctypes.byref(pfd))

    hrc = opengl32.wglCreateContext(hdc)
    opengl32.wglMakeCurrent(hdc, hrc)

    opengl32.glEnable(GL_DEPTH_TEST)
    opengl32.glClearColor(0.1, 0.1, 0.2, 1.0)

    opengl32.glMatrixMode(GL_PROJECTION)
    opengl32.glLoadIdentity()
    opengl32.glFrustum(-1.33, 1.33, -1, 1, 1, 100)
    opengl32.glMatrixMode(GL_MODELVIEW)

    return hdc

# ======================================
# Render Loop
# ======================================

angle = 0.0

GL_SMOOTH = 0x1D01
opengl32.glShadeModel.argtypes = [ctypes.c_uint]
opengl32.glShadeModel.restype = None
opengl32.glShadeModel(GL_SMOOTH)

opengl32.glColor3f.argtypes = [ctypes.c_float]*3
opengl32.glColor3f.restype = None


def render():
    global angle

    opengl32.glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    opengl32.glLoadIdentity()
    opengl32.glTranslatef(0.0, 0.0, -6.0)
    opengl32.glRotatef(angle, 1.0, 1.0, 0.0)

    opengl32.glBegin(GL_QUADS)

    # FRONT
    opengl32.glColor3f(1,0,0) ; opengl32.glVertex3f(-1,-1, 1)
    opengl32.glColor3f(0,1,0) ; opengl32.glVertex3f( 1,-1, 1)
    opengl32.glColor3f(0,0,1) ; opengl32.glVertex3f( 1, 1, 1)
    opengl32.glColor3f(1,1,0) ; opengl32.glVertex3f(-1, 1, 1)

    # BACK
    opengl32.glColor3f(1,0,1) ; opengl32.glVertex3f(-1,-1,-1)
    opengl32.glColor3f(0,1,1) ; opengl32.glVertex3f(-1, 1,-1)
    opengl32.glColor3f(1,1,1) ; opengl32.glVertex3f( 1, 1,-1)
    opengl32.glColor3f(0,0,0) ; opengl32.glVertex3f( 1,-1,-1)

    # LEFT
    opengl32.glColor3f(1,0.5,0) ; opengl32.glVertex3f(-1,-1,-1)
    opengl32.glColor3f(0,0.5,1) ; opengl32.glVertex3f(-1,-1, 1)
    opengl32.glColor3f(0.5,1,0) ; opengl32.glVertex3f(-1, 1, 1)
    opengl32.glColor3f(1,0,0.5) ; opengl32.glVertex3f(-1, 1,-1)

    # RIGHT
    opengl32.glColor3f(0.2,0.8,0.4) ; opengl32.glVertex3f(1,-1,-1)
    opengl32.glColor3f(0.8,0.2,0.6) ; opengl32.glVertex3f(1, 1,-1)
    opengl32.glColor3f(0.3,0.7,1.0) ; opengl32.glVertex3f(1, 1, 1)
    opengl32.glColor3f(1.0,0.4,0.2) ; opengl32.glVertex3f(1,-1, 1)

    # TOP
    opengl32.glColor3f(0.9,0.9,0.2) ; opengl32.glVertex3f(-1,1,-1)
    opengl32.glColor3f(0.2,0.9,0.9) ; opengl32.glVertex3f(-1,1, 1)
    opengl32.glColor3f(0.9,0.2,0.9) ; opengl32.glVertex3f( 1,1, 1)
    opengl32.glColor3f(0.4,0.4,0.4) ; opengl32.glVertex3f( 1,1,-1)

    # BOTTOM
    opengl32.glColor3f(0.6,0.1,0.1) ; opengl32.glVertex3f(-1,-1,-1)
    opengl32.glColor3f(0.1,0.6,0.1) ; opengl32.glVertex3f( 1,-1,-1)
    opengl32.glColor3f(0.1,0.1,0.6) ; opengl32.glVertex3f( 1,-1, 1)
    opengl32.glColor3f(0.8,0.8,0.8) ; opengl32.glVertex3f(-1,-1, 1)

    opengl32.glEnd()

    angle += 0.7

def main():
    hwnd = create_window()
    hdc = setup_opengl(hwnd)

    msg = wintypes.MSG()

    while True:
        while user32.PeekMessageW(ctypes.byref(msg), None, 0, 0, 1):
            user32.TranslateMessage(ctypes.byref(msg))
            user32.DispatchMessageW(ctypes.byref(msg))
            if msg.message == WM_QUIT:
                return

        render()
        gdi32.SwapBuffers(hdc)

if __name__ == "__main__":
    main()
