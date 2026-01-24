# Dylan Kenneth Eliot

"""
render opengl on posix, including rendering border.

This allows for rendering without x11 required. This means Mesa+vulkan/wayland can make use of raw tty's as canvases.
ModernGL only has to render the GL/GLSL frames that are then mapped to posix color maps.

Currently, this includes enough colorspace for full/true-color.

This also allows for a environment that is able to render displayable content.

Later on, mouse and keyboard controls can be added.

Currently useful with gotty for portable terminal game/terminal based renderer.
  Or useful for 3d environment to build the 3d portions of an AI layer that requires a opengl 'canvas' without vendor lockin.

"""


import sys
import time
import numpy as np
import moderngl
from pyrr import Matrix44

# Outer character window (inside border)
OUTER_W = 100
OUTER_H = 50

# Inner logical drawing area (smaller "character size" visually)
INNER_W = 30
INNER_H = 20

# Derived pixel resolution for inner braille area
PIXEL_W = INNER_W * 2
PIXEL_H = INNER_H * 4

LUMA_THRESHOLD = 1

SYS_CLEAR   = "\x1b[2J"
CURSOR_HOME = "\x1b[H"
HIDE_CURSOR = "\x1b[?25l"
SHOW_CURSOR = "\x1b[?25h"
RESET       = "\x1b[0m"

def rgb_fg_ansi(r, g, b):
    return f"\x1b[38;2;{r};{g};{b}m"

def braille_char_from_mask(mask):
    return chr(0x2800 + mask)

DOT_BIT_LUT = np.array([
    0, 3,
    1, 4,
    2, 5,
    6, 7,
], dtype=np.uint8)

vertices = np.array([
    -1,-1,-1, 255,0,0,
     1,-1,-1, 0,255,0,
     1, 1,-1, 0,0,255,
    -1, 1,-1, 255,255,0,
    -1,-1, 1, 0,255,255,
     1,-1, 1, 255,0,255,
     1, 1, 1, 255,255,255,
    -1, 1, 1, 0,0,0,
], dtype='f4')

indices = np.array([
    0,1,2, 2,3,0,
    4,5,6, 6,7,4,
    0,4,7, 7,3,0,
    1,5,6, 6,2,1,
    3,2,6, 6,7,3,
    0,1,5, 5,4,0,
], dtype='i4')

ctx = moderngl.create_standalone_context()

color_tex = ctx.texture((PIXEL_W, PIXEL_H), 3)
depth_rb = ctx.depth_renderbuffer((PIXEL_W, PIXEL_H))
fbo = ctx.framebuffer(color_attachments=[color_tex], depth_attachment=depth_rb)

prog = ctx.program(
    vertex_shader="""
        #version 330
        uniform mat4 mvp;
        in vec3 in_position;
        in vec3 in_color;
        out vec3 v_color;
        void main() {
            gl_Position = mvp * vec4(in_position, 1.0);
            v_color = in_color / 255.0;
        }
    """,
    fragment_shader="""
        #version 330
        in vec3 v_color;
        out vec4 f_color;
        void main() {
            f_color = vec4(v_color, 1.0);
        }
    """
)

vbo = ctx.buffer(vertices.tobytes())
ibo = ctx.buffer(indices.tobytes())
vao = ctx.vertex_array(prog, [(vbo, '3f 3f', 'in_position', 'in_color')], ibo)

# Aspect ratio based on inner pixel area
proj = Matrix44.perspective_projection(45.0, PIXEL_W / PIXEL_H, 0.1, 100.0)
lookat = Matrix44.look_at((4,3,3), (0,0,0), (0,1,0))
mvp_base = proj * lookat

angle = 0.0

# Precompute coords for inner braille grid
cell_pixel_coords = np.empty((INNER_H, INNER_W, 8, 2), dtype=np.int32)
for cy in range(INNER_H):
    for cx in range(INNER_W):
        px0 = cx * 2
        py0 = cy * 4
        idx = 0
        for y in range(4):
            py = py0 + y
            for x in range(2):
                px = px0 + x
                cell_pixel_coords[cy, cx, idx, 0] = px
                cell_pixel_coords[cy, cx, idx, 1] = py
                idx += 1

# Lines only for inner area
inner_lines = [""] * INNER_H

TOP_LEFT     = "┌"
TOP_RIGHT    = "┐"
BOTTOM_LEFT  = "└"
BOTTOM_RIGHT = "┘"
HORIZ        = "─"
VERT         = "│"

# Pad inner image to be centered in OUTER_W x OUTER_H
PAD_X = max(0, (OUTER_W - INNER_W) // 2)
PAD_Y_TOP = max(0, (OUTER_H - INNER_H) // 2)
PAD_Y_BOTTOM = OUTER_H - INNER_H - PAD_Y_TOP

sys.stdout.write(HIDE_CURSOR + SYS_CLEAR)
sys.stdout.flush()

try:
    while True:
        fbo.use()
        ctx.enable(moderngl.DEPTH_TEST)
        ctx.clear(0.0, 0.0, 0.0, 1.0)

        angle += 2.0
        rot = Matrix44.from_eulers((np.radians(angle), np.radians(angle), 0.0))
        mvp = mvp_base * rot
        prog['mvp'].write(mvp.astype('f4').tobytes())

        vao.render()

        data = fbo.read(components=3)
        img = np.frombuffer(data, dtype=np.uint8).reshape((PIXEL_H, PIXEL_W, 3))
        img = np.flipud(img)

        dot_bit_lut = DOT_BIT_LUT
        threshold = LUMA_THRESHOLD
        rgb_ansi = rgb_fg_ansi
        braille_chr = braille_char_from_mask
        reset = RESET
        img_local = img
        coords_all = cell_pixel_coords

        for cy in range(INNER_H):
            line_buf = []
            coords_row = coords_all[cy]
            for cx in range(INNER_W):
                coords = coords_row[cx]
                mask = 0
                active_r = 0
                active_g = 0
                active_b = 0
                active_count = 0

                for i in range(8):
                    px = coords[i, 0]
                    py = coords[i, 1]
                    r = int(img_local[py, px, 0])
                    g = int(img_local[py, px, 1])
                    b = int(img_local[py, px, 2])
                    luma = 0.2126 * r + 0.7152 * g + 0.0722 * b
                    if luma >= threshold:
                        lx = px & 1
                        ly = py & 3
                        bit = int(dot_bit_lut[lx + (ly << 1)])
                        mask |= (1 << bit)
                        active_r += r
                        active_g += g
                        active_b += b
                        active_count += 1

                if mask == 0 or active_count == 0:
                    line_buf.append(reset + " ")
                else:
                    inv = 1.0 / active_count
                    avg_r = int(active_r * inv)
                    avg_g = int(active_g * inv)
                    avg_b = int(active_b * inv)
                    seq = rgb_ansi(avg_r, avg_g, avg_b)
                    ch = braille_chr(mask)
                    line_buf.append(seq + ch + reset)

            inner_lines[cy] = "".join(line_buf)

        # Compose full 80x80 region with centered inner smaller image
        # Top/bottom padding rows are spaces
        empty_row = " " * OUTER_W
        padded_lines = []

        for _ in range(PAD_Y_TOP):
            padded_lines.append(empty_row)

        pad_left = " " * PAD_X
        pad_right = " " * (OUTER_W - INNER_W - PAD_X)
        for line in inner_lines:
            padded_lines.append(pad_left + line + pad_right)

        for _ in range(PAD_Y_BOTTOM):
            padded_lines.append(empty_row)

        # Add border around OUTER_W x OUTER_H region
        top_border    = TOP_LEFT + (HORIZ * OUTER_W) + TOP_RIGHT
        bottom_border = BOTTOM_LEFT + (HORIZ * OUTER_W) + BOTTOM_RIGHT
        bordered_lines = [top_border]
        for line in padded_lines:
            bordered_lines.append(VERT + line + VERT)
        bordered_lines.append(bottom_border)

        sys.stdout.write(CURSOR_HOME + "\n".join(bordered_lines))
        sys.stdout.flush()

        time.sleep(1/60)

except KeyboardInterrupt:
    pass
finally:
    sys.stdout.write(RESET + SHOW_CURSOR + "\n")
    sys.stdout.flush()
