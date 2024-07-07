# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This is a way to encode and decode data from an image.


"""

def ordinal_to_rgba(ordinal):
    ordinal = ordinal & 0xFFFFFFFF
    
    r = int(ordinal & 0xFF)
    g = int((ordinal >> 8) & 0xFF)
    b = int((ordinal >> 16) & 0xFF)
    a = int((ordinal >> 24) & 0xFF)

    return r, g, b, a

def rgba_to_ordinal(r, g, b, a):
    r = int(r) & 0xFF
    g = int(g) & 0xFF
    b = int(b) & 0xFF
    a = int(a) & 0xFF
    
    ordinal = (a << 24) | (b << 16) | (g << 8) | r
    
    return ordinal
