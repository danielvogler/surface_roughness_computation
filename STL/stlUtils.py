import struct

def unpack (f, sig, l):
    s = f.read (l)
    return struct.unpack(sig, s)

def read_triangle(f,triangles,points,normals):
    n = unpack(f,"<3f", 12)
    p1 = unpack(f,"<3f", 12)
    p2 = unpack(f,"<3f", 12)
    p3 = unpack(f,"<3f", 12)
    b = unpack(f,"<h", 2)

    normals.append(n)
    l = len(points)
    points.append(list(p1))
    points.append(list(p2))
    points.append(list(p3))
    triangles.append([l, l+1, l+2])
    #bytecount.append(b[0])

def read_triangle_i(f,triangles,points,normals,i):
    n = unpack(f,"<3f", 12)
    p1 = unpack(f,"<3f", 12)
    p2 = unpack(f,"<3f", 12)
    p3 = unpack(f,"<3f", 12)
    b = unpack(f,"<h", 2)

    normals[i] = n
    l = len(points)
    points[i*3] = p1
    points[i*3+1] = p2
    points[i*3+2] = p3
    triangles[i] = [l, l+1, l+2]
    #bytecount.append(b[0])


def read_length(f):
    length = struct.unpack("@i", f.read(4))
    return length[0]

def read_header(f):
    f.seek(f.tell()+80)
    
def importSTLfile(filename):
    triangles = []
    points = []
    normals = []
    try:
        f = open ( filename, "rb")

        read_header(f)
        l = read_length(f)
        print "l:",l
        try:
            triangles = [None] * l
            points    = [None] * (3*l)
            normals   = [None] * l
            for i in range(l):
                read_triangle_i(f,triangles,points,normals,i)
        except Exception, e:
            print "Exception",e[0]
        # print len(normals), len(points), len(triangles)

    except Exception, e:
        print e
    return triangles,points,normals
    


