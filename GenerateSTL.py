import numpy as np
import pandas as pd

def read_porous(filename, dimension = 2):
    if (dimension == 2):
        lattice = pd.read_csv(filename, delimiter=",", header=None, names=["X", "Y", "lattice"], index_col=None,skiprows=[0,1,2])
        return lattice

    # Real 3D program not supported yet
    if (dimension == 3):
        lattice = pd.read_csv(filename, delimiter=",", header=None, names=["X", "Y", "Z", "lattice"], index_col=None,skiprows=[0,1,2])
        return None

def create_solid_STL_quasi_3D(fo,lattice): 
    # vertices number = (M+1)*(N+1)
    # create quasi-3D porous media with 2D plt file
    M = int(lattice.loc[:,"X"].max() + 1)
    N = int(lattice.loc[:,"Y"].max() + 1)
    for i in range(M):
        for j in range(N):
            if (lattice.loc[i+j*M]["lattice"] == 1):
                faces = [0, 0, 0, 0]
                if (j == 0):
                    faces[0] = 1
                    if (i == 0):
                        faces[2] = 1
                        if (lattice.loc[i+1+j*M]["lattice"] == 0):
                            faces[3] = 1
                        if (lattice.loc[i+(j+1)*M]["lattice"] == 0):
                            faces[1] = 1
                    elif (i == M-1):
                        faces[3] = 1
                        if (lattice.loc[i+(j+1)*M]["lattice"] == 0):
                            faces[1] = 1
                        if (lattice.loc[i-1+j*M]["lattice"] == 0):
                            faces[2] = 1
                    else:
                        if (lattice.loc[i-1+j*M]["lattice"] == 0):
                            faces[2] = 1
                        if (lattice.loc[i+1+j*M]["lattice"] == 0):
                            faces[3] = 1
                        if (lattice.loc[i+(j+1)*M]["lattice"] == 0):
                            faces[1] = 1
                elif (j == N-1):
                    faces[1] = 1
                    if (i == 0):
                        faces[2] = 1
                        if (lattice.loc[i+(j-1)*M]["lattice"] == 0):
                            faces[0] = 1
                        if (lattice.loc[i+1+j*M]["lattice"] == 0):
                            faces[3] = 1
                    elif (i == M-1):
                        faces[3] = 1
                        if (lattice.loc[i+(j-1)*M]["lattice"] == 0):
                            faces[0] = 1
                        if (lattice.loc[i-1+j*M]["lattice"] == 0):
                            faces[2] = 1
                    else:
                        if (lattice.loc[i-1+j*M]["lattice"] == 0):
                            faces[2] = 1
                        if (lattice.loc[i+1+j*M]["lattice"] == 0):
                            faces[3] = 1
                        if (lattice.loc[i+(j-1)*M]["lattice"] == 0):
                            faces[0] = 1
                else:
                    if (i == 0):
                        faces[2] = 1
                        if (lattice.loc[i+1+j*M]["lattice"] == 0):
                            faces[3] = 1
                        if (lattice.loc[i+(j-1)*M]["lattice"] == 0):
                            faces[0] = 1
                        if (lattice.loc[i+(j+1)*M]["lattice"] == 0):
                            faces[1] = 1
                    elif (i == M-1):
                        faces[3] = 1
                        if (lattice.loc[i+(j-1)*M]["lattice"] == 0):
                            faces[0] = 1
                        if (lattice.loc[i+(j+1)*M]["lattice"] == 0):
                            faces[1] = 1
                        if (lattice.loc[i-1+j*M]["lattice"] == 0):
                            faces[2] = 1
                    else:
                        if (lattice.loc[i-1+j*M]["lattice"] == 0):
                            faces[2] = 1
                        if (lattice.loc[i+1+j*M]["lattice"] == 0):
                            faces[3] = 1
                        if (lattice.loc[i+(j-1)*M]["lattice"] == 0):
                            faces[0] = 1
                        if (lattice.loc[i+(j+1)*M]["lattice"] == 0):
                            faces[1] = 1

                success = create_faces(lattice, fo, faces, i, j)
                if success == False:
                    print("WriteFile ERROR at (%i,%i)" %(i,j))
                    return False

    return True

def create_faces (lattice, fo, faces, i, j, k=0, dx=1, dy=1, dz=1, scale=1e-3):
    # faces = [down, up, left, right]
    # one-hot code

    M = int(lattice.loc[:,"X"].max() + 1)    # i direction
    N = int(lattice.loc[:,"Y"].max() + 1)    # j direction
    
    if (i > M-1) or (j > N-1):
        return False

    if faces[0] == 1:   # down
        normal_vector = [0.0, -1.0, 0.0]
        v1 = np.array([i*dx, j*dy, k*dz]) * scale
        v2 = np.array([(i+1)*dx, j*dy, k*dz]) * scale
        v3 = np.array([(i+1)*dx, j*dy, (k+1)*dz]) * scale
        v4 = np.array([i*dx, j*dy, (k+1)*dz]) * scale

        fo.writelines("\t facet normal %f %f %f\n" %(normal_vector[0], normal_vector[1], normal_vector[2]))
        fo.writelines("\t \t outer loop\n")
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v1[0], v1[1], v1[2]))
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v2[0], v2[1], v2[2]))
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v3[0], v3[1], v3[2]))
        fo.writelines("\t \t endloop\n")
        fo.writelines("\t endfacet\n")

        fo.writelines("\t facet normal %f %f %f\n" %(normal_vector[0], normal_vector[1], normal_vector[2]))
        fo.writelines("\t \t outer loop\n")
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v1[0], v1[1], v1[2]))
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v3[0], v3[1], v3[2]))
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v4[0], v4[1], v4[2]))
        fo.writelines("\t \t endloop\n")
        fo.writelines("\t endfacet\n")

    if faces[1] == 1: # up
        normal_vector = [0.0, 1.0, 0.0]
        v1 = np.array([i*dx, (j+1)*dy, k*dz]) * scale
        v2 = np.array([(i+1)*dx, (j+1)*dy, k*dz]) * scale
        v3 = np.array([(i+1)*dx, (j+1)*dy, (k+1)*dz]) * scale
        v4 = np.array([i*dx, (j+1)*dy, (k+1)*dz]) * scale

        fo.writelines("\t facet normal %f %f %f\n" %(normal_vector[0], normal_vector[1], normal_vector[2]))
        fo.writelines("\t \t outer loop\n")
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v1[0], v1[1], v1[2]))
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v3[0], v3[1], v3[2]))
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v2[0], v2[1], v2[2]))
        fo.writelines("\t \t endloop\n")
        fo.writelines("\t endfacet\n")

        fo.writelines("\t facet normal %f %f %f\n" %(normal_vector[0], normal_vector[1], normal_vector[2]))
        fo.writelines("\t \t outer loop\n")
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v1[0], v1[1], v1[2]))
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v4[0], v4[1], v4[2]))
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v3[0], v3[1], v3[2]))
        fo.writelines("\t \t endloop\n")
        fo.writelines("\t endfacet\n")

    if faces[2] == 1: # left
        normal_vector = [-1.0, 0.0, 0.0]
        v1 = np.array([i*dx, j*dy, k*dz]) * scale
        v2 = np.array([i*dx, (j+1)*dy, k*dz]) * scale
        v3 = np.array([i*dx, (j+1)*dy, (k+1)*dz]) * scale
        v4 = np.array([i*dx, j*dy, (k+1)*dz]) * scale

        fo.writelines("\t facet normal %f %f %f\n" %(normal_vector[0], normal_vector[1], normal_vector[2]))
        fo.writelines("\t \t outer loop\n")
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v1[0], v1[1], v1[2]))
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v3[0], v3[1], v3[2]))
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v2[0], v2[1], v2[2]))
        fo.writelines("\t \t endloop\n")
        fo.writelines("\t endfacet\n")

        fo.writelines("\t facet normal %f %f %f\n" %(normal_vector[0], normal_vector[1], normal_vector[2]))
        fo.writelines("\t \t outer loop\n")
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v1[0], v1[1], v1[2]))
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v4[0], v4[1], v4[2]))
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v3[0], v3[1], v3[2]))
        fo.writelines("\t \t endloop\n")
        fo.writelines("\t endfacet\n")

    if faces[3] == 1: # right
        normal_vector = [1.0, 0.0, 0.0]
        v1 = np.array([(i+1)*dx, j*dy, k*dz]) * scale
        v2 = np.array([(i+1)*dx, (j+1)*dy, k*dz]) * scale
        v3 = np.array([(i+1)*dx, (j+1)*dy, (k+1)*dz]) * scale
        v4 = np.array([(i+1)*dx, j*dy, (k+1)*dz]) * scale

        fo.writelines("\t facet normal %f %f %f\n" %(normal_vector[0], normal_vector[1], normal_vector[2]))
        fo.writelines("\t \t outer loop\n")
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v1[0], v1[1], v1[2]))
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v2[0], v2[1], v2[2]))
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v3[0], v3[1], v3[2]))
        fo.writelines("\t \t endloop\n")
        fo.writelines("\t endfacet\n")

        fo.writelines("\t facet normal %f %f %f\n" %(normal_vector[0], normal_vector[1], normal_vector[2]))
        fo.writelines("\t \t outer loop\n")
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v1[0], v1[1], v1[2]))
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v3[0], v3[1], v3[2]))
        fo.writelines("\t \t \t vertex %f %f %f\n" %(v4[0], v4[1], v4[2]))
        fo.writelines("\t \t endloop\n")
        fo.writelines("\t endfacet\n")

    # front
    normal_vector_front = [0.0, 0.0, 1.0]
    v1 = np.array([i*dx, j*dy, (k+1)*dz]) * scale
    v2 = np.array([(i+1)*dx, j*dy, (k+1)*dz]) * scale
    v3 = np.array([(i+1)*dx, (j+1)*dy, (k+1)*dz]) * scale
    v4 = np.array([i*dx, (j+1)*dy, (k+1)*dz]) * scale

    fo.writelines("\t facet normal %f %f %f\n" %(normal_vector_front[0],normal_vector_front[1], normal_vector_front[2]))
    fo.writelines("\t \t outer loop\n")
    fo.writelines("\t \t \t vertex %f %f %f\n" %(v1[0], v1[1], v1[2]))
    fo.writelines("\t \t \t vertex %f %f %f\n" %(v2[0], v2[1], v2[2]))
    fo.writelines("\t \t \t vertex %f %f %f\n" %(v3[0], v3[1], v3[2]))
    fo.writelines("\t \t endloop\n")
    fo.writelines("\t endfacet\n")

    fo.writelines("\t facet normal %f %f %f\n" %(normal_vector_front[0], normal_vector_front[1],normal_vector_front[2]))
    fo.writelines("\t \t outer loop\n")
    fo.writelines("\t \t \t vertex %f %f %f\n" %(v1[0], v1[1], v1[2]))
    fo.writelines("\t \t \t vertex %f %f %f\n" %(v3[0], v3[1], v3[2]))
    fo.writelines("\t \t \t vertex %f %f %f\n" %(v4[0], v4[1], v4[2]))
    fo.writelines("\t \t endloop\n")
    fo.writelines("\t endfacet\n")

    # back
    normal_vector_back = [0.0, 0.0, -1.0]
    v1 = np.array([i*dx, j*dy, k*dz]) * scale
    v2 = np.array([(i+1)*dx, j*dy, k*dz]) * scale
    v3 = np.array([(i+1)*dx, (j+1)*dy, k*dz]) * scale
    v4 = np.array([i*dx, (j+1)*dy, k*dz]) * scale

    fo.writelines("\t facet normal %f %f %f\n" %(normal_vector_back[0],normal_vector_back[1], normal_vector_back[2]))
    fo.writelines("\t \t outer loop\n")
    fo.writelines("\t \t \t vertex %f %f %f\n" %(v1[0], v1[1], v1[2]))
    fo.writelines("\t \t \t vertex %f %f %f\n" %(v3[0], v3[1], v3[2]))
    fo.writelines("\t \t \t vertex %f %f %f\n" %(v2[0], v2[1], v2[2]))
    fo.writelines("\t \t endloop\n")
    fo.writelines("\t endfacet\n")

    fo.writelines("\t facet normal %f %f %f\n" %(normal_vector_back[0], normal_vector_back[1],normal_vector_back[2]))
    fo.writelines("\t \t outer loop\n")
    fo.writelines("\t \t \t vertex %f %f %f\n" %(v1[0], v1[1], v1[2]))
    fo.writelines("\t \t \t vertex %f %f %f\n" %(v4[0], v4[1], v4[2]))
    fo.writelines("\t \t \t vertex %f %f %f\n" %(v3[0], v3[1], v3[2]))
    fo.writelines("\t \t endloop\n")
    fo.writelines("\t endfacet\n")

    return True

def create_solid_STL_2D(fo, lattice, dx = 1, dy = 1, dz = 1, scale = 1e-3):
    M = int(lattice.loc[:,"X"].max() + 1)
    N = int(lattice.loc[:,"Y"].max() + 1)
    k = 0
    for i in range(M):
        for j in range(N):
            if (lattice.loc[i+j*M]["lattice"] == 1):
                normal_vector_back = [0.0, 0.0, -1.0]
                v1 = np.array([i*dx, j*dy, k*dz]) * scale
                v2 = np.array([(i+1)*dx, j*dy, k*dz]) * scale
                v3 = np.array([(i+1)*dx, (j+1)*dy, k*dz]) * scale
                v4 = np.array([i*dx, (j+1)*dy, k*dz]) * scale

                fo.writelines("\t facet normal %f %f %f\n" %(normal_vector_back[0],normal_vector_back[1], normal_vector_back[2]))
                fo.writelines("\t \t outer loop\n")
                fo.writelines("\t \t \t vertex %f %f %f\n" %(v1[0], v1[1], v1[2]))
                fo.writelines("\t \t \t vertex %f %f %f\n" %(v3[0], v3[1], v3[2]))
                fo.writelines("\t \t \t vertex %f %f %f\n" %(v2[0], v2[1], v2[2]))
                fo.writelines("\t \t endloop\n")
                fo.writelines("\t endfacet\n")

                fo.writelines("\t facet normal %f %f %f\n" %(normal_vector_back[0], normal_vector_back[1],normal_vector_back[2]))
                fo.writelines("\t \t outer loop\n")
                fo.writelines("\t \t \t vertex %f %f %f\n" %(v1[0], v1[1], v1[2]))
                fo.writelines("\t \t \t vertex %f %f %f\n" %(v4[0], v4[1], v4[2]))
                fo.writelines("\t \t \t vertex %f %f %f\n" %(v3[0], v3[1], v3[2]))
                fo.writelines("\t \t endloop\n")
                fo.writelines("\t endfacet\n")

    return True

def create_fluid_STL_2D(fo, lattice, dx = 1, dy = 1, dz = 1, scale = 1e-3):
    M = int(lattice.loc[:,"X"].max() + 1)
    N = int(lattice.loc[:,"Y"].max() + 1)
    k = 0
    for i in range(M):
        for j in range(N):
            if (lattice.loc[i+j*M]["lattice"] == 0):
                normal_vector_back = [0.0, 0.0, -1.0]
                v1 = np.array([i*dx, j*dy, k*dz]) * scale
                v2 = np.array([(i+1)*dx, j*dy, k*dz]) * scale
                v3 = np.array([(i+1)*dx, (j+1)*dy, k*dz]) * scale
                v4 = np.array([i*dx, (j+1)*dy, k*dz]) * scale

                fo.writelines("\t facet normal %f %f %f\n" %(normal_vector_back[0],normal_vector_back[1], normal_vector_back[2]))
                fo.writelines("\t \t outer loop\n")
                fo.writelines("\t \t \t vertex %f %f %f\n" %(v1[0], v1[1], v1[2]))
                fo.writelines("\t \t \t vertex %f %f %f\n" %(v3[0], v3[1], v3[2]))
                fo.writelines("\t \t \t vertex %f %f %f\n" %(v2[0], v2[1], v2[2]))
                fo.writelines("\t \t endloop\n")
                fo.writelines("\t endfacet\n")

                fo.writelines("\t facet normal %f %f %f\n" %(normal_vector_back[0], normal_vector_back[1],normal_vector_back[2]))
                fo.writelines("\t \t outer loop\n")
                fo.writelines("\t \t \t vertex %f %f %f\n" %(v1[0], v1[1], v1[2]))
                fo.writelines("\t \t \t vertex %f %f %f\n" %(v4[0], v4[1], v4[2]))
                fo.writelines("\t \t \t vertex %f %f %f\n" %(v3[0], v3[1], v3[2]))
                fo.writelines("\t \t endloop\n")
                fo.writelines("\t endfacet\n")

    return True

def Example_2D():
    fo = open("porous.stl","w")

    # Read 2D plt file
    filename = "PorMed_2D.plt"

    Dimension = 3

    block_string = "block" 
    fo.writelines("solid %s\n" %(block_string) )

    lattice = read_porous(filename)

    if (Dimension == 2):
        Done = create_fluid_STL_2D(fo, lattice)
    elif (Dimension == 3):
        Done = create_solid_STL_quasi_3D(fo, lattice)
    else:
        Done = False

    if Done == False:
        print("Creat STL file failure!")
    else:
        fo.writelines("endsolid\n")

    fo.close()

def Example_3D():
    fo = open("porous.stl","w")

    # Read 3D plt file
    filename = "PorMed_3D.plt"

    # Real 3D program not supported yet

    fo.close()

if __name__ == "__main__":
    # 2D example
    Example_2D()

    #########################
    # 3D example
    # Real 3D program not supported yet
    Example_3D()
    