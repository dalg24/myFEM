# Python code
# Author: Bruno Turcksin
# Date: 2010-04-10 22:22:23.568034

#----------------------------------------------------------------------------#
## Class data                                                               ##
#----------------------------------------------------------------------------#

"""Read the input file and create the file necessary for Triangle, Tetgen and
the cross section file."""

import os
import utils

class data :
    """Read the input file given at the script and create the file necessary
    for Triangle, Tetgen and the cross section file."""
    
    def __init__(self,file_path) :

# We read the input file
        file = open(file_path,'r')
        self.input = file.read()
        
# We read the dimension of the problem
        self.read_dimension()

# We read the geometry
        self.read_geometry()

# We read the cross sections ids
#        self.read_xs_id()

        file.close()

#----------------------------------------------------------------------------#

    def search(self,keyword) :
        """Search for the given keyword. The beginning and the end of the 
        interesting section are saved in self.begin and self.end."""

# We search for the keyword in the file
        begin = "# BEGIN "
        end = "# END "

# We look for the beginning of the section        
        self.begin = self.input.find(begin+keyword) 
        if self.begin == -1 :
            utils.abort("Cannot find "+begin+keyword)
        else :
            self.begin += len(begin+keyword) + 1

# We look for the end of the section
        self.end = self.input.find(end+keyword)
        if self.end == -1 :
            utils.abort("Cannot find "+end+keyword)

# We look for an error
        error = self.input[self.begin:self.end-1].find(end)
        if error != -1 :
            utils.abort("Nested sections.")

#----------------------------------------------------------------------------#

    def read_next(self) :
        """Read the next value in the file. Return the value and a boolean to
        know if the number was correctly read."""

# We read the next element
        value = ' '

# We skip the blank space
        while self.input[self.pos] == " ":
            self.pos += 1
        
# We read the value        
        while self.input[self.pos] != " " : 
            if (self.input[self.pos] != "\n") :
                value += self.input[self.pos]
                self.pos += 1
            else :
                break
        if value != " " :
            value = float(value)
            number = True
        else :
            number = False

        return value,number

#----------------------------------------------------------------------------#

    def read_geometry(self) :
        """Read the geometry of the problem and put all the coordinates of the
        vertices of the polygons in self.polygons. The size of the cells is 
        computed here."""

# We read the geometry of the problem
        self.search("GEOMETRY")
        self.pos = self.begin
        self.polygons = []
        self.n_vertices = []
        while self.pos < self.end :         
            n_vertices,number = self.read_next()
            n_vertices = int(n_vertices)
            self.n_vertices.append(n_vertices)
            counter = 0.0
            polygon = []
            for i in xrange(0,n_vertices) :
                coord = []
                for j in xrange(0,self.dimension) :
                    value,number = self.read_next()
                    coord.append(value)
                if number :
                    counter += 1.0
                    polygon.append(coord)
                else :
                    utils.abort("Problem while reading the geometry.")
            self.polygons.append(polygon)
            self.pos += 1

# We read the number of cells that the user wants
        self.search("CELLS")
        self.pos = self.begin
        value,number = self.read_next()
        if number :
            self.n_cells = value
        else :
            utils.abort("Problem while reading the numbers of cells")

# We compute the maximum area of the cells.
        """The area of the domain is approximated by the area of the smaller
        rectangle which contains the domain."""
        left = 10^100
        right = -10^100
        bottom = 10e100
        top = -10^100
        
        self.n_polygons = len(self.polygons)
        for i in xrange(0,self.n_polygons) :
            for j in xrange(0,self.n_vertices[i]) :
                if self.polygons[i][j][0] < left :
                    left = self.polygons[i][j][0]
                if self.polygons[i][j][0] > right :
                    right = self.polygons[i][j][0]
                if self.polygons[i][j][1] < bottom :
                    bottom = self.polygons[i][j][1]
                if self.polygons[i][j][1] > top :
                    top = self.polygons[i][j][1]

# We convert the list in a float
        top = float(top)
        bottom = float(bottom)
        right = float(right)
        left = float(left)

# We compute the area
        self.area = (top-bottom)*(right-left)

        if self.dimension == 3 :
            down = 10^100
            up = -10^100
            for i in xrange(0,self.n_polygons) :
                for j in xrange(0,self.n_vertices[i]) :
                    if self.polygons[i][j][2] < down :
                        down = self.polygons[0][i][2]
                    if self.polygons[i][j][2] > up :
                        up = self.polygons[0][i][2]
            self.volume = self.area*(up-down)

        self.area_cells = self.area/self.n_cells

#----------------------------------------------------------------------------#
       
    def read_dimension(self) :
        """Read the dimension of the problem : 2D or 3D."""

# We read the dimension of the problem
        self.search("DIMENSION")
        self.pos = self.begin
        value,number = self.read_next()
        if number :
            if value == 2 or value == 3 :
                self.dimension = int(value)
        else :
            utils.abort("Wrong dimension")

#----------------------------------------------------------------------------#

    def read_xs_id(self) :
        """Read the cross sections ids. They are put in self.xs_id"""

# We read the cross section
        self.search("XS IDS")
        self.pos = self.begin
        self.xs_id = []
        for i in xrange(0,self.n_polygons) :
            value,number = self.read_next()
            if number :
                self.xs_id.append(value)
            else :
                utils.abort("Problem while reading the cross sections id")

#----------------------------------------------------------------------------#

    def create_mesh(self) :
        """Create the mesh using a mesh generator : Triangle for 2D and 
        Tetgen for 3D."""

        if self.dimension == 2 :
            self.generate_triangle_input()
            os.system("triangle -pqnea"+str(self.area_cells)+\
                    " mesh_input.poly")
        else :
            utils.abort("The mesh has to be 2D.")

#----------------------------------------------------------------------------#

    def generate_triangle_input(self) :
        """Generate the input files needed by Triangle."""

        total_n_vertices = 0.0
        for i in xrange(0,self.n_polygons) :
            total_n_vertices += self.n_vertices[i]
# Open the input file for triangle
        input_file = open("mesh_input.poly","w")
# Write the first line of the file : number of vertices
        line = str(int(total_n_vertices))
        line += " 2 0 0\n"
        input_file.write(line)
# Write the vertices
        vertex = 0
        for i in xrange(0,self.n_polygons) :
            for j in xrange(0,self.n_vertices[i]) :
                line = str(vertex)+" "+str(self.polygons[i][j][0])+" "+\
                        str(self.polygons[i][j][1])+" "+str(i)+\
                        " "+str(self.n_polygons)+"\n"
                vertex += 1
                input_file.write(line)
# Write the number of segments
        line = str(sum(self.n_vertices))+"\n"
        input_file.write(line)
# Write the segments
        segment = 0
        for i in xrange(0,self.n_polygons) :
            for j in xrange(0,self.n_vertices[i]-1) :
                line = str(segment)+" "+str(segment)+" "+str(segment+1)+"\n"
                segment += 1
                input_file.write(line)
            line = str(segment)+" "+str(segment)+" "+\
                    str(segment-self.n_vertices[i]+1)+"\n"
            input_file.write(line)
            segment += 1
# Write the number of holes
        input_file.write("0\n")

        input_file.close()
