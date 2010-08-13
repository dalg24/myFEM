# Python code
# Author: Bruno Turcksin
# Date: 2010-04-10 22:22:23.568034

#----------------------------------------------------------------------------#
## Class data                                                               ##
#----------------------------------------------------------------------------#

"""Read the input file and create the file necessary for Triangle, Tetgen and
the cross section file."""

import utils

class data :
    """Read the input file given at the script and create the file necessary
    for Triangle, Tetgen and the cross section file."""
    
    def __init__(self,file_path) :

# We read the input file
        file = open(file_path,'r')
        self.input = file.read()

# We read the geometry
        self.read_geometry()

# We read the dimension of the problem
        self.read_dimension()

        file.close()

#----------------------------------------------------------------------------#

    def search(self,keyword) :
        """Search for the given keyword. The beginning and the end of the 
        interesting section are saved in self.start and self.finish."""

# We search for the keyword in the file
        begin = "# BEGIN "
        end = "# END "

# We look for the beginning of the section        
        self.start = self.input.find(begin+keyword) 
        if self.start == -1 :
            utils.abort("Cannot find "+begin+keyword)
        else :
            self.start += len(begin+keyword) + 1

# We look for the end of the section
        self.finish = self.input.find(end+keyword)
        if self.finish == -1 :
            utils.abort("Cannot find "+end+keyword)

# We look for an error
        error = self.input[self.start:self.finish-1].find(end)
        if error != -1 :
            utils.abort("Nested sections.")

#----------------------------------------------------------------------------#

    def read_geometry(self) :
        """Read the geometry of the problem and put all the coordinates of the
        vertices of the polygons in self.polygons. The id corresponding to the
        polygons are put in self.xs_id. The size of the cells is computed here
        ."""

# We read the geometry of the problem
        self.search("GEOMETRY")
        self.pos = self.start
        self.polygons = []
        self.xs_id = []
        while self.pos < self.finish :         
            n_edges,number = self.read_next()
            counter = 0.0
            polygon = []
            while self.input[self.pos] != "\n" :
                value,number = self.read_next()
                if number :
                    if counter != n_edges :
                        counter += 1.0
                        polygon.append(value)
                    else :
                        self.xs_id.append(value)
                else :
                    utils.abort("Problem while reading the geometry.")
            if counter != n_edges :
                utils.abort("The number of edges is wrong.")
            else :
                self.polygons.append(polygon)
            self.pos += 1

# We read the number of cells that the user wants
        self.search("CELLS")
        self.pos = self.start
        value,number = self.read_next()
        if number :
            self.n_cells = value
        else :
            utils.abort("Problem while reading the numbers of cells")

# We compute the maximum area of the cells

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
       
    def read_dimension(self) :
        """Read the dimension of the problem : 2D or 3D."""

# read the dimension of the problem
        self.search("DIMENSION")
        self.pos = self.start
        value,number = self.read_next()
        if number :
            if value == 2 or value == 3 :
                self.dimension = value
        else :
            print ("Wrong dimension.")
            print ("The program aborted")
            exit()

#----------------------------------------------------------------------------#

    def generate_triangle_input(self) :
        """Generate the input files needed by Triangle."""
    
