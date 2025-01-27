
from fipy import *
import pandas as pd

class CustomViewer (TSVViewer):
    
    def makedf(self, time, filename=None):
        """
        "plot" the coordinates and values of the variables to `filename`.
        If `filename` is not provided, "plots" to `stdout`.

        >>> from fipy.meshes import Grid1D
        >>> m = Grid1D(nx = 3, dx = 0.4)
        >>> from fipy.variables.cellVariable import CellVariable
        >>> v = CellVariable(mesh = m, name = "var", value = (0, 2, 5))
        >>> TSVViewer(vars = (v, v.grad)).plot() #doctest: +NORMALIZE_WHITESPACE
        x       var     var_gauss_grad_x
        0.2     0       2.5
        0.6     2       6.25
        1       5       3.75

        >>> from fipy.meshes import Grid2D
        >>> m = Grid2D(nx = 2, dx = .1, ny = 2, dy = 0.3)
        >>> v = CellVariable(mesh = m, name = "var", value = (0, 2, -2, 5))
        >>> TSVViewer(vars = (v, v.grad)).plot() #doctest: +NORMALIZE_WHITESPACE
        x       y       var     var_gauss_grad_x        var_gauss_grad_y
        0.05    0.15    0       10      -3.33333333333333
        0.15    0.15    2       10      5
        0.05    0.45    -2      35      -3.33333333333333
        0.15    0.45    5       35      5

        Parameters
        ----------
        filename : str
            If not `None`, the name of a file to save the image into.
        """

        mesh = self.vars[0].mesh
        dim = mesh.dim

        # neu
        
        headings = []
        
        for index in range(dim):
            headings.extend(self._axis[index])

        for var in self.vars:
            name = var.name
            if (isinstance(var, CellVariable) or isinstance(var, FaceVariable)) and var.rank == 1:
                for index in range(dim):
                    headings.extend(["%s_%s" % (name, self._axis[index])])
            else:
                headings.extend([name])
                
        
        # neu ende
        
        cellVars = [var for var in self.vars if isinstance(var, CellVariable)]
        faceVars = [var for var in self.vars if isinstance(var, FaceVariable)]

        if len(cellVars) > 0:
            values = mesh.cellCenters.globalValue
            for var in self.vars:
                if isinstance(var, CellVariable) and var.rank == 1:
                    values = numerix.concatenate((values, numerix.array(var.globalValue)))
                else:
                    values = numerix.concatenate((values, (numerix.array(var.globalValue),)))

            

        if len(faceVars) > 0:
            values = mesh.faceCenters.globalValue
            for var in self.vars:
                if isinstance(var, FaceVariable) and var.rank == 1:
                    values = numerix.concatenate((values, numerix.array(var.globalValue)))
                else:
                    values = numerix.concatenate((values, (numerix.array(var.globalValue),)))
                    
        df = pd.DataFrame(data= values.transpose(),columns=headings)
        
        df.insert(2,"time",time)

        return df

      