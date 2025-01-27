from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.variables.cellToFaceVariable import _CellToFaceVariable
from fipy.tools import numerix
from fipy.variables.cellVariable import CellVariable

class _OneSidedCellToFaceVariable(_CellToFaceVariable):
    def _calcValue_(self, alpha, id1, id2):
        cell1 = numerix.take(self.var, id1, axis=-1)
        cell2 = numerix.take(self.var, id2, axis=-1)
        return numerix.maximum(cell1,cell2)
    
    


class ModifiedCellVariable(CellVariable):
    
    @property
    def onesidedFaceValue(self):

        if not hasattr(self, '_onesidedFaceValue'):
            # from fungal_interactions.alternativeFaceValue import _OneSidedCellToFaceVariable
            self._onesidedFaceValue = _OneSidedCellToFaceVariable(self)

        return self._onesidedFaceValue

    faceValue = onesidedFaceValue
    