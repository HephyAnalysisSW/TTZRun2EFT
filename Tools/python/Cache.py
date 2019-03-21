""" Just a wrapper for a safer way to use ResultsDB! Always use protection!
"""

from Analysis.Tools.ResultsDB import ResultsDB

class Cache:
    def __init__( self, filename, name, columns, overwrite=False ):
        self.columns  = sorted(columns)
        self.filename = filename
        self.DB       = ResultsDB( filename, name, self.columns )

    def translateKey( self, key ):
        key = { k:str( val ) for k, val in sorted(key.iteritems()) }
        return key

    def contains( self, key ):
        key = self.translateKey( key )
        return self.DB.contains( key )

    def get( self, key ):
        key = self.translateKey( key )
        return self.DB.get( key ).val

    def getDicts( self, key ):
        key = self.translateKey( key )
        return self.DB.getDicts( key )

    def add( self, key, val, overwrite=True ):
        key = self.translateKey( key )
        return self.DB.add( key, str(val), overwrite )

    def remove( self, key ):
        key = self.translateKey( key )
        return self.DB.removeObjects( key )

