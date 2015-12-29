
class Scene(object):
    """
    Allow the user to specify some scene consisting of a set of
    point sources with specifed positions in arcseconds and fluxes.
    """

    def __init__(self, stars_x, stars_y, stars_f):
        """
        X Position (in arcseconds)
        Y Position (in arcseconds)
        Flux (in electrons/sec)
        positions are in pa=0 with North up and East to the left.
        Positive x increases to the East.
        """
        self.xpos = stars_x
        self.ypos = stars_y
        self.flux = stars_f

        return

