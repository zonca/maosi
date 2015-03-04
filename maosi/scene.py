
class Scene(object):
    """
    Allow the user to specify some scene consisting of a set of
    point sources with specifed pixel positions and fluxes.
    """

    def __init__(self, stars_x, stars_y, stars_f):
        """
        X Position (in pixels)
        Y Position (in pixels)
        Flux (in electrons/sec)
        """
        self.xpos = stars_x
        self.ypos = stars_y
        self.flux = stars_f

        return

