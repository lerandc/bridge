import json

class EELS():
    """Helper class for writing EELS cards for FEFF simulation input files.

    Attributes:
        spectrum (str): "ELNES" or "EXELFS", required.
        energy (float): beam energy (mrad), required.
        alpha (float): convergence semiangle of beam (mrad), required.
        beta (float): collection semiangle of EELS detector (mrad), required.
        Nra (Tuple[int, int]): integration grid parameters, in number of rings
                            and points per ring. nr rings with na*(2*i-1) points
                            in i'th ring. optional, default (50,1).
        Dxy (Tuple[float, float]): position of detector in scattering plane, 
                                in mrad along x, y axes. optioanl, default (0.0, 0.0).
        K (Tuple[float, float, float]): wave vector of incoming electron in crystal 
                            frame. arbitrary units, only direction used. 
                            required if aver is false. 
        kmesh (Tuple[float, float, float]): maximum k value, 
                                            size of output k grid far from edge,
                                            energy step of grid at edge.
                                            defaults (5, 0.05, 0.1) for ELNES,
                                            (8, 0.07, 0.0) for EXELFS.
        aver (bool): whether or not to calculate orientation averaged spectrum. 
                    if false, must set K. optional, default True.
        cross (bool): whether or not to use cross terms for scattering cross sections.
                    if false, uses only direct terms. optional, default True.
        relat (bool): whether or not to use relativistic formula for cross-section
                    calculations. optional, default True.
    """

    def __init__(self, spectrum, energy, alpha, beta, **kwargs):
        self._spectrum = None
        self._energy = None
        self._alpha = None
        self._beta = None
        self._Nra = (50,1)
        self._Dxy = (0.0, 0.0)
        self._K = None
        self._aver = True
        self._cross = True
        self._relat = True

        self.spectrum = spectrum
        self.energy = energy
        self.alpha = alpha 
        self.beta = beta

        if self.spectrum == "ELNES":
            self.kmesh = (5, 0.05, 0.1)
        else:
            self.kmesh = (8, 0.07, 0.0)

        for key, value in kwargs.items():
            setattr(self, key, value)

    def __repr__(self):
        attrs = ["spectrum", "energy", "alpha", "beta", "Nra", "Dxy", "K", 
                "kmesh", "aver", "cross", "relat"]

        repr_dict = {}
        for a in attrs:
            repr_dict[a] = getattr(self, a)

        return json.dumps(repr_dict)

    def __str__(self):
        out_str = self.spectrum + " " \
                + str(self.kmesh[0]) + " " + str(self.kmesh[1]) + " " + str(self.kmesh[2]) + "\n" \
                + str(self.energy) + " " \
                + str(int(self.aver)) + " " + str(int(self.cross)) + " " + str(int(self.relat)) + "\n"

        if not self.aver:
            out_str += str(self.K[0]) + " " + str(self.K[1]) + " " + str(self.K[2]) + "\n"

        out_str += str(self.beta) + " " + str(self.alpha) + "\n" \
                + str(self.Nra[0]) + " " + str(self.Nra[1]) + "\n" \
                + str(self.Dxy[0]) + " " + str(self.Dxy[1]) + "\n"

        return out_str

    @property
    def spectrum(self):
        return self._spectrum

    @spectrum.setter
    def spectrum(self, val):
        assert((val=="ELNES") or (val=="EXELFS"))
        self._spectrum = val
    
    @property
    def energy(self):
        return self._energy

    @energy.setter
    def energy(self, val):
        assert(float(val) > 0.0)
        self._energy = float(val)
        
    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, val):
        assert(float(val) > 0.0)
        self._alpha = float(val)

    @property
    def beta(self):
        return self._beta

    @beta.setter
    def beta(self, val):
        assert(float(val) > 0.0)
        self._beta = float(val)

    @property
    def Nra(self):
        return self._Nra

    @Nra.setter
    def Nra(self, val):
        val = [int(x) for x in val]
        assert(len(tuple(val))==2), "Must provide iterable of only two values for Nra."
        self._Nra = val

    @property
    def Dxy(self):
        return self._Dxy

    @Dxy.setter
    def Dxy(self, val):
        val = [float(x) for x in val]
        assert(len(tuple(val))==2), "Must provide iterable of only two values for Dxy."
        self._Dxy = val

    @property
    def K(self):
        return self._K

    @K.setter
    def K(self, val):
        val = [float(x) for x in val]
        assert(len(tuple(val))==3), "Must provide iterable of only two values for K."
        self._K = val

    @property
    def kmesh(self):
        return self._kmesh

    @kmesh.setter
    def kmesh(self, val):
        val = [float(x) for x in val]
        assert(len(tuple(val))==3), "Must provide iterable of only two values for kmesh."
        self._kmesh = val

    @property
    def aver(self):
        return self._aver
    
    @aver.setter
    def aver(self, val):
        self._aver = bool(val)
        
    @property
    def cross(self):
        return self._cross
    
    @cross.setter
    def cross(self, val):
        self._cross = bool(val)

    @property
    def relat(self):
        return self._relat
    
    @relat.setter
    def relat(self, val):
        self._relat = bool(val)

