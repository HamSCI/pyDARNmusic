class site:
    def __init__(self,radCode=None,dt=None):
        if dt is None or radCode is None:
            return
        import pydarn
        hdw_data    = pydarn.read_hdw_file(radCode,dt)    
        self.geolat = 0.0
        self.geolon = 0.0
        self.alt = 0.0
        self.boresite = 0.0
        self.bmsep = 0.0
        self.vdir = 0
        self.atten = 0.0
        self.tdiff = 0.0
        self.phidiff = 0.0
        self.interfer = [0.0, 0.0, 0.0]
        self.recrise = 0.0
        self.maxatten = 0
        self.maxgate = 0
        self.maxbeam = 0
        self.update_site_values(hdw_data)

    def __len__(self):
        """Object length """
        return 1

    def __str__(self):
        """Object string representation"""
        outstring = 'tval: {0} \
                    \ngeolat: {1:5.2f} \
                    \ngeolon: {2:5.2f} \
                    \nalt: {3:6.2f} \
                    \nboresite: {4:5.2f} \
                    \nbmsep: {5:5.2f} \
                    \nvdir: {6} \
                    \natten: {7:5.2f} \
                    \ntdiff: {8:6.4f} \
                    \nphidiff: {9:3.1f} \
                    \ninterfer: [{10:5.2f}, {11:5.2f}, {12:5.2f}] \
                    \nrecrise: {13:5.3f} \
                    \nmaxatten: {14} \
                    \nmaxgate: {15} \
                    \nmaxbeam: {16}'.format(self.tval, self.geolat,
                                            self.geolon, self.alt,
                                            self.boresite, self.bmsep,
                                            self.vdir, self.atten, self.tdiff,
                                            self.phidiff, self.interfer[0],
                                            self.interfer[1], self.interfer[2],
                                            self.recrise, self.maxatten,
                                            self.maxgate, self.maxbeam)
        return outstring


    def update_site_values(self, hdw_data):
        self.geolat = hdw_data.geographic.lat
        self.geolon = hdw_data.geographic.lon
        self.alt = hdw_data.geographic.alt
        self.boresite = hdw_data.boresight
        self.bmsep = hdw_data.beam_separation
        self.vdir = hdw_data.velocity_sign
        self.atten = hdw_data.rx_attenuator
        self.tdiff = hdw_data.tdiff
        self.phidiff = hdw_data.phase_sign
        self.__teminterfer = hdw_data.interferometer_offset
        self.interfer = [self.__teminterfer.x, self.__teminterfer.y, self.__teminterfer.z]
        self.recrise = hdw_data.rx_rise_time
        self.maxatten = hdw_data.attenuation_stages
        self.maxgate = hdw_data.gates
        self.maxbeam = hdw_data.beams