from pydarn import SuperDARNRadars
def getRadarByCode(radCode=None):
    radarsInfo = SuperDARNRadars.radars.items()
    for key, radar in radarsInfo:
        if radCode == radar.hardware_info.abbrev:
            return key 

def getRadarById(radId=None):
    radarsInfo = SuperDARNRadars.radars.items()
    for key, radar in radarsInfo:
        if radId == radar.hardware_info.stid:
            return radar.hardware_info.abbrev

