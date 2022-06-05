def getParamDict(param):
    """Get information about a parameter, including units, default ranges,
    and axis labels.

    Parameters
    ----------
    param : str
        name of parameter

    Returns
    -------
    paramDict : str
        dictionary containing information about the chosen parameter

    Example
    -------
        paramDict = getParamDict('w_l')

    written by Nathaniel Frissell, 2013-07

    """
    import numpy as np

    # Create empty dictionary.
    paramDict = {}

    if param == 'p_l' or param == 'power':
        paramDict['param'] = 'power'
        paramDict['label'] = r'$\lambda$ Power'
        paramDict['unit'] = 'dB'
        paramDict['range'] = (0, 30)
    elif param == 'p_s':
        paramDict['param'] = 'power'
        paramDict['label'] = r'$\sigma$ Power'
        paramDict['unit'] = 'dB'
        paramDict['range'] = (0, 30)
    elif param == 'v' or param == 'velocity':
        paramDict['param'] = 'velocity'
        paramDict['label'] = 'Velocity'
        paramDict['unit'] = 'm s^{-1}'
        paramDict['range'] = (-500, 500)
    elif param.find('vheight') >= 0:
        paramDict['param'] = 'height'
        paramDict['label'] = "h'"
        paramDict['unit'] = 'km'
        paramDict['range'] = (75.0, 900.0)
    elif param == 'w_l' or param == 'width':
        paramDict['param'] = 'width'
        paramDict['label'] = r'$\lambda$ Spectral Width'
        paramDict['unit'] = 'm s^{-1}'
        paramDict['range'] = (0, 100)
    elif param == 'w_s':
        paramDict['param'] = 'width'
        paramDict['label'] = r'$\sigma$ Spectral Width'
        paramDict['unit'] = 'm s^{-1}'
        paramDict['range'] = (0, 100)
    elif param.find('elv') >= 0:
        paramDict['param'] = 'elevation'
        paramDict['label'] = 'Elevation'
        paramDict['unit'] = 'degrees'
        paramDict['range'] = (10, 30)
    elif param == 'phi0':
        paramDict['param'] = 'phi'
        paramDict['label'] = r'$\phi_0$'
        paramDict['unit'] = 'radians'
        paramDict['range'] = (-np.pi, np.pi)
    return paramDict
