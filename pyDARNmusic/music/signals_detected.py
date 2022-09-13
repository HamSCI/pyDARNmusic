from ..utils.musicUtils import stringify_signal_list

class SigDetect(object):
    """Class to hold information about detected signals.

    Methods
    -------
    string
    reorder

    Written by Nathaniel A. Frissell, Fall 2013
    Updated by: Francis Tholley, 2022
    """
    def __init__(self):
        pass
    def string(self):
        """Method to convert a list of signal dictionaries into strings.

        Returns
        -------
        stringInfo : list of str
            String representation of the signal information.

        Written by Nathaniel A. Frissell, Fall 2013
        Updated by: Francis Tholley, 2022
        """
        return stringify_signal_list(self.info)
    def reorder(self):
        """Method to sort items in .info by signal maximum value (from the scaled kArr) and update nrSignals.

        Written by Nathaniel A. Frissell, Fall 2013
        Updated by: Francis Tholley, 2022
        """
        #Do the sorting...
        from operator import itemgetter
        newlist = sorted(self.info,key=itemgetter('max'),reverse=True)

        #Put in the order numbers...
        order = 1
        for item in newlist:
            item['order'] = order
            order = order + 1

        #Save the list to the dataObj...
        self.info = newlist

        #Update the nrSigs
        self.nrSigs = len(newlist)

