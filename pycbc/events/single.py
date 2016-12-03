""" utilities for assigning FAR to single detector triggers
"""

def live_single_far_line(ifo, triggers):
    from pycbc.events import newsnr
    from pycbc.io.live import SingleForGraceDB
    if len(triggers['snr']) == 0:
        return None 

    i = triggers['snr'].argmax()
    if triggers['snr'][i] > 13.0:
        d = {key: triggers[key][i] for key in triggers}
        d['stat'] = newsnr(triggers['snr'][i], triggers['chisq'][i])
        d['ifar'] = .01
        return SingleForGraceDB(ifo, d)
    else:
        return None
        
