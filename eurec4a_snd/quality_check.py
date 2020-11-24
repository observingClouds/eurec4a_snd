"""
Quality check functions
"""
import numpy as np

def TU_sensor(snd, logging):
    """
    Meteomodem soundings have occasionally
    a mismatch between T, Td and RH
    """
    idx = np.where(snd.Temperature == snd.Dewpoint+273.15)
    _snd = snd.iloc[idx]
    idx_pd = _snd.index
    if np.any(snd.loc[idx_pd, 'Humidity'] != 100):
        logging.warning('Humidity mismatch, setting Td to nan')
        snd.loc[idx_pd, 'Dewpoint'] = np.nan
        # snd.loc[idx_pd, 'Temperature'] = np.nan
    return snd