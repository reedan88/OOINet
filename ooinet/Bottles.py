import re
import numpy as np
import pandas as pd
from ooinet import haversine as hs

class QualityFlags():
    """QARTOD QC-flag definitions"""
    GOOD = 1
    UNKNOWN = 2
    SUSPECT = 3
    BAD = 4
    MISSING = 9

FLAGS = QualityFlags

def check_fill(flag):
    """Check if an OOI discrete sample flag is a fill value (-9999999)"""
    if pd.isna(flag):
        return True
    elif str(flag) == "-9999999":
        return True
    elif "1" not in str(flag):
        return True
    else:
        return False


def parse_flag(flag):
    """Parse an OOI discrete sample flag. Returns fill or nan when appropriate."""
    locs=[]
    for match in re.finditer("1", flag[::-1], re.S):
        locs.append(match.span()[0])
    return locs


def interp_ctd_flag(flag):
    """Interpret OOI discrete sample CTD flags to QARTOD QC-flags."""
    # First filter for fill
    if check_fill(flag):
        return flag
    else:
        parsed_flag = parse_flag(flag)
        max_bit = max(parsed_flag)
        if max_bit == 1:
            return QualityFlags.MISSING
        elif max_bit == 2:
            return QualityFlags.GOOD
        elif max_bit == 3:
            return QualityFlags.SUSPECT
        elif max_bit == 4:
            return QualityFlags.BAD
        else:
            return QualityFlags.UNKNOWN


def interp_discrete_flag(flag):
    """Interpret OOI discrete sample Bottle flags to QARTOD QC-flags."""
    # First filter for fill values
    if check_fill(flag):
        return flag
    else:
        parsed_flag = parse_flag(flag)
        max_bit = max(parsed_flag)
        if max_bit == 1:
            return QualityFlags.MISSING
        elif max_bit == 2:
            return QualityFlags.GOOD
        elif max_bit == 3:
            return QualityFlags.SUSPECT
        elif max_bit == 4:
            return QualityFlags.BAD
        else:
            return QualityFlags.UNKNOWN

def interp_replicate_flag(flag):
    """
    Interpret OOI discrete sample Bottle replicate flags. Returns a boolean
    if a sample is a duplicate/replicate sample.
    """
    # First filter for fill values
    if check_fill(flag):
        return flag
    else:
        parsed_flag = parse_flag(flag)
        max_bit = max(parsed_flag)
        if max_bit == 3 or max_bit == 4:
            return True
        else:
            return False


def interp_niskin_flag(flag):
    """Interpret OOI discrete Niskin Bottle flags to QARTOD QC-flags."""
    if check_fill(flag):
        return flag
    else:
        parsed_flag = parse_flag(flag)
        max_bit = max(parsed_flag)
        if max_bit == 1:
            return QualityFlags.MISSING
        elif max_bit == 2:
            return QualityFlags.GOOD
        elif max_bit == 3 or max_bit == 4 or max_bit == 5:
            return QualityFlags.SUSPECT
        else:
            return QualityFlags.UNKNOWN


def convert_times(x):
    """Parse times in OOI discrete summary spreadsheets to python datetimes."""
    if type(x) is str:
        x = x.replace(" ","")
        x = pd.to_datetime(x, utc=False)
    else:
        pass
    return x


def not_statistically_sigificant(x):
    """
    Replace OOI discrete nutrient sample values that are not statistically
    significant with zero.
    """
    if type(x) is str:
        if "<" in x:
            x = 0
    return x


def clean_data(bottleData):
    """
    Process, clean, and convert the OOI Discrete Sampling summary sheets.
    
    This function takes the Discrete Sample summary sheets provided by OOI
    and cleans up the spreadsheets, converts data types to be more useable,
    and intrepts the bit flag-maps into QARTOD-type flags.
    
    Parameters
    ----------
    bottleData: (pandas.DataFrame)
        A dataframe containing the loaded OOI Discrete Sampling summary data.
        
    Returns
    -------
    bottleData: (pandas.DataFrame)
        The discrete sampling data cleaned up and standardized.
    """
    # Replace -9999999 with NaNs
    bottleData = bottleData.replace(to_replace="-9999999", value=np.nan)
    bottleData = bottleData.replace(to_replace=-9999999, value=np.nan)
    
    # Convert times from strings to pandas datetime objects
    bottleData["Start Time [UTC]"] = bottleData["Start Time [UTC]"].apply(lambda x: convert_times(x))
    bottleData["CTD Bottle Closure Time [UTC]"] = bottleData["CTD Bottle Closure Time [UTC]"].apply(lambda x: convert_times(x))

    # Convert any values with a "<", which indicates a value not statistically significant from zero, with zero
    bottleData = bottleData.applymap(not_statistically_sigificant)
    
    # Interpret the quality flags to QARTOD flag values
    for col in bottleData.columns:
        if "Flag" in col:
            if "CTD" in col and "File" not in col:
                bottleData[col] = bottleData[col].apply(lambda x: interp_ctd_flag(x))
            elif "Discrete" in col:
                bottleData[col] = bottleData[col].apply(lambda x: interp_discrete_flag(x))
            elif "Replicate" in col:
                bottleData[col] = bottleData[col].apply(lambda x: interp_replicate_flag(x))
            elif "Niskin" in col:
                bottleData[col] = bottleData[col].apply(lambda x: interp_niskin_flag(x))
            else:
                pass
            
    return bottleData


def convert_oxygen(bottleData):
    """Convert oxygen from ml/l to umol/kg"""
    oxy = bottleData["Discrete Oxygen [mL/L]"]
    
    # Get relevant parameters
    SP = bottleData[["CTD Salinity 1 [psu]", "CTD Salinity 2 [psu]"]].mean(axis=1)
    T = bottleData[["CTD Temperature 1 [deg C]", "CTD Temperature 2 [deg C]"]].mean(axis=1)
    P = bottleData["CTD Pressure [db]"]
    LON = bottleData["Start Longitude [degrees]"]
    LAT = bottleData["Start Latitude [degrees]"]
    
    # Calculate absolute salinity & conservative temperature
    SA = gsw.SA_from_SP(SP, P, LON, LAT)
    CT = gsw.CT_from_t(SA, T, P)
    
    # Calculate Density
    RHO = gsw.rho(SA, CT, P)
    
    # Now convert from ml/l to umol/kg
    mvO2 = 22.392*1000 # ml/mole O2
    mole = (oxy/mvO2)*(1000/1)*(1/RHO)*1e6
    
    return mole


def findNearest(bottleData, buoyLoc, maxDist):
    """Find the bottle sample values within a maximum distance from the buoy
    
    Parameters
    ----------
    bottleData: (pd.DataFrame -> strings or floats)
        A tuple of (latitude, longitude) values in decimal degrees of the bottle sample location
    buoyLoc: (tuple -> floats)
        A tuple of (latitude, longitude) values in decimal degrees of the buoy location
    maxDist: (float)
        Maximum distance in km away for a sample location from the buoy location
    
    Returns
    -------
    mask: (boolean)
        Returns True or False boolean if sampleLoc < maxDist from buoyLoc
    """
    # Get the startLat/startLon as floats
    startLat = bottleData["Start Latitude [degrees]"].apply(lambda x: float(x))
    startLon = bottleData["Start Longitude [degrees]"].apply(lambda x: float(x))
    
    # Calculate the distance
    distance = []
    for lat, lon in zip(startLat, startLon):
        sampleLoc = (lat, lon)
        distance.append(hs.haversine(sampleLoc, buoyLoc))
    
    # Filter the results
    return [d <= maxDist for d in distance]


def findSamples(bottleData, buoyLoc, buoyDepth, maxDist, depthTol):
    """Find the bottle sample values within a maximum distance from the buoy
    
    Parameters
    ----------
    bottleData: (pd.DataFrame -> strings or floats)
        A tuple of (latitude, longitude) values in decimal degrees of the bottle sample location
    buoyLoc: (tuple -> floats)
        A tuple of (latitude, longitude) values in decimal degrees of the buoy location
    buoyDepth: (float)
        Deployment depth of the instrument
    maxDist: (float)
        Maximum distance in km away for a sample location from the buoy location
    depthTol: (float)
        Maximum depth difference to select samples from the buoyDepth
    
    Returns
    -------
    mask: (boolean)
        Returns True or False boolean if sampleLoc < maxDist from buoyLoc
    """
    # Filter for the nearest samples
    nearest = findNearest(bottleData, buoyLoc, maxDist)
    bottleData = bottleData[nearest]
    
    # Filter based on depth
    depthMin = buoyDepth - depthTol
    depthMax = buoyDepth + depthTol
    bottleData = bottleData[(bottleData["CTD Depth [m]"] >= depthMin) & (bottleData["CTD Depth [m]"] <= depthMax)]
    
    return bottleData
