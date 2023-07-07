import re
import numpy as np
import pandas as pd
import gsw
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


def clean_data(bottle_data):
    """
    Process, clean, and convert the OOI Discrete Sampling summary sheets.
    
    This function takes the Discrete Sample summary sheets provided by OOI
    and cleans up the spreadsheets, converts data types to be more useable,
    and intrepts the bit flag-maps into QARTOD-type flags.
    
    Parameters
    ----------
    bottle_data: (pandas.DataFrame)
        A dataframe containing the loaded OOI Discrete Sampling summary data.
        
    Returns
    -------
    bottle_data: (pandas.DataFrame)
        The discrete sampling data cleaned up and standardized.
    """
    # Replace -9999999 with NaNs
    bottle_data = bottle_data.replace(to_replace="-9999999", value=np.nan)
    bottle_data = bottle_data.replace(to_replace=-9999999, value=np.nan)
    
    # Convert times from strings to pandas datetime objects
    bottle_data["Start Time [UTC]"] = bottle_data["Start Time [UTC]"].apply(lambda x: convert_times(x))
    bottle_data["CTD Bottle Closure Time [UTC]"] = bottle_data["CTD Bottle Closure Time [UTC]"].apply(lambda x: convert_times(x))

    # Convert any values with a "<", which indicates a value not statistically significant from zero, with zero
    bottle_data = bottle_data.applymap(not_statistically_sigificant)
    
    # Interpret the quality flags to QARTOD flag values
    for col in bottle_data.columns:
        if "Flag" in col:
            if "CTD" in col and "File" not in col:
                bottle_data[col] = bottle_data[col].apply(lambda x: interp_ctd_flag(x))
            elif "Discrete" in col:
                bottle_data[col] = bottle_data[col].apply(lambda x: interp_discrete_flag(x))
            elif "Replicate" in col:
                bottle_data[col] = bottle_data[col].apply(lambda x: interp_replicate_flag(x))
            elif "Niskin" in col:
                bottle_data[col] = bottle_data[col].apply(lambda x: interp_niskin_flag(x))
            else:
                pass
            
    return bottle_data


def convert_oxygen(bottle_data):
    """Convert oxygen from ml/l to umol/kg"""
    oxy = bottle_data["Discrete Oxygen [mL/L]"]
    
    # Get relevant parameters
    SP = bottle_data[["CTD Salinity 1 [psu]", "CTD Salinity 2 [psu]"]].mean(axis=1)
    T = bottle_data[["CTD Temperature 1 [deg C]", "CTD Temperature 2 [deg C]"]].mean(axis=1)
    P = bottle_data["CTD Pressure [db]"]
    LON = bottle_data["Start Longitude [degrees]"]
    LAT = bottle_data["Start Latitude [degrees]"]
    
    # Calculate absolute salinity & conservative temperature
    SA = gsw.SA_from_SP(SP, P, LON, LAT)
    CT = gsw.CT_from_t(SA, T, P)
    
    # Calculate Density
    RHO = gsw.rho(SA, CT, P)
    
    # Now convert from ml/l to umol/kg
    mvO2 = 22.392*1000 # ml/mole O2
    mole = (oxy/mvO2)*(1000/1)*(1/RHO)*1e6
    
    return mole


def find_nearest(bottle_data, buoy_loc, max_dist):
    """Find the bottle sample values within a maximum distance from the buoy
    
    Parameters
    ----------
    bottle_data: (pd.DataFrame -> strings or floats)
        A tuple of (latitude, longitude) values in decimal degrees of the bottle sample location
    buoyLoc: (tuple -> floats)
        A tuple of (latitude, longitude) values in decimal degrees of the buoy location
    max_dist: (float)
        Maximum distance in km away for a sample location from the buoy location
    
    Returns
    -------
    mask: (boolean)
        Returns True or False boolean if sampleLoc < maxDist from buoyLoc
    """
    # Get the startLat/startLon as floats
    start_lat = bottle_data["Start Latitude [degrees]"].apply(lambda x: float(x))
    start_lon = bottle_data["Start Longitude [degrees]"].apply(lambda x: float(x))
    
    # Calculate the distance
    distance = []
    for lat, lon in zip(start_lat, start_lon):
        sample_loc = (lat, lon)
        distance.append(hs.haversine(sample_loc, buoy_loc))
    
    # Filter the results
    return [d <= max_dist for d in distance]


def find_samples(bottle_data, buoy_loc, buoy_depth, max_dist, depth_tol):
    """Find the bottle sample values within a maximum distance from the buoy
    
    Parameters
    ----------
    bottle_data: (pd.DataFrame -> strings or floats)
        A tuple of (latitude, longitude) values in decimal degrees of the bottle sample location
    buoy_loc: (tuple -> floats)
        A tuple of (latitude, longitude) values in decimal degrees of the buoy location
    buoy_depth: (float)
        Deployment depth of the instrument
    max_dist: (float)
        Maximum distance in km away for a sample location from the buoy location
    depth_tol: (float)
        Maximum depth difference to select samples from the buoy_depth
    
    Returns
    -------
    mask: (boolean)
        Returns True or False boolean if sampleLoc < maxDist from buoyLoc
    """
    # Filter for the nearest samples
    nearest = find_nearest(bottle_data, buoy_loc, max_dist)
    bottle_data = bottle_data[nearest]
    
    # Filter based on depth
    depth_min = buoy_depth - depth_tol
    depth_max = buoy_depth + depth_tol
    bottle_data = bottle_data[(bottle_data["CTD Depth [m]"] >= depth_min) & (bottle_data["CTD Depth [m]"] <= depth_max)]
    
    return bottle_data