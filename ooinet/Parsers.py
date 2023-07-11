import os
import re
import csv
import pytz
import numpy as np
import pandas as pd
import datetime as dt

class METBK:
    """Contains parser for METBK raw data"""
    def __init__(self):
        self.DATA = {
            'TIMESTAMP': [],
            'BAROMETRIC_PRESSURE': [],
            'RELATIVE_HUMIDITY': [],
            'AIR_TEMPERATURE': [],
            'LONGWAVE_IRRADIANCE': [],
            'PRECIPITATION': [],
            'SEA_SURFACE_TEMPERATURE': [],
            'SEA_SURFACE_CONDUCTIVITY': [],
            'SHORTWAVE_IRRADIANCE': [],
            'WIND_EASTWARD': [],
            'WIND_NORTHWARD': [],
        }
        
        self.DATA_INDEX = {
            'TIMESTAMP': 0,
            'BAROMETRIC_PRESSURE': 1,
            'RELATIVE_HUMIDITY': 2,
            'AIR_TEMPERATURE': 3,
            'LONGWAVE_IRRADIANCE': 4,
            'PRECIPITATION': 5,
            'SEA_SURFACE_TEMPERATURE': 6,
            'SEA_SURFACE_CONDUCTIVITY': 7,
            'SHORTWAVE_IRRADIANCE': 8,
            'WIND_EASTWARD': 9,
            'WIND_NORTHWARD': 10,
        }
        
        self.DATA_PATTERN = (r'(-*\d+\.\d+|NaN)' +  # BPR
                       '\s*(-*\d+\.\d+|NaN)' +  # RH %
                       '\s*(-*\d+\.\d+|NaN)' +  # RH temp
                       '\s*(-*\d+\.\d+|NaN)' +  # LWR
                       '\s*(-*\d+\.\d+|NaN)' +  # PRC
                       '\s*(-*\d+\.\d+|NaN)' +  # ST
                       '\s*(-*\d+\.\d+|NaN)' +  # SC
                       '\s*(-*\d+\.\d+|NaN)' +  # SWR
                       '\s*(-*\d+\.\d+|NaN)' +  # We
                       '\s*(-*\d+\.\d+|NaN)' +  # Wn
                       '.*?' + '\n')  # throw away batteries

        self.SAMPLE_TIMESTAMP_PATTERN = (r'\d{4}/\d{2}/\d{2}' +  # Date in yyyy/mm/dd
                            '\s*\d{2}:\d{2}:\d{2}.\d+') # Time in HH:MM:SS.fff  
        
    def parse_data(self, raw_data):
        """
        Parses a line of METBK data into the individual sensor components

        Parameters
        ----------
        raw_data: str
            The opened data froma raw data file that has been read line-by-line
        
        Returns
        -------
        self.DATA: dict
            A dictionary of the parsed raw data stored into the applicable measurements
        """
        # Interate line-by-line through the raw data
        for line in raw_data:
            if line is not None:
                # Check that the line contains data
                try:
                    # If the last entry is a number, it is a line of data
                    float(line.split()[-1])
                    # Replace Na with NaN
                    line = re.sub(r'Na ', 'NaN ', line)
                    # Next, match the timestamp
                    timestamp = re.findall(self.SAMPLE_TIMESTAMP_PATTERN, line)
                    # Remove the timestamp from the data string
                    line = re.sub(timestamp[0], '', line)
                    # Get the data
                    raw_data = re.findall(self.DATA_PATTERN, line)[0]
                except:
                    # Check if there is a parsable timestamp
                    timestamp = re.findall(self.SAMPLE_TIMESTAMP_PATTERN, line)
                    if len(timestamp) != 0:
                        # Create an empty array of all NaNs
                        raw_data = ['NaN']*10
                    else:
                        # There is no useful info on the line
                        line = None
                
                # Append the timestamp to the start of the list
                if line is not None:
                    raw_data = list(raw_data)
                    raw_data.insert(0, timestamp[0])

                    # Now put the data into the data dictionary
                    for key in self.DATA.keys():
                        # Get the index of the data
                        index = self.DATA_INDEX.get(key)
                        self.DATA[key].append(raw_data[index])


    def process_data(self):
        """
        Process the parsed METBK data into a dataframe with derived variables
        
        This function takes in the parsed METBK data, converts it to a dataframe,
        indexes via time, converts data types from strings, and then resamples
        the data into 10 minute averages. Then, it derives the practical salinity
        from the conductivity and temperature, adjusts the barometric pressure to
        sea-surface-equivalent, and derives the absolute wind speed (m/s) and 
        wind direction (degrees) from the north and east wind vector components
        
        Parameters
        ----------
        self.METBK_DATA: dict
            The parsed METBK data stored in a dictionary
            
        Returns
        -------
        df: pandas.DataFrame
            A pandas DataFrame with the METBK data with the METBK_DATA dict keys as
            column headers, indexed by time, resampled to 10-minute averages, and 
            the derived variables.
        """
        # First, stick it into a dataframe
        df = pd.DataFrame(self.METBK_DATA)

        # Next, convert types
        df["TIMESTAMP"] = df["TIMESTAMP"].apply(lambda x: pd.to_datetime(x))
        df.set_index(keys=["TIMESTAMP"], inplace=True)
        df = df.applymap(float)

        # Bin into 10-minute increments
        df = df.resample('10T').mean()

        # Calculate practical salinity
        C = df["SEA_SURFACE_CONDUCTIVITY"]
        T = df["SEA_SURFACE_TEMPERATURE"]
        P = 1
        df["SEA_SURFACE_PRACTICAL_SALINITY"] = calculate_practical_salinity(C, T, P)

        # Adjust the barometric pressure to sea-level
        df["SEA_LEVEL_PRESSURE"] = adjust_pressure_to_sea_level(
            df["BAROMETRIC_PRESSURE"], 
            df["AIR_TEMPERATURE"], 
            4.05)

        # Calculate the wind speed
        df['WIND_SPEED'] = calculate_wind_speed(
            df['WIND_EASTWARD'],
            df['WIND_NORTHWARD'])

        # Calculate the wind direction
        df['WIND_DIRECTION'] = calculate_wind_direction(
            df['WIND_EASTWARD'],
            df['WIND_NORTHWARD'])
        
        # Adjust for wind directionsthat are outside 0-360
        df['WIND_DIRECTION'] = df["WIND_DIRECTION"].apply(
                        lambda x: x+360 if x < 0 else x)
        

        # Save the data into a new data structure
        self.DATAFRAME = df
    
    
class WAVSS():
    
    def __init__(self):
        
        self.DATA = {
            'TIMESTAMP': [],
            'INSTRUMENT_DATE': [],
            'INSTRUMENT_TIME': [],
            'INSTRUMENT_SERIAL': [],
            'BUOY_ID': [],
            'LATITUDE': [],
            'LONGITUDE': [],
            'N_ZERO_CROSSINGS': [],
            'AVERAGE_WAVE_HEIGHT': [],
            'MEAN_SPECTRAL_PERIOD': [],
            'MAXIMUM_WAVE_HEIGHT': [],
            'SIGNIFICANT_WAVE_HEIGHT': [],
            'SIGNIFICANT_PERIOD': [],
            'AVERAGE_HEIGHT_10TH_HIGHEST': [],
            'AVERAGE_PERIOD_10TH_HIGHEST': [],
            'MEAN_WAVE_PERIOD': [],
            'PEAK_PERIOD': [],
            'TP5': [],
            'HMO': [],
            'MEAN_DIRECTION': [],
            'MEAN_SPREAD': []
        }
        
        self.DATA_INDEX = {
            'TIMESTAMP': 0,
            'INSTRUMENT_DATE': 2,
            'INSTRUMENT_TIME': 3,
            'INSTRUMENT_SERIAL': 4,
            'BUOY_ID': 5,
            'LATITUDE': 6,
            'LONGITUDE': 7,
            'N_ZERO_CROSSINGS': 8,
            'AVERAGE_WAVE_HEIGHT': 9,
            'MEAN_SPECTRAL_PERIOD': 10,
            'MAXIMUM_WAVE_HEIGHT': 11,
            'SIGNIFICANT_WAVE_HEIGHT': 12,
            'SIGNIFICANT_PERIOD': 13,
            'AVERAGE_HEIGHT_10TH_HIGHEST': 14,
            'AVERAGE_PERIOD_10TH_HIGHEST': 15,
            'MEAN_WAVE_PERIOD': 16,
            'PEAK_PERIOD': 17,
            'TP5': 18,
            'HMO': 19,
            'MEAN_DIRECTION': 20,
            'MEAN_SPREAD': 21
        }
        
        
    
    def parse_data(self, raw_data):
        """
        Parse the raw_data into the different measurements
        
        Parameters
        ----------
        raw_data: str
            A string of each line of data from the instrument raw data file
            
        Returns
        -------
        self.DATA: dict
            A dictionary of the parsed raw data stored into the applicable
            measurements
            
        """
        for line in raw_data:
            
            # Check that its a wave_statistics measurement
            if '$TSPWA' not in line:
                continue

            # Dump everything after the "*"
            line = re.sub(r'\*.*', '', line, flags=re.DOTALL)

            # Split the data
            line = re.split(r' \$|,', line)

            # Check that it is a full data record. If not, return none
            if len(line) != 22:
                continue

            # Parse the raw data into the data dictionary based on index
            for key in self.DATA_INDEX.keys():
                # Get the index of a particular measurement
                index = self.DATA_INDEX.get(key)

                # Put the parsed raw data into the data dictionary
                self.DATA[key].append(line[index])
                
                
    def process_data(self):
        """
        Process the parsed WAVSS data into a dataframe resampled to 10-minute avg

        This function takes in the parsed WAVSS data, converts it to a dataframe,
        indexes via time, converts data types from strings, and then resamples
        the data into 10 minute averages. 

        Parameters
        ----------
        self.DATA: dict
            The parsed METBK data stored in a dictionary

        Returns
        -------
        self.DATAFRAME: pandas.DataFrame
            A pandas DataFrame with the METBK data with the METBK_DATA dict keys as
            column headers, indexed by time, resampled to 10-minute averages, and 
            the derived variables.
        """
        # First, stick it into a dataframe
        df = pd.DataFrame(self.DATA)
        df.drop(columns=["BUOY_ID", "LATITUDE", "LONGITUDE"], inplace=True)

        # Next, convert types
        df["TIMESTAMP"] = df["TIMESTAMP"].apply(lambda x: pd.to_datetime(x))
        df.set_index(keys=["TIMESTAMP"], inplace=True)
        df = df.applymap(float)

        # Bin into 10-minute increments
        df = df.resample('10T').mean()

        # Save the results to a dataframe
        self.DATAFRAME = df