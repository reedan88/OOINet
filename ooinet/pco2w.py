# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.6.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

from instrument import Instrument
import numpy as np
import xarray as xr


class PCO2W(Instrument):
    """Attributes and methods for reprocessing PCO2W data.
    
    This code is originally written by C. Wingard
    (chris.wingard@oregonstate.edu) and is adapted here.
    """
    
    def __init__(self):
        """Reformat PCO2W datasets downloaded from OOINet."""
        self.attrs = {
            'unique_id': {
                'long_name': 'Instrument Unique ID',
                'comment': ('One byte checksum summary of the instrument serial number, name, calibration date and firmware '
                            'version serving as a proxy for a unique ID. While identified as the instrument unique ID, it is '
                            'possible for more than one instrument to have the same checksum summary. Thus, the uniqueness '
                            'of this value should be considered with a certain degree of caution.'),
                # 'units': ''    # deliberately left blank, no units for this value
            },
            'duplicates': {
                'long_name': 'Duplicate Measurements Array',
                'comment': 'Dimensional indexing array created for the duplicate measurements collected during sampling.',
                # 'units': ''    # deliberately left blank, no units for this value
            },
            'dark_reference': {
                'long_name': 'Dark LED Reference Intensity',
                'comment': ('Dark LED reference intensity. Dark signal and reference intensities range between 50 and 200 '
                            'counts. Values outside of that range would indicate an issue with the instrument electronics. '
                            'Obtained from the light_measurements variable in the Data Portal sourced data file.'),
                'units': 'counts',
                'ooinet_variable_name': 'light_measurements'
            },
            'dark_signal': {
                'long_name': 'Dark LED Signal Intensity',
                'comment': ('Dark LED signal intensity. Dark signal and reference intensities range between 50 and 200 '
                            'counts. Values outside of that range would indicate an issue with the instrument electronics. '
                            'Obtained from the light_measurements variable in the Data Portal sourced data file.'),
                'units': 'counts',
                'ooinet_variable_name': 'light_measurements'
            },
            'reference_434': {
                'long_name': 'Reference Intensity at 434 nm',
                'comment': ('Optical absorbance reference intensity at 434 nm. Reference and signal intensities range '
                            'between 0 and 4096. Values should be greater than ~1500. Lower intensities will result in '
                            'higher noise in the absorbance and pCO2 measurements. Obtained from the light_measurements'
                            'variable in the Data Portal sourced data file.'),
                'units': 'counts',
                'ooinet_variable_name': 'light_measurements'
            },
            'signal_434': {
                'long_name': 'Signal Intensity at 434 nm',
                'comment': ('Optical absorbance signal intensity at 434 nm. Reference and signal intensities range '
                            'between 0 and 4096. Values should be greater than ~1500. Lower intensities will result in '
                            'higher noise in the absorbance and pCO2 measurements. Obtained from the light_measurements'
                            'variable in the Data Portal sourced data file.'),
                'units': 'counts',
                'ooinet_variable_name': 'light_measurements'
            },
            'reference_620': {
                'long_name': 'Reference Intensity at 620 nm',
                'comment': ('Optical absorbance reference intensity at 620 nm. Reference and signal intensities range '
                            'between 0 and 4096. Values should be greater than ~1500. Lower intensities will result in '
                            'higher noise in the absorbance and pCO2 measurements. Obtained from the light_measurements'
                            'variable in the Data Portal sourced data file.'),
                'units': 'counts',
                'ooinet_variable_name': 'light_measurements'
            },
            'signal_620': {
                'long_name': 'Signal Intensity at 620 nm',
                'comment': ('Optical absorbance signal intensity at 620 nm. Reference and signal intensities range between 0 '
                            'and 4096. Values should be greater than ~1500. Lower intensities will result in higher noise in '
                            'the absorbance and pCO2 measurements. Obtained from the light_measurements variable in the Data '
                            'Portal sourced data file.'),
                'units': 'counts',
                'ooinet_variable_name': 'light_measurements'
            },
            'absorbance_blank_434': {
                'long_name': 'Blank Optical Absorbance Ratio at 434 nm',
                'comment': ('The Optical Absorbance ratio at 434 nm collected during the blank cycle (measured against in-situ '
                            'seawater in the absence of reagent) and used to calculate the PCO2WAT data product. Values are '
                            'updated approximately every 2-3 days.'),
                'units': 'counts',
                'data_product_identifier': 'CO2ABS1-BLNK_L0',
                'stream': 'pco2w_abc_dcl_instrument_blank'
            },
            'absorbance_blank_620': {
                'long_name': 'Blank Optical Absorbance Ratio at 620 nm',
                'comment': ('The Optical Absorbance ratio at 620 nm collected during the blank cycle (measured against in-situ '
                            'seawater in the absence of reagent) and used to calculate the PCO2WAT data product. Values are '
                            'updated approximately every 2-3 days.'),
                'units': 'counts',
                'data_product_identifier': 'CO2ABS2-BLNK_L0',
                'stream': 'pco2w_abc_dcl_instrument_blank'
            },
            'absorbance_ratio_434': {
                'long_name': 'Optical Absorbance Ratio at 434 nm',
                'comment': 'The optical absorbance ratio at 434 nm collected during the measurement cycle.',
                'units': 'counts',
                'data_product_identifier': 'CO2ABS1-SAMP_L0'
            },
            'absorbance_ratio_620': {
                'long_name': 'Optical Absorbance Ratio at 620 nm',
                'comment': 'The optical absorbance ratio at 620 nm collected during the measurement cycle.',
                'units': 'counts',
                'data_product_identifier': 'CO2ABS2-SAMP_L0'
            },
            'raw_thermistor': {
                'long_name': 'Raw Thermistor Temperature',
                'comment': 'Thermistor resistivity measured in counts at the end of the measurement cycle.',
                'data_product_identifier': 'CO2THRM_L1',
                'units': 'counts'
            },
            'thermistor_temperature': {
                'long_name': 'Thermistor Temperature',
                'comment': ('Thermistor temperature refers to the internal instrument temperature of the pCO2 sensor, as '
                            'measured by the thermistor.'),
                'units': 'degrees_Celsius',
                'data_product_identifier': 'CO2THRM_L1',
                'ancillary_variables': 'raw_thermistor'
            },
            'raw_battery_voltage': {
                'long_name': 'Raw Battery Voltage',
                'comment': ('Raw internal battery voltage measured in counts. May actually reflect external voltage if '
                            'external power is applied'),
                'units': 'counts'
            },
            'battery_voltage': {
                'long_name': 'Battery Voltage',
                'comment': 'Internal battery voltage. May actually reflect external voltage if external power is applied',
                'units': 'V',
                'ancillary_variables': 'raw_battery_voltage'
            },
            'pco2_seawater': {
                'long_name': 'pCO2 Seawater',
                'standard_name': 'partial_pressure_of_carbon_dioxide_in_sea_water',
                'comment': ('Partial Pressure of CO2 in Seawater provides a measure of the amount of CO2 and HCO3 in seawater. '
                            'Specifically, it refers to the pressure that would be exerted by CO2 if all other gases were '
                            'removed. Partial pressure of a gas dissolved in seawater is understood as the partial pressure in '
                            'air that the gas would exert in a hypothetical air volume in equilibrium with that seawater. '
                            'NOTE: the following calibration coefficients in the parameterFunctionMap are deprecated: "ea434", '
                            '"eb434", "ea620", "eb620". Those arguments are still present in the function declaration, but '
                            'they are unused and can be supplied with an arbitrary fill value.'),
                'data_product_identifier': 'PCO2WAT_L1',
                'units': 'uatm',
                'ancillary_variables': ('absorbance_ratio_434 absorbance_blank_434 absorbance_ratio_620 absorbance_blank_620 '
                                        'thermistor_temperature'),
            },
        }
    
    def pco2w_instrument(self, ds):
        """
        Takes PCO2W data recorded by the instruments internally and cleans up the
        data set to make it more user-friendly. Primary task is renaming parameters
        and dropping some that are of limited use. Additionally, re-organize some
        of the variables to permit better assessments of the data.
        :param ds: initial PCO2W data set recorded by the instrument and
            downloaded from OOI via the M2M system
        :return: cleaned up and reorganized data set
        """
        # drop some of the variables:
        #   record_type == not used
        #   record_time == internal_timestamp == time, redundant so can remove
        #   absorbance_ratio_*_qc_results == incorrectly set tests, ignoring
        #   absorbance_ratio_*_qc_executed == incorrectly set tests, ignoring
        ds = ds.reset_coords()
        ds = ds.drop(['record_type', 'record_time', 'internal_timestamp',
                      'absorbance_ratio_434_qc_results', 'absorbance_ratio_434_qc_executed',
                      'absorbance_ratio_620_qc_results', 'absorbance_ratio_620_qc_executed'])

        # rename some of the variables for better clarity
        rename = {
            'voltage_battery': 'raw_battery_voltage',
            'thermistor_raw': 'raw_thermistor',
            'pco2w_thermistor_temperature': 'thermistor_temperature',
            'pco2w_thermistor_temperature_qc_executed': 'thermistor_temperature_qc_executed',
            'pco2w_thermistor_temperature_qc_results': 'thermistor_temperature_qc_results',
        }
        ds = ds.rename(rename)

        # now we need to reset the light array to named variables that will be more meaningful and useful in
        # the final data files
        light = ds.light_measurements.astype('int32')
        dark_reference = light[:, [0, 8]].values    # dark reference
        dark_signal = light[:, [1, 9]].values       # dark signal
        reference_434 = light[:, [2, 10]].values    # reference signal, 434 nm
        signal_434 = light[:, [3, 11]].values       # signal intensity, 434 nm
        reference_620 = light[:, [4, 12]].values    # reference signal, 620 nm
        signal_620 = light[:, [5, 13]].values       # signal intensity, 620 nm

        # create a data set with the duplicate measurements for each variable
        data = xr.Dataset({
            'dark_reference': (['time', 'duplicates'], dark_reference),
            'dark_signal': (['time', 'duplicates'], dark_signal),
            'reference_434': (['time', 'duplicates'], reference_434),
            'signal_434': (['time', 'duplicates'], signal_434),
            'reference_620': (['time', 'duplicates'], reference_620),
            'signal_620': (['time', 'duplicates'], signal_620)
        }, coords={'time': ds['time'],  'duplicates': np.arange(0, 2).astype('int32')})
        ds = ds.drop(['spectrum', 'light_measurements'])

        # merge the data sets back together
        ds = ds.merge(data)

        # calculate the battery voltage
        ds['battery_voltage'] = ds['raw_battery_voltage'] * 15. / 4096.

        # reset some of the data types
        data_types = ['deployment', 'raw_thermistor', 'raw_battery_voltage',
                      'absorbance_blank_434', 'absorbance_blank_620', 'absorbance_ratio_434',
                      'absorbance_ratio_620']
        for v in data_types:
            ds[v] = ds[v].astype('int32')

        data_types = ['thermistor_temperature', 'pco2_seawater']
        for v in data_types:
            ds[v] = ds[v].astype('float32')

        # reset some attributes
        for key, value in self.attrs.items():
            for atk, atv in value.items():
                if key in ds.variables:
                    ds[key].attrs[atk] = atv

        # add the original variable name as an attribute, if renamed
        for key, value in rename.items():
            ds[value].attrs['ooinet_variable_name'] = key

        return ds
    
    def quality_checks(self, ds):
        """
        Assessment of the raw data and the calculated seawater pCO2 for quality
        using a susbset of the QARTOD flags to indicate the quality. QARTOD
        flags used are:
            1 = Pass
            3 = Suspect or of High Interest
            4 = Fail
        Suspect flags are set based on ranges provided by the vendor. The final
        flag represents the worst case assessment of the data quality.
        
        Parameters
        ----------
        ds: (xarray.Dataset)
            An xarray dataset of the PCO2W data
            
        Returns
        -------
        ds: (xarray.Dataset)
            The PCO2W dataset returned with a new data variable 'qc_flags'
            which contains information about
        """
        qc_flag = ds.time.astype("int32") * 0 + 1

        # Check the dark reference & signal values
        refDark = (ds.dark_reference.mean(dim="duplicates") < 50) | (ds.dark_reference.mean(dim="duplicates") > 200)
        qc_flags[refDark] = 3
        sigDark = (ds.dark_signal.mean(dim="duplicates") < 50) | (ds.dark_signal.mean(dim="duplicates") > 200)
        qc_flags[sigDark] = 3

        # Check the 434 reference & signal values for suspect values
        ref434 = (ds.reference_434.mean(dim="duplicates") < 1500)
        qc_flags[ref434] = 3
        sig434 = (ds.signal_434.mean(dim="duplicates") < 1500)
        qc_flags[sig434] = 3

        # Check the 620 nm reference & signal values for suspect values
        ref620 = (ds.reference_620.mean(dim="duplicates") < 1500)
        qc_flags[ref620] = 3
        sig620 = (ds.signal_620.mean(dim="duplicates") < 1500)
        qc_flags[ref620] = 3

        # Check the 434 reference & signal values for bad values
        ref434 = (ds.reference_434.mean(dim="duplicates") < 0) | (ds.reference_434.mean(dim="duplicates") > 4096)
        qc_flags[ref434] = 4
        sig434 = (ds.signal_434.mean(dim="duplicates") < 0) | (ds.signal_434.mean(dim="duplicates") > 4096)
        qc_flags[sig434] = 4

        # Check the 620 reference & signal values for bad values
        ref620 = (ds.reference_620.mean(dim="duplicates") < 0) | (ds.reference_620.mean(dim="duplicates") > 4096)
        qc_flags[ref620] = 4
        sig620 = (ds.signal_620.mean(dim="duplicates") < 0) | (ds.signal_620.mean(dim="duplicates") > 4096)
        qc_flags[sig620] = 4
        
        # Add some attributes to the quality flags
        qc_flags.attrs = {
            "long_name": "Quality Flags",
            "comment": ("Assessment of the raw light intensity measurment data for quality using a " + 
                       "subset of the QARTOD flags to indicate quality. QARTOD flags used are: 1 = Pass; " +
                       "3 = Suspect; 4 = Fail. Suspect and Fail flags are determined using ranges provided " +
                       "by the vendor.")
        }
        
        # Now add it to the dataframe
        ds.qc_flags = qc_flags

        return ds   
        
