import re
import datetime
import numpy as np
import pandas as pd
import xarray as xr
from pyOOI import M2M
from pyOOI.utils import convert_time, ntp_seconds_to_datetime, unix_epoch_time

def process_file(ds):
    """
    Function to download one of the NetCDF files as an xarray data set, convert
    to time as the appropriate dimension instead of obs, and drop the
    extraneous timestamp variables (these were originally not intended to be
    exposed to users and have lead to some confusion as to their meaning). The
    ID and provenance variables are better off obtained directly from the M2M
    system via a different process. Having them included imposes unnecessary
    constraints on the processing, so they are also removed.

    :param catalog_file: Unique file, referenced by a URL relative to the
        catalog, to download and then convert into an xarray data set.
    :param gc: Boolean flag to indicate whether the file is from the Gold
        Copy THREDDS server (gc = True), or the user's M2M THREDDS catalog
        (gc = False, default).
    :return: downloaded data in an xarray dataset.
    """
    # addresses error in how the *_qartod_executed variables are set
    qartod_pattern = re.compile(r'^.+_qartod_executed$')
    for v in ds.variables:
        if qartod_pattern.match(v):
            # the shape of the QARTOD executed variables should compare to the provenance variable
            if ds[v].shape != ds['provenance'].shape:
                ds = ds.drop_vars(v)

    # convert the dimensions from obs to time and get rid of obs and other variables we don't need
    ds = ds.swap_dims({'obs': 'time'})
    ds = ds.reset_coords()
    keys = ['obs', 'id', 'provenance', 'driver_timestamp', 'ingestion_timestamp',
            'port_timestamp', 'preferred_timestamp']
    for key in keys:
        if key in ds.variables:
            ds = ds.drop_vars(key)

    # since the CF decoding of the time is failing, explicitly reset all instances where the units are
    # seconds since 1900-01-01 to the correct CF units and convert the values to datetime64[ns] types
    time_pattern = re.compile(r'^seconds since 1900-01-01.*$')
    ntp_date = np.datetime64('1900-01-01')
    for v in ds.variables:
        if 'units' in ds[v].attrs.keys():
            if isinstance(ds[v].attrs['units'], str):  # because some units use non-standard characters...
                if time_pattern.match(ds[v].attrs['units']):
                    del(ds[v].attrs['_FillValue'])  # no fill values for time!
                    ds[v].attrs['units'] = 'seconds since 1900-01-01T00:00:00Z'
                    np_time = ntp_date + (ds[v] * 1e9).astype('timedelta64[ns]')
                    ds[v] = np_time

    # sort by time
    ds = ds.sortby('time')

    # clear-up some global attributes we will no longer be using
    keys = ['DODS.strlen', 'DODS.dimName',
            'DODS_EXTRA.Unlimited_Dimension', '_NCProperties', 'feature_Type']
    for key in keys:
        if key in ds.attrs:
            del(ds.attrs[key])

    if ds.encoding['unlimited_dims']:
        del ds.encoding['unlimited_dims']

    # resetting cdm_data_type from Point to Station and the featureType from point to timeSeries
    ds.attrs['cdm_data_type'] = 'Station'
    ds.attrs['featureType'] = 'timeSeries'

    # update some of the global attributes
    ds.attrs['acknowledgement'] = 'National Science Foundation'
    ds.attrs['comment'] = 'Data collected from the OOI M2M API and reworked for use in locally stored NetCDF files.'

    return ds


def add_annotation_qc_flag(ds, annotations):
    """Add the annotation qc flags to a dataset as a data variable.

    From the annotations, add the QARTOD flags to the dataset for
    each relevant data variable in the annotations.

    Parameters
    ----------
    ds: (xarray.DataSet)
        The xarray dataset containing the OOI data for a given
        reference designator-method-stream
    annotations: (pandas.DataFrame)
        A dataframe with contains the annotations to add to the
        dataset

    Returns
    -------
    ds: (xarray.DataSet)
        The input xarray dataset with the annotation qc flags
        added as a named variable to the dataset.
    """
    # First, filter only for annotations which apply to the dataset
    stream = ds.attrs["stream"]
    stream_mask = annotations["stream"].apply(lambda x: True if x == stream
                                              or x is None else False)
    annotations = annotations[stream_mask]

    # Second, explode the annotations so each parameter is hit for each
    # annotation
    annotations = annotations.explode(column="parameters")

    # Third, get the unique parameters and their associated variable name
    stream_annos = {}
    for pid in annotations["parameters"].unique():
        if np.isnan(pid):
            param_name = "rollup"
        else:
            param_info = M2M.get_api(M2M.URLS["preload"]
                                       + "/" + str(pid))
            param_name = param_info["name"]
        stream_annos.update({param_name: pid})

    # ----------------------------------------------------------------------
    # Next, get the flags associated with each parameter or all parameters
    flags_dict = {}

    for key in stream_annos.keys():
        # Get the pid and associated name
        pid_name = key
        pid = stream_annos.get(key)

        # Get the annotations associated with the pid
        if np.isnan(pid):
            pid_annos = annotations[annotations["parameters"].isna()]
        else:
            pid_annos = annotations[annotations["parameters"] == pid]

        pid_annos = pid_annos.sort_values(by="qcFlag")

        # Create an array of flags to begin setting the qc-values
        pid_flags = pd.Series(np.zeros(ds.time.values.shape),
                              index=ds.time.values)

        # For each index, set the qcFlag for each respective time period
        for ind in pid_annos.index:
            beginDT = pid_annos["beginDT"].loc[ind]
            endDT = pid_annos["endDT"].loc[ind]
            qcFlag = pid_annos["qcFlag"].loc[ind]
            # Convert the time to actual datetimes
            beginDT = convert_time(beginDT)
            if endDT is None or np.isnan(endDT):
                endDT = datetime.datetime.now()
            else:
                endDT = convert_time(endDT)
            # Set the qcFlags for the given time range
            pid_flags[(pid_flags.index > beginDT) & (pid_flags.index < endDT)] = qcFlag

        # Save the results
        flags_dict.update({pid_name: pid_flags})

    # --------------------
    # Create a rollup flag
    rollup_flags = flags_dict.get("rollup")
    for key in flags_dict:
        flags = np.max([rollup_flags, flags_dict.get(key)], axis=0)
        rollup_flags = pd.Series(flags, index=rollup_flags.index)
    # Replace the "All" with the rollup results
    flags_dict["rollup"] = rollup_flags

    # --------------------
    # Add the flag results to the dataset for key in flags_dict
    for key in flags_dict.keys():
        # Generate a variable name
        var_name = "_".join((key.lower(), "annotations", "qc", "results"))

        # Next, build the attributes dictionary
        if key.lower() == "rollup":
            comment = "These qc flags are a rollup summary which represents a Human-in-the-loop (HITL) assessment of the data quality for all applicable data variables in the dataset."
        else:
            comment = f"These qc flags represent a Human-in-the-loop (HITL) assessment of the data quality for the specific data variable {key}."
        long_name = f"{key} qc_flag"
        attrs = {
            "comment": comment,
            "long_name": long_name
        }

        # Now add to the dataset
        flags = xr.DataArray(flags_dict.get(key), dims="time", attrs=attrs)
        ds[var_name] = flags

    return ds