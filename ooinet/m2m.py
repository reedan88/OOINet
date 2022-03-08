import re
import os
import time
import sys
import requests
import datetime
import numpy as np
import pandas as pd
import xarray as xr
from xml.dom import minidom
from urllib.request import urlopen
from urllib.request import urlretrieve
from urllib3.util import Retry
from requests.adapters import HTTPAdapter
from tqdm import tqdm
import dask
from dask.diagnostics import ProgressBar
from bs4 import BeautifulSoup


class M2M():
    """Request data from ooinet.oceanobservatories.org (OOINet) API by M2M.

    Attributes
    ----------
    username: (str)
        The OOINet supplied username. Requires registration. May be found in
        your user profile.
    token: (str)
        The OOINet supplied authentication token. Requires registration. May be
        found in your user profile.
    urls: (dict)
        Preset dictionary of urls for accessing data and metadata from OOINet.
    """

    def __init__(self, USERNAME, TOKEN):
        """Init M2M object.

        Params
        ------
        USERNAME: (str)
            The OOINet supplied username. Requires registration. May be found
            in your user profile.
        TOKEN: (str)
            The OOINet supplied authentication token. Requires registration.
            May be found in your user profile.
        """
        self.USERNAME = USERNAME
        self.TOKEN = TOKEN
        self.URLS = {
            'data': 'https://ooinet.oceanobservatories.org/api/m2m/12576/sensor/inv',
            'anno': 'https://ooinet.oceanobservatories.org/api/m2m/12580/anno/find',
            'vocab': 'https://ooinet.oceanobservatories.org/api/m2m/12586/vocab/inv',
            'asset': 'https://ooinet.oceanobservatories.org/api/m2m/12587',
            'deploy': 'https://ooinet.oceanobservatories.org/api/m2m/12587/events/deployment/inv',
            'preload': 'https://ooinet.oceanobservatories.org/api/m2m/12575/parameter',
            'cal': 'https://ooinet.oceanobservatories.org/api/m2m/12587/asset/cal',
            'fileServer':"https://opendap.oceanobservatories.org/thredds/fileServer/",
            'dodsC': "https://opendap.oceanobservatories.org/thredds/dodsC/",
            "goldCopy": "https://thredds.dataexplorer.oceanobservatories.org/thredds/catalog/ooigoldcopy/public/",
            "goldCopy_fileServer": "https://thredds.dataexplorer.oceanobservatories.org/thredds/fileServer/",
            "goldCopy_dodsC": "https://thredds.dataexplorer.oceanobservatories.org/thredds/dodsC/"
        }
        self.ENCODINGS = {
            "time": {"_FillValue": None},
            "lon": {"_FillValue": None},
            "lat": {"_FillValue": None},
            "z": {"_FillValue": None}
        }
        # Initialize a requests session
        self.SESSION = requests.Session()
        retry = Retry(connect=5, backoff_factor=0.5)
        adapter = HTTPAdapter(max_retries=retry)
        self.SESSION.mount("https://", adapter)
        
        # Set the REFDES attribute to None
        self.REFDES = None

    def _get_api(self, url, params=None):
        """Request the given url from OOINet."""
        r = self.SESSION.get(url, params=params, auth=(self.USERNAME, self.TOKEN))
        if r.status_code == requests.codes.ok:
            return r.json()
        else:
            return None

    def _ntp_seconds_to_datetime(self, ntp_seconds):
        """Convert OOINet timestamps to unix-convertable timestamps."""
        # Specify some constant needed for timestamp conversions
        ntp_epoch = datetime.datetime(1900, 1, 1)
        unix_epoch = datetime.datetime(1970, 1, 1)
        ntp_delta = (unix_epoch - ntp_epoch).total_seconds()

        return datetime.datetime.utcfromtimestamp(ntp_seconds - ntp_delta)

    def _convert_time(self, ms):
        if ms is None:
            return None
        else:
            return datetime.datetime.utcfromtimestamp(ms/1000)

    def _unix_epoch_time(self, date_time):
        """Convert a datetime to unix epoch microseconds."""
        # Convert the date time to a string
        date_time = int(pd.to_datetime(date_time).strftime("%s"))*1000
        return date_time

    def get_metadata(self, refdes):
        """Request metadata.

        Get the OOI Metadata for a specific instrument specified by its
        associated reference designator.

        Parameters
        ----------
        refdes: (str)
            OOINet standardized reference designator in the form of
            <array>-<node>-<instrument>.

        Returns
        -------
        results: (pandas.DataFrame)
            A dataframe with the relevant metadata of the given reference
            designator.
        """
        # First, construct the metadata request url
        array, node, instrument = refdes.split("-", 2)
        metadata_request_url = "/".join((self.URLS["data"], array, node,
                                         instrument, "metadata"))

        # Request the metadata
        metadata = self._get_api(metadata_request_url)

        # Parse the metadata
        metadata = self.parse_metadata(metadata)

        # Add in the reference designator
        metadata["refdes"] = refdes

        # Return the metadata
        return metadata

    def parse_metadata(self, metadata):
        """Parse metadata to dataframe.

        Parse the metadata dictionary for an instrument returned by OOI into
        a pandas dataframe.
        """
        # Put the two keys into separate dataframes
        metadata_times = pd.DataFrame(metadata["times"])
        metadata_parameters = pd.DataFrame(metadata["parameters"])

        # Merge the two into a single dataframe
        results = metadata_parameters.merge(metadata_times, left_on="stream",
                                            right_on="stream")
        results.drop_duplicates(inplace=True)

        # Return the results
        return results

    def get_deployments(self, refdes, deploy_num="-1", results=pd.DataFrame()):
        """Request deployment information for a reference designator.

        Get the deployment information for an instrument. Defaults to all
        deployments for a given instrument (reference designator) unless one is
        supplied.

        Parameters
        ----------
        refdes: (str)
            The reference designator for the instrument for which to request
            deployment information.
        deploy_num: (str)
            Optional to include a specific deployment number. Otherwise
            defaults to -1 which is all deployments.
        results: (pandas.DataFrame)
            Optional. Useful for recursive applications for gathering
            deployment information for multiple instruments.

        Returns
        -------
        results: (pandas.DataFrame)
            A table of the deployment information for the given instrument
            (reference designator) with deployment number, deployed water
            depth, latitude, longitude, start of deployment, end of deployment,
            and cruise IDs for the deployment and recovery.
        """
        # First, build the request
        array, node, instrument = refdes.split("-", 2)
        deploy_url = "/".join((self.URLS["deploy"], array, node, instrument,
                               deploy_num))

        # Next, get the deployments from the deploy url. The API returns a list
        # of dictionary objects with the deployment data.
        deployments = self._get_api(deploy_url)

        # Now, iterate over the deployment list and get the associated data for
        # each individual deployment
        while len(deployments) > 0:
            # Get a single deployment
            deployment = deployments.pop()

            # Process the dictionary data
            # Deployment Number
            deploymentNumber = deployment.get("deploymentNumber")

            # Location info
            location = deployment.get("location")
            depth = location["depth"]
            lat = location["latitude"]
            lon = location["longitude"]

            # Sensor info
            sensor = deployment.get("sensor")
            uid = sensor["uid"]
            assetId = sensor["assetId"]

            # Start and end times of the deployments
            startTime = self._convert_time(deployment.get("eventStartTime"))
            stopTime = self._convert_time(deployment.get("eventStopTime"))

            # Cruise IDs of the deployment and recover cruises
            deployCruiseInfo = deployment.get("deployCruiseInfo")
            recoverCruiseInfo = deployment.get("recoverCruiseInfo")
            if deployCruiseInfo is not None:
                deployID = deployCruiseInfo["uniqueCruiseIdentifier"]
            else:
                deployID = None
            if recoverCruiseInfo is not None:
                recoverID = recoverCruiseInfo["uniqueCruiseIdentifier"]
            else:
                recoverID = None

            # Put the data into a pandas dataframe
            data = np.array([[deploymentNumber, uid, assetId, lat, lon, depth,
                              startTime, stopTime, deployID, recoverID]])
            columns = ["deploymentNumber", "uid", "assetId",  "latitude",
                       "longitude", "depth", "deployStart", "deployEnd",
                       "deployCruise", "recoverCruise"]
            df = pd.DataFrame(data=data, columns=columns)

            # Generate the table results
            results = results.append(df)

        return results
    
    def _reformat_calInfo(calInfo, deployNum, uid):
        """Reformats calibration information returned from OOINet."""
        calDataFrame = pd.DataFrame()
        # Reformat the calibration information
        for cal in calInfo["calibration"]:
            for calData in cal["calData"]:
                calDataFrame = calDataFrame.append({
                    "deploymentNumber": int(deployNum),
                    "uid": uid,
                    "calCoef": calData["eventName"],
                    "calDate": self._convert_time(calData["eventStartTime"]),
                    "value": calData["value"],
                    "calFile": calData["dataSource"]
                }, ignore_index=True)
                
        return calDataFrame
    
    def _reformat_calInfo(self, calInfo, deployNum, uid):
        """Internal method to reformat calibration info"""
        calDataFrame = pd.DataFrame()
        # Reformat the calibration info
        for cal in calInfo["calibration"]:
            for calData in cal["calData"]:
                calDataFrame = calDataFrame.append({
                    "deploymentNumber": int(deployNum),
                    "uid": uid,
                    "calCoef": calData["eventName"],
                    "value": calData["value"],
                    "calFile": calData["dataSource"]
                }, ignore_index=True)
        return calDataFrame  
    
    def get_calibrations(self, refdes, deployments):
        """Get calibrations for deployments for a given reference designator.
        
        Return the effective calibrations for a given reference designator
        for all of the deployments.
        
        Parameters
        ----------
        refdes: (str)
            The reference designator for the instrument for which to request
            vocab information.
        deployments: (pandas.DataFrame)
            A dataframe with the given deployments for a reference designator
            as returned by the get_deployments method.
            
        Returns
        -------
        calibrations: (pandas.DataFrame)
            A table of the calibration information for the given reference
            designator for all of the deployments specified by the 
            deployments dataframe.
        """
        # Create a dask delayed object to gather the calibrations for the
        # instrument from each deployment
        cals = []
        for ind, row in deployments.iterrows():

            # Get the relevant request data
            deployNum, uid, deployStart, deployEnd = row["deploymentNumber"], row["uid"], row["deployStart"], row["deployEnd"]

            # Reformat the deployStart and deployEnd
            if pd.isna(deployEnd):
                deployEnd = datetime.datetime.now()

            # Construct the request parameters
            params = {
                "refdes": refdes,
                "uid": uid,
                "beginDT": deployStart.strftime("%Y-%m-%dT%H:%M:%S.000Z"),
                "endDT": deployEnd.strftime("%Y-%m-%dT%H:%M:%S.000Z")
            }

            # Request the data
            calInfo = dask.delayed(self._get_api)(self.URLS["cal"], params)

            # Reformat the data
            calDataFrame =  dask.delayed(self._reformat_calInfo)(calInfo, deployNum, uid)

            cals.append(calDataFrame)
        
        # Execute and get the calibrations for each deployment
        with ProgressBar():
            print(f"Fetching calibrations for {refdes}:")
            cals = dask.compute(cals)

        # Append individual calibrations into a single dataframe
        calibrations = pd.DataFrame()
        for cal in cals:
            calibrations = calibrations.append(cal, ignore_index=True)
        return calibrations

    def get_vocab(self, refdes):
        """Get OOI vocabulary.

        Return the OOI vocabulary for a given url endpoint. The vocab results
        contains info about the reference designator, names of the

        Parameters
        ----------
        refdes: (str)
            The reference designator for the instrument for which to request
            vocab information.

        Returns
        -------
        results: (pandas.DataFrame)
            A table of the vocab information for the given reference
            designator.
        """
        # First, construct the vocab request url
        array, node, instrument = refdes.split("-", 2)
        vocab_url = "/".join((self.URLS["vocab"], array, node, instrument))

        # Next, get the vocab data
        data = self._get_api(vocab_url)

        # Put the returned vocab data into a pandas dataframe
        vocab = pd.DataFrame()
        vocab = vocab.append(data)

        # Finally, return the results
        return vocab

    def get_datasets(self, search_url, datasets=pd.DataFrame(), **kwargs):
        """Search OOINet for available datasets for a url."""
        # Check if the method is attached to the url
        flag = ("inv" == search_url.split("/")[-4])
        # inst = re.search("[0-9]{2}-[023A-Z]{6}[0-9]{3}", search_url)
        # inst = re.search("[0-9]{2}-", search_url)

        # This means you are at the end-point
        if flag is True:
            # Get the reference designator info
            array, node, instrument = search_url.split("/")[-3:]
            refdes = "-".join((array, node, instrument))

            # Get the available deployments
            deploy_url = "/".join((self.URLS["deploy"], array, node,
                                   instrument))
            deployments = self._get_api(deploy_url)

            # Put the data into a dictionary
            info = pd.DataFrame(data=np.array([[array, node, instrument,
                                                refdes, search_url, deployments
                                                ]]),
                                columns=["array", "node", "instrument",
                                         "refdes", "url", "deployments"])
            # add the dictionary to the dataframe
            datasets = datasets.append(info, ignore_index=True)

        else:
            endpoints = self._get_api(search_url)

            while len(endpoints) > 0:

                # Get one endpoint
                new_endpoint = endpoints.pop()

                # Build the new request url
                new_search_url = "/".join((search_url, new_endpoint))

                # Get the datasets for the new given endpoint
                datasets = self.get_datasets(new_search_url, datasets)

        # Once recursion is done, return the datasets
        return datasets

    def search_datasets(self, array=None, node=None, instrument=None,
                        English_names=False):
        """Search OOINet for datasets.

        Parameters
        ----------
        array: (str)
            OOI abbreviation for a particular buoy on an array (e.g. Pioneer
            Central Surface Mooring = CP01CNSM)
        node: (str)
            Partial or full OOI abbreviation for a node on a buoy to search for
            (e.g. Multi-Function Node = MFD)
        instrument: (str)
            Partial or full OOI abbreviation for a particular instrument type
            to search for (e.g. CTD)
        English_names: (bool)
            Set to True if the descriptive names associated with the given
            array/node/instrument are wanted.

        Returns
        -------
        datasets: (pandas.DataFrame)
            A dataframe of all the OOI datasets which match the given search
            terms. If no search terms are entered, will return every dataset
            available in OOINet (slow).
        """
        # Build the request url
        dataset_url = f'{self.URLS["data"]}/{array}/{node}/{instrument}'

        # Truncate the url at the first "none"
        dataset_url = dataset_url[:dataset_url.find("None")-1]

        print(dataset_url)
        # Get the datasets
        datasets = self.get_datasets(dataset_url)

        # Now, it node is not None, can filter on that
        if node is not None:
            mask = datasets["node"].apply(lambda x: True if node
                                          in x else False)
            datasets = datasets[mask]

        # If instrument is not None
        if instrument is not None:
            mask = datasets["instrument"].apply(lambda x: True if instrument
                                                in x else False)
            datasets = datasets[mask]

        # Check if they want the English names for the associated datasets
        if English_names:
            vocab = {
                "refdes": [],
                "array_name": [],
                "node_name": [],
                "instrument_name": []
            }

            # Iterate through the given reference designators
            for refdes in datasets["refdes"]:
                # Request the vocab for the given reference designator
                refdes_vocab = self.get_vocab(refdes)

                # Check if it returns an empty dataframe - then fill with NaNs
                if len(refdes_vocab) == 0:
                    vocab["refdes"].append(refdes)
                    vocab["array_name"].append("No record")
                    vocab["node_name"].append("No record")
                    vocab["instrument_name"].append("No record")

                # If it isn't empty - Parse the refdes-specific vocab
                else:
                    vocab["refdes"].append(refdes)
                    vocab["array_name"].append(
                        refdes_vocab["tocL1"].iloc[0] + " " +
                        refdes_vocab["tocL2"].iloc[0])
                    vocab["node_name"].append(refdes_vocab["tocL3"].iloc[0])
                    vocab["instrument_name"].append(
                        refdes_vocab["instrument"].iloc[0])

            # Merge the results with the datasets
            vocab = pd.DataFrame(vocab)
            datasets = datasets.merge(vocab, left_on="refdes",
                                      right_on="refdes")
            # Sort the datasets
            columns = ["array", "array_name", "node", "node_name",
                       "instrument", "instrument_name", "refdes", "url",
                       "deployments"]
            datasets = datasets[columns]

        return datasets

    def get_datastreams(self, refdes):
        """Retrieve methods and data streams for a reference designator."""
        # Build the url
        array, node, instrument = refdes.split("-", 2)
        method_url = "/".join((self.URLS["data"], array, node, instrument))

        # Build a table linking the reference designators, methods, and data
        # streams
        stream_df = pd.DataFrame(columns=["refdes", "method", "stream"])
        methods = self._get_api(method_url)
        for method in methods:
            if "bad" in method:
                continue
            stream_url = "/".join((method_url, method))
            streams = self._get_api(stream_url)
            stream_df = stream_df.append({
                "refdes": refdes,
                "method": method,
                "stream": streams
            }, ignore_index=True)

        # Expand so that each row of the dataframe is unique
        stream_df = stream_df.explode('stream').reset_index(drop=True)

        # Return the results
        return stream_df

    def get_annotations(self, refdes, **kwargs):
        """Retrieve data annotations for a given reference designator.

        Parameters
        ----------
        refdes: (str)
            The reference designator which to query OOINet for the
            associated annotation data.

        Kwargs
        ------
        method: (str)
            An OOINet method associated with the reference designator.
            Limits annotations for the given reference designator
            to only that method.
        stream: (str)
            An OOINet stream associated with the reference designator.
            Limits annotations for the given reference designator
            to only that stream.
        beginDT: (str)
            Limit the data request to only data after this date. Date
            should be formatted using OOINet epoch time (can be
            calculated using _)
        endDT: (str)
            Limit the data request to only data before this date.
        """
        # Need to build a parameters dictionary to pass to requests
        params = {"refdes": refdes}
        for key in kwargs:
            val = kwargs.get(key)
            # Convert datetimes to unix epoch
            if key == "beginDT":
                val = self._unix_epoch_time(val)
            elif key == "endDT":
                val = self._unix_epoch_time(val)
            else:
                pass
            params.update({key: val})

        # Get the annotations as a json and put into a dataframe
        anno_data = self._get_api(self.URLS["anno"], params=params)
        anno_data = pd.DataFrame(anno_data)

        # Convert the flags to QARTOD flags
        codes = {
            None: 0,
            'pass': 1,
            'not_evaluated': 2,
            'suspect': 3,
            'fail': 4,
            'not_operational': 9,
            'not_available': 9,
            'pending_ingest': 9
        }
        anno_data['qcFlag'] = anno_data['qcFlag'].map(codes).astype('category')

        return anno_data

    def add_annotation_qc_flag(self, ds, annotations):
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
                param_info = self._get_api(self.URLS["preload"]
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
                beginDT = self._convert_time(beginDT)
                if endDT is None or np.isnan(endDT):
                    endDT = datetime.datetime.now()
                else:
                    endDT = self._convert_time(endDT)
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

    def get_parameter_data_levels(self, metadata):
        """Get parameters processing levels.

        Get the data levels (processing level) associated with the parameters
        for a given reference designator.

        Parameters
        ----------
        metadata: (pandas.DataFrame)
            The dataframe of metadata returned by get_metadata which contains
            the metadata for a given reference designator.

        Returns
        -------
        pid_dict: (dict)
            A dictionary with the data levels for each parameter id (Pid)
        """
        pdIds = np.unique(metadata["pdId"])
        pid_dict = {}
        for pid in pdIds:
            # Build the preload url
            preload_url = "/".join((self.URLS["preload"], pid.strip("PD")))
            # Query the preload data
            preload_data = self._get_api(preload_url)
            data_level = preload_data.get("data_level")
            # Update the results dictionary
            pid_dict.update({pid: data_level})

        return pid_dict

    def filter_parameter_ids(self, pdId, pid_dict):
        """Filter for processed data products."""
        # Check if pdId should be kept
        data_level = pid_dict.get(pdId)
        if data_level == 1:
            return True
        else:
            return False

    def get_thredds_url(self, refdes, method, stream, goldCopy=False, **kwargs):
        """
        Return the url for the THREDDS server for the desired dataset(s).

        Parameters
        ----------
        refdes: (str)
            Reference designator for the instrument
        method: (str)
            The method (i.e. telemetered) for the given reference designator
        stream: (str)
            The stream associated with the reference designator and method

        Kwargs
        ------
        beginDT: (str)
            Limit the data request to only data after this date.
        endDT: (str)
            Limit the data request to only data before this date.
        format: (str)
            e.g. "application/netcdf" (the default)
        include_provenance (str):
            'true' returns a text file with the provenance information
        include_annotations: (str)
            'true' returns a separate text file with annotations for the data
            within the given date range

        Returns
        -------
        thredds_url: (str)
            A url to the OOI Thredds server which contains the desired datasets
        """
        # Set the refdes as an attribute
        if self.REFDES is None:
            self.REFDES = refdes

        # Build the data request url
        if goldCopy is True:
            goldCopy_url = self.URLS["goldCopy"] + "-".join((self.REFDES, method, stream)) + "/catalog.html"
            return goldCopy_url
        else:
            array, node, instrument = refdes.split("-", 2)
            data_request_url = "/".join((self.URLS["data"], array, node,
                                     instrument, method, stream))

        # Ensure proper datetime format for the request
        if 'beginDT' in kwargs.keys():
            kwargs['beginDT'] = pd.to_datetime(kwargs['beginDT']).strftime(
                '%Y-%m-%dT%H:%M:%S.%fZ')
        if 'endDT' in kwargs.keys():
            kwargs['endDT'] = pd.to_datetime(kwargs['endDT']).strftime(
                '%Y-%m-%dT%H:%M:%S.%fZ')

        # Build the query
        if len(kwargs) == 0:
            kwargs = {"require_deployment":True}
        else:
            kwargs.update({
                "require_deployment":True
            })
        params = kwargs

        # Request the data
        r = self.SESSION.get(data_request_url, params=params, auth=(self.USERNAME, self.TOKEN))
        
        if r.status_code == requests.codes.ok:
            data_urls = r.json()
        else:
            print(r.reason)
            return None
        
        # NEW: Wait until OOINet has finished processing to fulfill the request
        print(f"Waiting for {refdes}-{method}-{stream} to process.")
        status_url = [url for url in data_urls["allURLs"] if re.match(r'.*async_results.*', url)][0]
        status_url = status_url + "/status.txt"
        with tqdm(total=400, desc="Waiting", file=sys.stdout) as bar:
            for i in range(400):
                r = self.SESSION.get(status_url)
                bar.update()
                bar.refresh()
                if r.status_code == requests.codes.ok:
                    bar.n = 400
                    bar.last_print_n = 400
                    break
                else:
                    time.sleep(3)

        # The asynchronous data request is contained in the 'allURLs' key,
        # in which we want to find the url to the thredds server
        for d in data_urls['allURLs']:
            if 'thredds' in d:
                thredds_url = d

        return thredds_url

    def get_thredds_catalog(self, thredds_url):
        """
        Get the dataset catalog for the requested data stream.

        Parameters
        ----------
        thredds_url (str): the THREDDS server url for the
            requested data stream

        Returns
        -------
        catalog (list): the THREDDS catalog of datasets for
            the requested data stream
        """
        # ==========================================================
        # Parse out the dataset_id from the thredds url
        page = self.SESSION.get(thredds_url).text
        soup = BeautifulSoup(page, "html.parser")
        pattern = re.compile('.*\\.nc$')
        catalog = sorted([node.get('href') for node in soup.find_all('a', text=pattern)])

        return catalog

    def parse_catalog(self, catalog, exclude=[]):
        """Parse THREDSS catalog for netCDF files.

        Parses the THREDDS catalog for the netCDF files. The exclude
        argument takes in a list of strings to check a given catalog
        item against and, if in the item, not return it.

        Parameters
        ----------
        catalog: (list)
            The THREDDS catalog of datasets for the requested data stream
        exclude: (list)
            Keywords to filter files out of the THEDDS catalog

        Returns
        -------
        datasets: (list)
            A list of netCDF datasets which contain the associated .nc datasets
        """
        datasets = [citem for citem in catalog if citem.endswith('.nc')]
        if type(exclude) is not list:
            raise ValueError('arg exclude must be a list')
        for ex in exclude:
            if type(ex) is not str:
                raise ValueError(f'Element {ex} of exclude must be a string.')
            datasets = [dset for dset in datasets if ex not in dset]
        return datasets

    def download_netCDF_files(self, catalog, goldCopy=False, saveDir=None):
        """Download netCDF files for given netCDF datasets.

        Downloads the netCDF files returned by parse_catalog. If no path is
        specified for the save directory, will download the files to the
        current working directory.

        Parameters
        ----------
        datasets: (list)
            The netCDF datasets to download
        saveDir: (str)
            The full path to the directory where to download the netCDF files.
        goldCopy: (boolean)
            If you should use the gold copy thredds server (which is static)
        """
        ## Step 1 - parse the catalog to get the correct web address to download the files
        if goldCopy is True:
            goldCopy = self.URLS["goldCopy_fileServer"]
        else:
            fileServer = self.URLS["fileServer"]
        netCDF_files = [re.sub("catalog.html\?dataset=", fileServer, file) for file in catalog]

        # Specify and make the relevant save directory
        if saveDir is not None:
            # Make the save directory if it doesn't exists
            if not os.path.exists(saveDir):
                os.makedirs(saveDir)
        else:
            saveDir = os.getcwd()

        # Download and save the netCDF files from the HTTPServer
        # to the save directory
        download = []
        for file in netCDF_files:
            filename = file.split("/")[-1]
            download.append(dask.delayed(urlretrieve)(file, "/".join((saveDir, filename))))
            
        with ProgressBar():
            print(f"Downloading files to {saveDir}")
            dask.compute(download)


    def _check_files(self, netCDF_files):
        """Internal method checks netCDF files if they are bad or empty."""
        
        def check_length(data, filename):
            length = len(data["time"])
            if length > 2:
                return filename
            else:
                pass

        # Okay, try redoing this approach to not break
        lengths = []
        for file in netCDF_files:
            # First, open the file
            data = dask.delayed(xr.open_dataset)(file)
            # Second, check the length of the primary index
            length = dask.delayed(check_length)(data, file)
            lengths.append(length)
        
        print("Checking and removing bad files: ")
        with ProgressBar():
            good_files = dask.compute(*lengths)
            
        return list(good_files)
    
    def _reprocess_dataset(self, ds):
        """Reprocess the netCDF dataset to conform to CF-standards.

        Parameters
        ----------
        ds: (xarray.DataSet)
            An opened xarray dataset of the netCDF file.

        Returns
        -------
        ds: (xarray.DataSet)
            Reprocessed xarray DataSet
        """
        # Remove the *_qartod_executed variables
        qartod_pattern = re.compile(r"^.+_qartod_executed.+$")
        for v in ds.variables:
            if qartod_pattern.match(v):
                # the shape of the QARTOD executed should compare to the
                # provenance variable
                if ds[v].shape[0] != ds["provenance"].shape[0]:
                    ds = ds.drop_vars(v)

        # Reset the dimensions and coordinates
        ds = ds.swap_dims({"obs": "time"})
        ds = ds.reset_coords()
        keys = ["obs", "id", "provenance", "driver_timestamp",
                "ingestion_timestamp", 'port_timestamp', 'preferred_timestamp']
        for key in keys:
            if key in ds.variables:
                ds = ds.drop_vars(key)
        ds = ds.sortby('time')

        # clear-up some global attributes we will no longer be using
        keys = ['DODS.strlen', 'DODS.dimName',
                'DODS_EXTRA.Unlimited_Dimension', '_NCProperties',
                'feature_Type']
        for key in keys:
            if key in ds.attrs:
                del(ds.attrs[key])

        # Fix the dimension encoding
        if ds.encoding['unlimited_dims']:
            del ds.encoding['unlimited_dims']

        # resetting cdm_data_type from Point to Station and the featureType
        # from point to timeSeries
        ds.attrs['cdm_data_type'] = 'Station'
        ds.attrs['featureType'] = 'timeSeries'

        # update some of the global attributes
        ds.attrs['acknowledgement'] = 'National Science Foundation'
        ds.attrs['comment'] = 'Data collected from the OOI M2M API and reworked for use in locally stored NetCDF files.'

        return ds
    
    def _trim_datasets(self, ds):
        """Internal method trims datasets to avoid overlapping time dimensions."""
        # --------------------------------
        # First, get the deployments
        deployments = self.get_deployments(self.REFDES)

        # --------------------------------
        # Second, get the deployment times
        deployments = deployments.sort_values(by="deploymentNumber")
        deployments = deployments.set_index(keys="deploymentNumber")
        # Shift the start times by (-1) 
        deployEnd = deployments["deployStart"].shift(-1)
        # Find where the deployEnd times are earlier than the deployStart times
        mask = deployments["deployEnd"] > deployEnd
        # Wherever the deployEnd times occur after the shifted deployStart times, replace those deployEnd times
        deployments["deployEnd"][mask] = deployEnd[mask]
        deployments["deployEnd"] = deployments["deployEnd"].apply(lambda x: pd.to_datetime(x))

        # ---------------------------------
        # With the deployments info, can write a preprocess function to filter 
        # the data based on the deployment number
        depNum = np.unique(ds["deployment"])
        deployInfo = deployments.loc[depNum]
        deployStart = deployInfo["deployStart"].values[0]
        deployEnd = deployInfo["deployEnd"].values[0]

        # Select the dataset data which falls within the specified time range
        ds = ds.sel(time=slice(deployStart, deployEnd))

        return ds
    
    def _preprocess(self, ds):
        """Internal preprocess routine to feed into xarray.open_mfdataset"""
        # --------------------------------
        ds = self._reprocess_dataset(ds)
        ds = self._trim_datasets(ds)
        return ds

    def load_netCDF_datasets(self, catalog, goldCopy=False):
        """Open the netCDF files directly from the THREDDS opendap server.

        Parameters
        ----------
        catalog: (list)
            A list of the netCDF datasets from the THREDDS catalog to open
            into an xarray DataSet. Note: Data will NOT be loaded into 
            memory.

        Returns
        -------
        ds: (xarray.DataSet)
            A xarray DataSet of the concatenated and reprocessed netCDF
            datasets into a single xarray DataSet object.

        """
        # Outline how to process the code
        # -------------------------------
        # First, need to rename the filenames with the proper heading
        if goldCopy is True:
            dodsC = self.URLS["goldCopy_dodsC"]
        else:
            dodsC = self.URLS["dodsC"]
        netCDF_files = [re.sub("catalog.html\?dataset=", dodsC, file) for file in catalog]
        
        # -------------------------------
        # Second, 
        netCDF_files = [f+"#fillmismatch" for f in netCDF_files]
        
        # -------------------------------
        # Third, check and remove any files which are malformed
        # and remove the bad ones
        # netCDF_files = self._check_files(netCDF_files)

        # Load the datasets into a concatenated xarray DataSet
        with ProgressBar():
            print("\n"+f"Loading netCDF_files for {self.REFDES}:")
            ds = xr.open_mfdataset(netCDF_files, preprocess=self._preprocess, parallel=True)

        # Add in the English name of the dataset
        refdes = "-".join(ds.attrs["id"].split("-")[:4])
        vocab = self.get_vocab(refdes)
        ds.attrs["Location_name"] = " ".join((vocab["tocL1"].iloc[0],
                                              vocab["tocL2"].iloc[0],
                                              vocab["tocL3"].iloc[0]))

        return ds


