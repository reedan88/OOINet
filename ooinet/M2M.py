import os
import re
import time
import netrc
import requests
import datetime
import logging
import numpy as np
import pandas as pd
import xarray as xr
from queue import Queue
from bs4 import BeautifulSoup
from urllib.request import urlretrieve
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Any

# Import shared utilies (kept exactly as requested)
from ooinet.utils import ntp_seconds_to_datetime, convert_time, unix_epoch_time
from ooinet.Download import setup_download_dir, download_file, DownloadWorker

# configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Initialize credentials
try:
    nrc = netrc.netrc()
    AUTH = nrc.authenticators('ooinet.oceanobservatories.org')
    if AUTH is None:
        raise RuntimeError(
            'No entry found for machine ``ooinet.oceanobservatories.org`` in the .netrc file'
        )
    login, _, password = AUTH  # keep original order but only use login/password
except FileNotFoundError as e:
    raise OSError(e, os.strerror(e.errno), os.path.expanduser('~'))

# Set the URL endpoints for ooinet.oceanobservatories.org
URLS = {
    'data': 'https://ooinet.oceanobservatories.org/api/m2m/12576/sensor/inv',
    'anno': 'https://ooinet.oceanobservatories.org/api/m2m/12580/anno/find',
    'vocab': 'https://ooinet.oceanobservatories.org/api/m2m/12586/vocab/inv',
    'asset': 'https://ooinet.oceanobservatories.org/api/m2m/12587',
    'deploy': 'https://ooinet.oceanobservatories.org/api/m2m/12587/events/deployment/inv',
    'preload': 'https://ooinet.oceanobservatories.org/api/m2m/12575/parameter',
    'cal': 'https://ooinet.oceanobservatories.org/api/m2m/12587/asset/cal',
    'fileServer': "https://opendap.oceanobservatories.org/thredds/fileServer/",
    'dodsC': "https://opendap.oceanobservatories.org/thredds/dodsC/",
    "goldCopy": "https://thredds.dataexplorer.oceanobservatories.org/thredds/catalog/ooigoldcopy/public/",
    "goldCopy_fileServer": "https://thredds.dataexplorer.oceanobservatories.org/thredds/fileServer/",
    "goldCopy_dodsC": "https://thredds.dataexplorer.oceanobservatories.org/thredds/dodsC/"
}

# Set the connections
SESSION = requests.Session()
adapter = requests.adapters.HTTPAdapter(max_retries=0)
SESSION.mount("https://", adapter)
DEFAULT_TIMEOUT = 10  # seconds


def get_api(url: str, params: Optional[Dict[str, Any]] = None, timeout: int = DEFAULT_TIMEOUT) -> Optional[Dict]:
    """
    Perform an authenticated GET request to the OOINet API and return parsed JSON.

    Parameters
    ----------
    url : str
        The URL to request.
    params : dict or None, optional
        Querystring parameters to include in the request.
    timeout : int, optional
        Request timeout in seconds.

    Returns
    -------
    dict or None
        Parsed JSON response if status_code == 200, otherwise None.

    Raises
    ------
    requests.RequestException
        Re-raises network-related exceptions after logging.
    """
    try:
        r = SESSION.get(url, params=params, auth=(login, password), timeout=timeout)
        r.raise_for_status()
        return r.json()
    except requests.HTTPError as e:
        logger.error("HTTP error requesting %s: %s", url, e)
    except requests.RequestException as e:
        logger.error("Network error requesting %s: %s", url, e)
        raise
    except ValueError as e:
        logger.error("JSON decode error for %s: %s", url, e)
    return None


def get_datasets(search_url: str, datasets: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """
    Recursively walk an OOINet endpoint (ending with 'inv') and build a dataframe of datasets.

    Parameters
    ----------
    search_url : str
        Endpoint to start walking; expected to be an OOINet 'inv' endpoint or a subtree.
    datasets : pandas.DataFrame or None, optional
        DataFrame to append results to. If None, an empty DataFrame will be used.

    Returns
    -------
    pandas.DataFrame
        DataFrame with discovered datasets. Columns include array, node, instrument, refdes, url, deployments.
    """
    if datasets is None:
        datasets = pd.DataFrame()
    
    # Check if the method is attached to the url
    flag = ("inv" == search_url.split("/")[-4]) if len(search_url.split("/")) >= 4 else False

    # This means you are at the end-point
    if flag is True:
        # Get the reference designator info
        array, node, instrument = search_url.split("/")[-3:]
        refdes = "-".join((array, node, instrument))

        # Get the available deployments
        deploy_url = "/".join((URLS["deploy"], array, node, instrument))
        deployments = get_api(deploy_url) or []

        # Put the data into a dictionary
        info = {
            "array": array,
            "node": node,
            "instrument": instrument,
            "refdes": refdes,
            "url": search_url,
            "deployments": [deployments],
        }
        
        # Convert it to a dataframe
        info_df = pd.DataFrame(info)
        
        # add the dictionary to the dataframe
        datasets = pd.concat([datasets, info_df], ignore_index=True)

    else:
        endpoints = get_api(search_url)
        
        # endpoints expected to be a list of endpoint names
        while endpoints:
            new_endpoint = endpoints.pop()
            new_search_url = "/".join((search_url, new_endpoint))
            datasets = get_datasets(new_search_url, datasets=datasets)

    # Once recursion is done, return the datasets
    return datasets


def search_datasets(array: Optional[str] = None, node: Optional[str] = None,
                    instrument: Optional[str] = None, English_names: bool = False) -> pd.DataFrame:
    """
    Search OOINet for datasets matching given array/node/instrument.

    Parameters
    ----------
    array : str or None
        OOI array abbreviation, or None to search all.
    node : str or None
        Partial or full node abbreviation to filter.
    instrument : str or None
        Partial or full instrument abbreviation to filter.
    English_names : bool, optional
        If True, include human-readable names (via get_vocab).

    Returns
    -------
    pandas.DataFrame
        DataFrame of matching datasets.
    """
    # Build the request url; if any piece is None, format will include 'None' so strip at first 'None'
    dataset_url = f'{URLS["data"]}/{array}/{node}/{instrument}'
    
    # Truncate at the first "None"
    none_idx = dataset_url.find("None")
    if none_idx != -1:
        dataset_url = dataset_url[:none_idx - 1]

    # Display search to user
    logger.info("Searching %s", dataset_url)
    datasets = get_datasets(dataset_url)

    # filter by node/instrument if provided
    if node is not None and not datasets.empty:
        datasets = datasets[datasets["node"].str.contains(node, na=False)]

    if instrument is not None and not datasets.empty:
        datasets = datasets[datasets["instrument"].str.contains(instrument, na=False)]

    # Check if they want the English names for the associated datasets
    if English_names and not datasets.empty:
        vocab_rows = []
        for refdes in datasets["refdes"]:
            refdes_vocab = get_vocab(refdes)
            if refdes_vocab is None or refdes_vocab.empty:
                vocab_rows.append({
                    "refdes": refdes,
                    "array_name": "No record",
                    "node_name": "No record",
                    "instrument_name": "No record"
                })
            else:
                vocab_rows.append({
                    "refdes": refdes,
                    "array_name": f'{refdes_vocab["tocL1"].iloc[0]} {refdes_vocab["tocL2"].iloc[0]}',
                    "node_name": refdes_vocab["tocL3"].iloc[0],
                    "instrument_name": refdes_vocab["instrument"].iloc[0]
                })
        vocab_df = pd.DataFrame(vocab_rows)
        datasets = datasets.merge(vocab_df, on="refdes", how="left")
        columns = ["array", "array_name", "node", "node_name",
                   "instrument", "instrument_name", "refdes", "url", "deployments"]
        # keep only columns that exist
        columns = [c for c in columns if c in datasets.columns]
        datasets = datasets[columns]

    return datasets


def get_vocab(refdes: str) -> pd.DataFrame:
    """
    Return OOI vocabulary (human readable names) for a reference designator.

    Parameters
    ----------
    refdes : str
        Reference designator in the form '<array>-<node>-<instrument>'.

    Returns
    -------
    pandas.DataFrame
        DataFrame of vocabulary information. Empty DataFrame if none.
    """
    # First, construct the vocab request url
    array, node, instrument = refdes.split("-", 2)
    vocab_url = "/".join((URLS["vocab"], array, node, instrument))

    # Next, get the vocab data
    data = get_api(vocab_url)

    # Put the returned vocab data into a pandas dataframe
    vocab = pd.DataFrame()
    vocab = pd.concat([vocab, pd.DataFrame(data)])

    # Finally, return the results
    return vocab


def reformat_calInfo(calInfo: Dict, deployNum: int, uid: str) -> pd.DataFrame:
    """
    Reformat calibration info JSON into a DataFrame.

    Parameters
    ----------
    calInfo : dict
        JSON structure returned by calibration API.
    deployNum : int
        Deployment number to annotate each row.
    uid : str
        Asset UID.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns: deploymentNumber, uid, calCoef, value, calFile
    """
    rows = []
    for cal in calInfo.get("calibration", []):
        for calData in cal.get("calData", []):
            rows.append({
                "deploymentNumber": int(deployNum),
                "uid": uid,
                "calCoef": calData.get("eventName"),
                "value": calData.get("value"),
                "calFile": calData.get("dataSource")
            })
    return pd.DataFrame(rows)


def get_calibrations_by_refdes(refdes: str, deployments: pd.DataFrame) -> Optional[pd.DataFrame]:
    """
    Get calibration coefficients for all deployments of a reference designator.

    Parameters
    ----------
    refdes : str
        Reference designator in the form '<array>-<node>-<instrument>'.
    deployments : pandas.DataFrame
        Deployments DataFrame as returned by get_deployments().

    Returns
    -------
    pandas.DataFrame or None
        Combined calibrations across deployments or None if none found.
    """
    
    # Create a dask delayed object to gather the calibrations for the
    # instrument from each deployment
    cals_list = []
    for _, row in deployments.iterrows():

        # Get the relevant request data
        deployNum = row.get("deploymentNumber")
        uid = row.get("uid")
        deployStart = row.get("deployStart")
        deployEnd = row.get("deployEnd") or datetime.datetime.now()

        # Construct the request parameters
        params = {
            "refdes": refdes,
            "uid": uid,
            "beginDT": pd.to_datetime(deployStart).strftime("%Y-%m-%dT%H:%M:%S.000Z") if deployStart is not None else None,
            "endDT": pd.to_datetime(deployEnd).strftime("%Y-%m-%dT%H:%M:%S.000Z")
        }

        # Remove None parameters
        params = {k: v for k, v in params.items() if v is not None}

        # Request the data
        calInfo = get_api(URLS["cal"], params) or {}
        if calInfo:
            df = reformat_calInfo(calInfo, deployNum, uid)
            if not df.empty:
                cals_list.append(df)

    # Append individual calibrations into a single dataframe
    if not cals_list:
        calibrations = None
    else:
        calibrations = pd.concat(cals_list, ignore_index=True)
    
    return calibrations


def get_calibrations_by_uid(uid: str, params: Optional[Dict[str, Any]] = None) -> pd.DataFrame:
    """
    Get calibrations for a specific instrument by its UID.

    Return the calibration coefficients for a given instrument
    in OOI for all of its deployments by its UID. The instrument
    may have been deployed at different locations (i.e. reference
    designators).

    Parameters
    ----------
    uid: str
        The unique id (UID) for the instrument for which to request
        vocab information.
    params: (dict), optional
        Optional parameters to pass to the calibration request to limit
        the data returned. Optional keywords are:
            beginDT: earliest datestring ("%Y-%m-%dT%H:%M:%S.000Z") to limit search
            endDT: latest datestring ("%Y-%m-%dT%H:%M:%S.000Z") to limit search

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns uid, calCoef, calDate, value, calFile sorted by calDate and calCoef.
    """
    cal_url = URLS["cal"]
    query = dict(params) if params else {}
    query["uid"] = uid

    # Make the calibration data request
    calInfo = get_api(cal_url, params=query) or {}

    # Put the data into a pandas dataframe, sorted by calibration date and coefficient name
    rows = []
    for c in calInfo.get("calibration", []):
        for cc in c.get("calData", []):
            rows.append({
                "uid": cc.get("assetUid"),
                "calCoef": cc.get("eventName"),
                "calDate": convert_time(cc.get("eventStartTime")),
                "value": cc.get("value"),
                "calFile": cc.get("dataSource")
            })
    if not rows:
        calibrations = pd.DataFrame(columns=["uid", "calCoef", "calDate", "value", "calFile"])
    else:
        calibrations = pd.DataFrame(rows)
        calibrations.sort_values(by=["calDate", "calCoef"], inplace=True)
        calibrations.reset_index(drop=True, inplace=True)

    return calibrations


def get_deployments(refdes: str, deploy_num: str = "-1") -> pd.DataFrame:
    """
    Request deployment information for a reference designator.

    Get the deployment information for an instrument. Defaults to all
    deployments for a given instrument (reference designator) unless one is
    supplied.

    Parameters
    ----------
    refdes: str
        The reference designator for the instrument for which to request
        deployment information.
    deploy_num: str, optional
        Optional to include a specific deployment number. Otherwise
        defaults to -1 which is all deployments.

    Returns
    -------
    pandas.DataFrame
        A table of the deployment information for the given instrument
        (reference designator) with deployment number, deployed water
        depth, latitude, longitude, start of deployment, end of deployment,
        and cruise IDs for the deployment and recovery.
    """
    # First, build the request
    array, node, instrument = refdes.split("-", 2)
    deploy_url = "/".join((URLS["deploy"], array, node, instrument, deploy_num))

    # Next, get the deployments from the deploy url. The API returns a list
    # of dictionary objects with the deployment data.
    deployments = get_api(deploy_url) or []

    # Now, iterate over the deployment list and get the associated data for
    # each individual deployment
    rows = []
    for deployment in deployments:
        # Deployment Number
        deploymentNumber = deployment.get("deploymentNumber")

        # Location info
        location = deployment.get("location", {})
        depth = location["depth"]
        lat = location["latitude"]
        lon = location["longitude"]

        # Sensor info
        sensor = deployment.get("sensor", {})
        uid = sensor.get("uid")
        assetId = sensor.get("assetId")

        # Start and end times of the deployments
        startTime = convert_time(deployment.get("eventStartTime"))
        stopTime = convert_time(deployment.get("eventStopTime"))

        # Cruise IDs of the deployment and recover cruises
        deployCruiseInfo = deployment.get("deployCruiseInfo") if deployCruiseInfo else None
        recoverCruiseInfo = deployment.get("recoverCruiseInfo") if recoverCruiseInfo else None
        
        # Put the data into a dictionary
        rows.append({
            "deploymentNumber": deploymentNumber,
            "uid": uid,
            "assetId": assetId,
            "latitude": lat,
            "longitude": lon,
            "depth": depth,
            "deployStart": startTime,
            "deployEnd": stopTime,
            "deployCruise": deployID,
            "recoverCruise": recoverID
        })

    # Generate the results table
    if not rows:
        results = pd.DataFrame(columns=["deploymentNumber", "uid", "assetId", "latitude", "longitude",
                                        "depth", "deployStart", "deployEnd", "deployCruise", "recoverCruise"])
    else:
        results = pd.DataFrame(rows)
        # Sort the deployments by deployment number and reset the index
        results = results.sort_values(by="deploymentNumber")
        results = results.reset_index(drop=True)

    return results


def get_datastreams(refdes: str) -> pd.DataFrame:
    """
    Retrieve available methods and streams for a reference designator.

    Parameters
    ----------
    refdes : str
        Reference designator '<array>-<node>-<instrument>'.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns refdes, method, stream (each stream is a single row).
    """
    # Build the url
    array, node, instrument = refdes.split("-", 2)
    method_url = "/".join((URLS["data"], array, node, instrument))
    methods = get_api(method_url) or []

    # Build a table linking the reference designators, methods, and data streams
    rows = []
    for method in methods:
        if "bad" in method:
            continue
        stream_url = "/".join((method_url, method))
        streams = get_api(stream_url) or []
        for s in streams:
            rows.append({"refdes": refdes, "method": method, "stream": s})

    # Create the DataFrame table
    if not rows:
        stream_df = pd.DataFrame(columns=["refdes", "method", "stream"])
    else:
        stream_df = pd.DataFrame(rows)

    return stream_df



def parse_metadata(metadata: Dict) -> pd.DataFrame:
    """
    Parse metadata JSON for an instrument into a DataFrame.

    Parameters
    ----------
    metadata : dict
        JSON metadata as returned by the OOINet metadata endpoint.

    Returns
    -------
    pandas.DataFrame
        DataFrame combining 'times' and 'parameters' keyed on 'stream'.
    """
    # Put the two keys into separate dataframes
    metadata_times = pd.DataFrame(metadata.get("times", []))
    metadata_parameters = pd.DataFrame(metadata.get("parameters", []))

    # Merge the two into a single dataframe
    if metadata_times.empty and metadata_parameters.empty:
        results = None
    else:
        
    results = metadata_parameters.merge(metadata_times, left_on="stream",
                                        right_on="stream")
    results.drop_duplicates(inplace=True)

    # Return the results
    return results


def get_metadata(refdes: str) -> pd.DataFrame:
    """Request metadata.

    Get the OOI Metadata for a specific instrument specified by its
    associated reference designator.

    Parameters
    ----------
    refdes: str
        OOINet standardized reference designator in the form of
        <array>-<node>-<instrument>.

    Returns
    -------
    pandas.DataFrame
        A dataframe with the relevant metadata of the given reference
        designator.
    """
    # First, construct the metadata request url
    array, node, instrument = refdes.split("-", 2)
    metadata_request_url = "/".join((URLS["data"], array, node, instrument, "metadata"))

    # Request the metadata
    metadata = get_api(metadata_request_url) or {}

    # Parse the metadata
    parsed = parse_metadata(metadata)
    if parsed is None:
        parsed = pd.DataFrame()

    # Add in the reference designator
    parsed["refdes"] = refdes

    # Return the metadata
    return parsed


def get_parameter_data_levels(metadata):
    """
    Determine processing/data levels for parameter IDs in metadata via the preload API.

    Parameters
    ----------
    metadata: pandas.DataFrame
        The dataframe of metadata returned by get_metadata() which contains
        the metadata for a given reference designator.

    Returns
    -------
    dict
        Mapping of pdId -> data_level (if available)
    """
    pdIds = np.unique(metadata["pdId"])
    pid_dict = {}
    for pid in pdIds:
        # Build the preload url
        preload_url = "/".join((URLS["preload"], pid.strip("PD")))
        # Query the preload data
        preload_data = get_api(preload_url)
        data_level = preload_data.get("data_level")
        # Update the results dictionary
        pid_dict.update({pid: data_level})

    return pid_dict

def filter_science_parameters(metadata: pd.DataFrame) -> pd.DataFrame:
    """
    Filter metadata to only science-level (processed) parameters.

    Parameters
    ----------
    metadata : pandas.DataFrame
        DataFrame returned by get_metadata.

    Returns
    -------
    pandas.DataFrame
        Filtered metadata containing only science-level parameters.
    """
    if metadata is None or metadata.empty:
        return pd.DataFrame()
        
    data_levels = get_parameter_data_levels(metadata)

    def is_science(pdId):
        level = data_levels.get(pdId)
        return (level is not None) and (level > 0)

    mask = metadata["pdId"].apply(is_science)
    return metadata[mask]


def science_parameters(metadata: pd.DataFrame, remove: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Group and return science-level parameters, optionally removing streams by keywords.

    Parameters
    ----------
    metadata : pandas.DataFrame
        The metadata DataFrame from get_metadata().
    remove : list of str or None, optional
        List of keywords used to filter out streams (case-insensitive substring matching).

    Returns
    -------
    pandas.DataFrame
        Grouped DataFrame indexed by refdes, method, stream with unique values aggregated.
    """
    if metadata is None or metadata.empty:
        return pd.DataFrame()

    # Filter for only the metadata parameters associated with science-variables
    science_vars = filter_science_parameters(metadata)
    if remove:
        # case-insensitive check against the 'particleKey' column (if present)
        for keyword in remove:
            keyword_low = keyword.lower()
            if "particleKey" in science_vars.columns:
                science_vars = science_vars[~science_vars["particleKey"].str.lower().str.contains(keyword_low, na=False)]
            else:
                # if particleKey doesn't exist, attempt to filter on stream names
                science_vars = science_vars[~science_vars["stream"].str.lower().str.contains(keyword_low, na=False)]
    
    # Group the variables togther by the data stream
    grouped = (science_vars.groupby(by=["refdes", "method", "stream"])
               .agg(lambda x: pd.unique(x.values.ravel()).tolist())
               .reset_index())
    
    # flatten single element lists
    grouped = grouped.applymap(lambda x: x[0] if isinstance(x, list) and len(x) == 1 else x)
    
    return grouped

    
def get_annotations(refdes, **kwargs):
    """Retrieve data annotations for a given reference designator.

    Parameters
    ----------
    refdes: str
        The reference designator which to query OOINet for the
        associated annotation data.

    Kwargs
    ------
    method: str
        An OOINet method associated with the reference designator.
        Limits annotations for the given reference designator
        to only that method.
    stream: str
        An OOINet stream associated with the reference designator.
        Limits annotations for the given reference designator
        to only that stream.
    beginDT: str
        Limit the data request to only data after this date. Date
        should be formatted using OOINet epoch time (can be
        calculated using _)
    endDT: str
        Limit the data request to only data before this date.
    """
    # Need to build a parameters dictionary to pass to requests
    params = {"refdes": refdes}
    for key, val in kwargs.items():
        if key in ("beginDT", "endDT"):
            params[key] = unix_epoch_time(val)
        else:
            params[key] = val

    # Get the annotations as a json and put into a dataframe
    anno_data = get_api(URLS["anno"], params=params) or []
    anno_ddf = pd.DataFrame(anno_data)

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
    if not anno_df.empty and 'qcFlag' in anno_df.columns:
        anno_df['qcFlag'] = anno_df['qcFlag'].map(codes).astype('category')

    return anno_df


def get_thredds_url(refdes: str, method: str, stream: str, goldCopy: bool = False,
                    timeout_sec: int = 300, **kwargs) -> Optional[str]:
    """
    Return the THREDDS dataset URL for a given refdes/method/stream.

    Parameters
    ----------
    refdes : str
        Reference designator '<array>-<node>-<instrument>'.
    method : str
        Data method (e.g. 'telemetered').
    stream : str
        Stream name.
    goldCopy : bool, optional
        If True, return the goldCopy catalog URL (static).
    timeout_sec : int, optional
        Maximum seconds to wait for async dataset creation before raising.

    Kwargs
    ------
    beginDT: str or datetime-like
        Limit the data request to only data after this date.
    endDT: str or datetime-like
        Limit the data request to only data before this date.
    format: str
        e.g. "application/netcdf" (the default)
    include_provenance: str
        'true' returns a text file with the provenance information
    include_annotations: str
        'true' returns a separate text file with annotations for the data
        within the given date range

    Returns
    -------
    str or None
        THREDDS URL when ready, otherwise None (or raises on network error).
    """
    # Build the data request url
    if goldCopy is True:
        goldCopy_url = URLS["goldCopy"] + "-".join((refdes, method, stream)) + "/catalog.html"
        return goldCopy_url
    else:
        array, node, instrument = refdes.split("-", 2)
        data_request_url = "/".join((URLS["data"], array, node, instrument, method, stream))

    # Ensure proper datetime format for the request
    params = dict(kwargs)
    for key in ("beginDT", "endDT"):
        if key in params:
            params[key] = pd.to_datetime(params[key]).strftime('%Y-%m-%dT%H:%M:%S.%fZ')

    # Ensure require_deployment included
    params.setdefault("require_deployment", True)

    # Request the data
    logger.info("Waiting for request to process for %s %s %s", refdes, method, stream)
    urls = get_api(data_request_url, params=params) or {}
    all_urls = urls.get("allURLs", [])
    
    # find status_url
    status_candidates = [u for u in all_urls if re.search(r'async_results', u)]
    if not status_candidates:
        # If no async_results, try to find thredds url directly
        thredds_candidates = [u for u in all_urls if 'thredds' in u]
        thredds_url = thredds_candidates[0] if thredds_candidates else None
    else:
        status_url = status_candidates[0].rstrip("/") + "/status.txt"
        start = time.time()
        while True:
            try:
                status_resp = SESSION.get(status_url, timeout=DEFAULT_TIMEOUT)
                if status_resp.status_code == requests.codes.ok:
                    break
            except requests.RequestException:
                logger.debug("Waiting for async result at %s", status_url)
            elapsed = time.time() - start
            if elapsed > timeout_sec:
                logger.error("Timed out waiting for THREDDS URL after %s seconds", timeout_sec)
                raise TimeoutError(f"Timed out waiting for THREDDS URL (status {status_url})")
            time.sleep(2)

        # The asynchronous data request is contained in the 'allURLs' key,
        # in which we want to find the url to the thredds server
        for d in all_urls:
            if 'thredds' in d:
                thredds_url = d
            else:
                thredds_url = None

    return thredds_url


def get_thredds_catalog(thredds_url: str) -> List[str]:
    """
    Get the dataset catalog for the requested data stream.

    Parameters
    ----------
    thredds_url: str
        the THREDDS server url for the requested data stream

    Returns
    -------
    list of str
        The THREDDS catalog of datasets for the requested data stream
    """
    # Parse out the dataset_id from the thredds url
    try:
        page = SESSION.get(thredds_url, timeout=DEFAULT_TIMEOUT).text
    except requests.RequestException as e:
        logger.error("Error fetching THREDDS catalog %s: %s", thredds_url, e)
        return []
    
    soup = BeautifulSoup(page, "html.parser")
    pattern = re.compile(r'.*\.nc$')
    catalog = sorted([node.get('href') for node in soup.find_all('a', text=pattern)])
    
    # Return the catalog
    return catalog


def download_netCDF_files(catalog: List[str], goldCopy: bool = False, saveDir: Optional[str] = None, verbose: bool = True):
    """
    Download netCDF files for given netCDF datasets.

    Downloads the netCDF files returned by parse_catalog. If no path is
    specified for the save directory, will download the files to the
    current working directory.

    Parameters
    ----------
    catalog: list of str
        The netCDF catalog which to download
    goldCopy: bool, optional
        If you should use the gold copy thredds server (which is static)
    saveDir: str or None, optional
        The path to the directory where to download the netCDF files. If no path
        is specified will download to your current directory
    verbose: bool, optional
        If True, show status of download.
    """
    # Step 1 - parse the catalog to get the correct web address to download the files
    fileServer = URLS["goldCopy_fileServer"] if goldCopy else URLS["fileServer"]
    netCDF_files = [re.sub(r"catalog.html\?dataset=", fileServer, file) for file in catalog]

    # Step 2 - Specify and make the relevant save directory
    setup_download_dir(saveDir)
    target_dir = Path(saveDir or os.getcwd())

    # Step 3 - Initialize workers
    # 3A. Identify number availabe cores and divide in half
    cores = os.cpu_count() or 1
    workers = max(1, cores // 2)
    # 3B. Setup a Queue
    queue = Queue()
    # 3C. Start workers
    for _ in range(workers):
        w = DownloadWorker(queue)
        w.daemon = True
        w.start()

    # Step 4. Execute the processes and download the files
    if verbose:
        logger.info("----- Downloading files -----")
    for file_url in netCDF_files:
        queue.put((str(target_dir), file_url))
    # Execute
    queue.join()


def clean_catalog(catalog, stream, deployments=None):
    """
    Clean a THREDDS catalog list by removing empty datasets and filtering by stream or deployments.

    Parameters
    ----------
    catalog : list of str
        Raw catalog entries (hrefs).
    stream : str
        Target stream substring to keep as primary datasets.
    deployments : pandas.DataFrame or None, optional
        If provided, filters datasets to only include datasets with deployment numbers appearing in this DataFrame.

    Returns
    -------
    tuple (datasets, ancillary)
    datasets : list of str
        Filtered list of datasets matching the requested stream and deployments.
    ancillary : list of str
        Remaining ancillary dataset entries.
    """
    # Next, check that the netCDF datasets are not empty by getting the timestamps in the
    # datasets and remove datasets with identical t1 == t2 timestamps (empty) 
    valid_datasets = []
        try:
            timestamps = dset.split("_")[-1].replace(".nc", "").split("-")
            t1, t2 = timestamps[0], timestamps[-1]
            if t1 != t2:
                valid_datasets.append(dset)
        except Exception:
            # if pattern unexpected, keep the dataset (safer)
            valid_datasets.append(dset)
            
    # Next, check if you want to filter for certain deployments
    if deployments is not None:
        # Check that the deployments are a list
        dep_list = list(deployments["deploymentNumber"].astype(int).tolist())
        filtered = []
        for dset in valid_datasets:
            match = re.search(r"deployment(\d{4})", dset)
            if not match:
                continue
            depNum = int(match.group(1))
            if depNum in dep_list:
                filtered.append(dset)
        valid_datasets = filtered      
        
    # Finally, determine if the dataset is either for the given instrument
    # or an ancillary instrument which supplies and input variable
    datasets = []
    ancillary = []
    for dset in valid_datasets:
        basename = dset.split("/")[-1]
        if stream in basename:
            datasets.append(dset)
        else:
            ancillary.append(dset)
                       
    return datasets, ancillary