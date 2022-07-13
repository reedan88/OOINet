import os
import datetime
import pandas as pd
from urllib.request import urlopen, urlretrieve

def ntp_seconds_to_datetime(ntp_seconds):
    """Convert OOINet timestamps to unix-convertable timestamps."""
    # Specify some constant needed for timestamp conversions
    ntp_epoch = datetime.datetime(1900, 1, 1)
    unix_epoch = datetime.datetime(1970, 1, 1)
    ntp_delta = (unix_epoch - ntp_epoch).total_seconds()

    return datetime.datetime.utcfromtimestamp(ntp_seconds - ntp_delta)


def convert_time(ms):
    """Calculate UTC timestamp from OOI milliseconds"""
    if ms is None:
        return None
    else:
        return datetime.datetime.utcfromtimestamp(ms/1000)


def unix_epoch_time(date_time):
    """Convert a datetime to unix epoch microseconds."""
    # Convert the date time to a string
    date_time = int(pd.to_datetime(date_time).strftime("%s"))*1000
    return date_time


def setup_download_dir(save_dir=None):
    """Setup the directory to download files to. If no path/directory specified,
    defaults to current working directory.
    """
    if saveDir is not None:
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
    else:
        save_dir = os.getcwd()
        
        
def download_file(download_dir, link):
    """Download a given link/file to the given download directory"""
    link = os.path.basename(link)
    # Check that the download directory
    setup_download_dir(download_dir)
    # Download the file to the download directory
    download_path = "/".join((download_dir, link))
    urlretrieve(link, download_path)