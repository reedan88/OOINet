import re
import collections
import numpy as np
import xarray as xr


class Instrument():
    """Super class which contains shared methods amongst instrument types."""

    def update_dataset(self, ds, depth):
        """Update a data set with global and variable metadata attributes.

        Updates a data set with global and variable level metadata attributes
        and sets appropriate dimensions and coordinate axes.

        Parameters
        ----------
        ds: (xarray.dataset)
            Dataset to update
        depth: (float or int)
            Instrument deployment depth

        Returns
        -------
        ds: (xarray.dataset)
            The updated data set
        """
        # add a default station identifier as a coordinate variable to the
        # data set
        ds.coords['station'] = 0
        ds = ds.expand_dims('station', axis=None)
        ds['station'].attrs = dict({
            'cf_role': 'timeseries_id',
            'long_name': 'Station Identifier',
            'comment': ds.attrs['subsite'].upper()
        })

        # determine if the latitude and longitude are set as global attribute
        # or a variable, and parse accordingly
        if 'lat' in ds.variables:
            lat = ds.lat.values[0][0]
            lon = ds.lon.values[0][0]
            ds.drop(['lat', 'lon'])
        else:
            lat = ds.attrs['lat']
            lon = ds.attrs['lon']
            del(ds.attrs['lat'])
            del(ds.attrs['lon'])

        # add the geospatial coordinates using the station identifier from
        # above as the dimension
        geo_coords = xr.Dataset({
            'lat': ('station', [lat]),
            'lon': ('station', [lon]),
            'z': ('station', [depth])
        }, coords={'station': [0]})

        geo_attrs = dict({
            'station': {
                'cf_role': 'timeseries_id',
                'long_name': 'Station Identifier',
                'comment': ds.attrs['subsite'].upper()
            },
            'lon': {
                'long_name': 'Longitude',
                'standard_name': 'longitude',
                'units': 'degrees_east',
                'axis': 'X',
                'comment': 'Deployment location'
            },
            'lat': {
                'long_name': 'Latitude',
                'standard_name': 'latitude',
                'units': 'degrees_north',
                'axis': 'Y',
                'comment': 'Deployment location'
            },
            'z': {
                'long_name': 'Depth',
                'standard_name': 'depth',
                'units': 'm',
                'comment': 'Instrument deployment depth',
                'positive': 'down',
                'axis': 'Z'
            }
        })
        for v in geo_coords.variables:
            geo_coords[v].attrs = geo_attrs[v]

        # merge the geospatial coordinates into the data set
        ds = ds.merge(geo_coords)

        # update coordinate attributes for all variables
        for v in ds.variables:
            if v not in ['time', 'lat', 'lon', 'z', 'station']:
                ds[v].attrs['coordinates'] = 'time lon lat z'

        # update some variable attributes to get somewhat closer to IOOS
        # compliance, more importantly convert QC variables to bytes and set
        # the attributes to define the flag masks and meanings.
        ds['deployment'].attrs['long_name'] = 'Deployment Number'
        qc_pattern = re.compile(r'^.+_qc_.+$')
        executed_pattern = re.compile(r'^.+_qc_executed$')
        results_pattern = re.compile(r'^.+_qc_results$')
        flag_masks = np.array([1, 2, 4, 8, 16, 32, 64, 128], dtype=np.uint8)
        for v in ds.variables:
            if qc_pattern.match(v):     # update QC variables
                ds[v] = (('station', 'time'), [[np.uint8(x) for
                                                x in ds[v].values[0]]])
                ds[v].attrs['long_name'] = re.sub('Qc', 'QC',
                                                  re.sub('_', ' ', v.title()))

                if executed_pattern.match(v):   # *_qc_executed variables
                    ds[v].attrs['flag_masks'] = flag_masks
                    ds[v].attrs['flag_meanings'] = (
                        'global_range_test local_range_test spike_test ' +
                        'poly_trend_test stuck_value_test gradient_test ' +
                        'propogate_flags')
                    ds[v].attrs['comment'] = 'Automated QC tests executed for the associated named variable.'

                    ancillary = re.sub('_qc_executed', '', v)
                    ds[v].attrs['ancillary_variables'] = ancillary
                    if 'standard_name' in ds[ancillary].attrs:
                        ds[v].attrs['standard_name'] = ds[ancillary].attrs[
                            'standard_name'] + ' qc_tests_executed'

                if results_pattern.match(v):    # *_qc_results variables
                    ds[v].attrs['flag_masks'] = flag_masks
                    ds[v].attrs['flag_meanings'] = (
                        'global_range_test_passed local_range_test_passed ' +
                        'spike_test_passed poly_trend_test_passed ' +
                        'stuck_value_test_passed gradient_test_passed ' +
                        'all_tests_passed')
                    ds[v].attrs['comment'] = ('QC result flags are set to true (1) if the test passed. Otherwise, if ' +
                                              'the test failed or was not executed, the flag is set to false (0).')

                    ancillary = re.sub('_qc_results', '', v)
                    ds[v].attrs['ancillary_variables'] = ancillary
                    if 'standard_name' in ds[ancillary].attrs:
                        ds[v].attrs['standard_name'] = ds[ancillary].attrs[
                            'standard_name'] + ' qc_tests_results'

        # convert the time values from a datetime64[ns] object to a floating
        # point number with the time in seconds
        ds['time'] = self.dt64_epoch(ds.time)
        ds['time'].attrs = dict({
            'long_name': 'Time',
            'standard_name': 'time',
            'units': 'seconds since 1970-01-01 00:00:00 0:00',
            'axis': 'T',
            'calendar': 'gregorian'
        })

        # return the data set for further work
        return ds

    def dt64_epoch(self, dt64):
        """Convert datetime64 to epoch time stamp.

        Convert a panda or xarray date/time value represented as a datetime64
        object (nanoseconds since 1970) to a float, representing an epoch time
        stamp (seconds since 1970-01-01).

        Parameters
        ----------
        dt64: (pd.datetime64)
            Pandas or xarray datetime64 object

        Returns
        -------
            Epoch time as seconds since 1970-01-01
        """
        epts = dt64.values.astype(float) / 10.0 ** 9
        return epts

    def dict_update(self, source, overrides):
        """Update nested dictionary or similar mapping.

        From https://stackoverflow.com/a/30655448. Modifies ``source`` in
        place. Replaces original dict_update used by poceans-core, also pulled
        from the same thread.
        """
        for key, value in overrides.items():
            if isinstance(value, collections.Mapping) and value:
                returned = self.dict_update(source.get(key, {}), value)
                source[key] = returned
            else:
                source[key] = overrides[key]
        return source
