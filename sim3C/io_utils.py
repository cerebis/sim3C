"""
meta-sweeper - for performing parametric sweeps of simulated
metagenomic sequencing experiments.
Copyright (C) 2016 "Matthew Z DeMaere"

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import bz2
import gzip
import json

import yaml

# default buffer for incremental read/write
DEF_BUFFER = 16384


def open_output(fname, mode, compress=None, gzlevel=6):
    """
    Open a text stream for reading or writing. Compression can be enabled
    with either 'bzip2' or 'gzip'. Additional option for gzip compression
    level. Compressed filenames are only appended with suffix if not included.

    :param fname: file name of output
    :param mode: read or write
    :param compress: gzip, bzip2
    :param gzlevel: gzip level (default 6)
    :return:
    """
    if compress == 'bzip2':
        if not fname.endswith('.bz2'):
            fname += '.bz2'
        # bz2 missing method to be wrapped by BufferedWriter. Just directly
        # supply a buffer size
        return bz2.open(fname, mode)
    elif compress == 'gzip':
        if not fname.endswith('.gz'):
            fname += '.gz'
        return gzip.open(fname, mode, compresslevel=gzlevel)
    else:
        return open(fname, mode)


def multicopy_tostream(fname, *ostreams, **kwargs):
    """
    Copy an input file to multiple output streams.
    :param fname: input file name
    :param ostreams: output streams
    :param kwargs: optional parameters: write_mode (default 'w'), compress [gzip, bzip2] default: None
    :return:
    """
    bufsize = DEF_BUFFER if 'bufsize' not in kwargs else kwargs['bufsize']

    with open(fname, 'r') as in_h:
        done = False
        while not done:
            buf = in_h.read(bufsize)
            if not buf:
                done = True
            for oi in ostreams:
                oi.write(buf)


def write_to_stream(stream, data, fmt='plain'):
    """
    Write an object out to a stream, possibly using a serialization format
    different to default string representation.

    :param stream: open stream to twrite
    :param data: object instance
    :param fmt: plain, json or yaml
    """
    if fmt == 'yaml':
        yaml.dump(data, stream, default_flow_style=False)
    elif fmt == 'json':
        json.dump(data, stream, indent=1)
    elif fmt == 'plain':
        stream.write('{0}\n'.format(data))


def read_from_stream(stream, fmt='yaml'):
    """
    Load an object instance from a serialized format. How, in terms of classes
    the object is represented will depend on the serialized information. For
    generic serialized formats, this is more than likely to involve dictionaries
    for classes with properties.

    :param stream: open stream to read
    :param fmt: yaml or json
    :return: loaded object
    """
    if fmt == 'yaml':
        return yaml.load(stream)
    elif fmt == 'json':
        return json_load_byteified(stream)


"""
Code below taken from Stack Exchange question.
http://stackoverflow.com/questions/956867/how-to-get-string-objects-instead-of-unicode-ones-from-json-in-python

JSON loading with UTF-8 encoding.

The following functions returns JSON results where Unicode strings are converted to UTF-8. Potentially an
unnecessary step as Python will handle referencing these strings transparently, but I wish to keep exchanged
data tables in a single encoding.

Attribution: Mirec Miskuf
"""


def json_loads_byteified(json_text):
    return _byteify(
        json.loads(json_text, object_hook=_byteify),
        ignore_dicts=True
    )


def json_load_byteified(file_handle):
    return _byteify(
        json.load(file_handle, object_hook=_byteify),
        ignore_dicts=True
    )


def _byteify(data, ignore_dicts=False):
    # if this is a unicode string, return its string representation
    if isinstance(data, str):
        return data.encode('utf-8')
    # if this is a list of values, return list of byteified values
    if isinstance(data, list):
        return [_byteify(item, ignore_dicts=True) for item in data]
    # if this is a dictionary, return dictionary of byteified keys and values
    # but only if we haven't already byteified it
    if isinstance(data, dict) and not ignore_dicts:
        return {_byteify(key, ignore_dicts=True): _byteify(value, ignore_dicts=True) for key, value in data.items()}
    # if it's anything else, return it in its original form
    return data
