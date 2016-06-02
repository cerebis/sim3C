import json
import yaml


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
    if isinstance(data, unicode):
        return data.encode('utf-8')
    # if this is a list of values, return list of byteified values
    if isinstance(data, list):
        return [_byteify(item, ignore_dicts=True) for item in data]
    # if this is a dictionary, return dictionary of byteified keys and values
    # but only if we haven't already byteified it
    if isinstance(data, dict) and not ignore_dicts:
        return {
            _byteify(key, ignore_dicts=True): _byteify(value, ignore_dicts=True)
            for key, value in data.iteritems()
            }
    # if it's anything else, return it in its original form
    return data
