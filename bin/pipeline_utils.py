import yaml


def write_to_stream(stream, data):
    yaml.dump(data, stream, default_flow_style=False)


def read_from_stream(stream):
    return yaml.load(stream)


def write_data(file_name, data):
    """
    Simple YAML output format for scoring or other pipeline associated data.

    :param file_name: the file to write output
    :param data: the data object
    """
    with open(file_name, 'w') as stream_out:
        write_to_stream(stream_out, data)


def read_data(file_name):
    """
    Read simple YAML format file of pipeline associated data.

    :param file_name: the file from which to read input
    :return: loaded data object
    """
    with open(file_name, 'r') as stream_in:
        return read_from_stream(stream_in)
