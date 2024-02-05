__version__ = '0.7'

__copyright__ = """Copyright (C) 2019 Matthew Z DeMaere
This is free software.  You may redistribute copies of it under the terms of
the GNU General Public License v3 <hhttps://www.gnu.org/licenses/gpl-3.0.en.html>.
There is NO WARRANTY, to the extent permitted by law.
"""


def version_stamp(full):
    """
    Create a string indicating the version and possibly extended details such as copyright
    :param full: when True add extended details (multi-line)
    :return: a version stamp string
    """
    if full:
        return 'Sim3C version {}\n{}'.format(__version__, __copyright__)
    else:
        return '{}'.format(__version__)


def date_stamp():
    """
    :return: Return a datetime stamp in the format YYYY-mm-dd hh:mm:ss.f
    """
    from datetime import datetime
    _now = datetime.now()
    return _now.strftime('%Y-%m-%d %H:%M:%S.%f')


def runtime_info():
    """
    :return: Return runtime info (version and date stamps) as a dict
    """
    return {'sim3C_version': f'{__version__}', 'run_timestamp': date_stamp()}
