# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""config
=========

The configuration is a JSON-formatted file.  See ``_config_example``.

"""
import os
import json
import sbsearch.config

__all__ = ['Config']

_config_example = '''
{
  "database": "/path/to/zchecker.db",
  "log": "/path/to/zchecker.log",
  "user": "IRSA account user name",
  "password": "IRSA account password",
  "cutout path": "/path/to/cutout/directory",
  "cutout size": "5arcmin",
  "stack path": "/path/to/stack/directory"
}
'''


class Config(sbsearch.config.Config):
    """ZChecker configuration.

    Controls database location, log file location, IRSA user
    credentials, object cutout locations, and cutout download size.
    Parameters are stored as object keys: `Config['user']`,
    `Config['log']`, etc., except the IRSA credentials, which are
    retrieved via `Config.auth`.

    Parameters
    ----------
    **kwargs
        Override parameters with these values.

    **kwargs
      Additional or updated configuration parameters and values.

    """

    DEFAULT_FILES = ['./zchecker.config',
                     os.path.expanduser('~/.config/zchecker.config')]

    def __init__(self, **kwargs):
        config = json.loads(_config_example)
        config.update(**kwargs)
        super().__init__(**config)

    @property
    def auth(self):
        return {
            "user": self['user'],
            "password": self['password']
        }
