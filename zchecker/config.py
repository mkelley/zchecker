# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""config
=========

The configuration is a JSON-formatted file:

{
  "database": "/path/to/zchecker.db",
  "log": "/path/to/zchecker.log",
  "user": "IRSA account user name",
  "password": "IRSA account password",
  "cutout path": "/path/to/cutout/target/directory"
}

"""
# Configuration file format should match the description in
# scripts/zchecker help.

class Config:
    """ZChecker configuration.

    Controls database location, log file location, IRSA user
    credentials, and object cutout locations.  Parameters are stored
    as object keys: `Config['user']`, `Config['log']`, etc., except
    the IRSA credentials, which are retrieved via `Config.auth`.

    Parameters
    ----------
    filename : string
      The file to load, or `None` to load the default location
      '~/.config/zchecker.config'.

    **kwargs
      Additional or updated configuration parameters and values.

    """

    def __init__(self, filename=None, **kwargs):
        import os
        import json
        
        if filename is None:
            filename = os.path.expanduser('~/.config/zchecker.config')

        with open(filename) as f:
            self.config = json.load(f)

        self.config.update(kwargs)

        self.auth = {
            'user': self.config.pop('user'),
            'password': self.config.pop('password')
        }

    def __getitem__(self, k):
        return self.config[k]
