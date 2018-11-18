# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
from contextlib import contextmanager
from subprocess import CalledProcessError
from .exceptions import DownloadError


def query(params, auth, logger=None):
    import requests
    from astropy.io import ascii

    # https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci?WHERE=obsdate>'2017-12-05 12:00'+AND+obsdate<'2017-12-06 12:00'+AND+ccdid=6+AND+qid=4&ct=html

    print('Querying IRSA...')
    r = requests.get(
        'https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci',
        auth=(auth['user'], auth['password']),
        params=params)

    if logger:
        logger.debug(r.url)

    print('Done.  {} lines returned.'.format(r.text.count('\n')))

    try:
        return ascii.read(r.text)
    except ascii.InconsistentTableError as e:
        print(r.text)
        raise e


class IRSA:
    """Context manager for IRSA connections.

    Parameters
    ----------
    path : string
      Directory in which to save cookies.txt file.

    auth : dictionary
      IRSA 'user' and 'password'.

    """

    def __init__(self, path, auth):
        self.path = path
        self.auth = auth

    def __enter__(self):
        url = (
            "https://irsa.ipac.caltech.edu/account/signon/login.do?"
            "josso_cmd=login&josso_username={}&josso_password={}"
        ).format(self.auth['user'], self.auth['password'])

        self._wget(url, '/dev/null', save_cookies=True)

        return self

    def __exit__(self, *args):
        self._wget("https://irsa.ipac.caltech.edu/account/signon/logout.do",
                   '/dev/null', save_cookies=True)

    def download(self, url, fn, clean_failed=True, logger=None):
        """Download from IRSA.


        Parameters
        ----------
        url : string
            The full URL of the file.

        fn : string
            The local file name of the downloaded data.

        clean_failed : bool, optional
            Delete files left over from failed downloads?

        logger : Logger, optional
            Log errors to ``logger``.


        Raises
        ------
        DownloadError

        """

        try:
            self._wget(url, fn)
        except CalledProcessError as e:
            # if logger:
            #    logger.error('Error downloading {} from {}: {}'.format(
            #        fn, url, str(e)))
            if os.path.exists(fn) and clean_failed:
                os.unlink(fn)
            raise DownloadError from e
        except KeyboardInterrupt as e:
            os.unlink(fn)
            raise e

    def _wget(self, url, fn, save_cookies=False):
        import sys
        from subprocess import check_call

        args = ['wget']
        if save_cookies:
            args.append('--save-cookies={}/cookies.txt'.format(self.path))
        else:
            args.append('--load-cookies={}/cookies.txt'.format(self.path))
        args.extend(['-O', fn, '-o', '/dev/null', url])

        check_call(args, stdout=sys.stderr)
