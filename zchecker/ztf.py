# Licensed under a 3-clause BSD style license - see LICENSE.rst
def query(params, auth):
    import requests
    from astropy.io import ascii
    
    # https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci?WHERE=obsdate>'2017-12-05 12:00'+AND+obsdate<'2017-12-06 12:00'+AND+ccdid=6+AND+qid=4&ct=html

    print('Querying IRSA...')
    r = requests.get(
        'https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci',
        auth=(auth['user'], auth['password']),
        params=params)
    print(r.url)

    print('Done.  {} lines returned.'.format(r.text.count('\n')))

    try:
        return ascii.read(r.text)
    except ascii.InconsistentTableError as e:
        print(r.text)
        raise e
