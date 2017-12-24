# Licensed under a 3-clause BSD style license - see LICENSE.rst
def setup(enabled):
    import sys
    import logging
    from datetime import datetime

    logger = logging.Logger('ZChecker')
    logger.setLevel(logging.DEBUG)

    # This test allows logging to work when it is run multiple times
    # from ipython
    if len(logger.handlers) == 0:
        formatter = logging.Formatter('%(levelname)s: %(message)s')

        if enabled:
            console = logging.StreamHandler(sys.stdout)
            console.setLevel(logging.DEBUG)
            console.setFormatter(formatter)
            logger.addHandler(console)

            logfile = logging.FileHandler('zchecker.log')
            logfile.setLevel(logging.INFO)
            logfile.setFormatter(formatter)
            logger.addHandler(logfile)

    logger.info('#' * 70)
    logger.info(datetime.now().isoformat())
    logger.info('Command line: ' + ' '.join(sys.argv[1:]))

    return logger
