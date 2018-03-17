# Licensed under a 3-clause BSD style license - see LICENSE.rst
def setup(enabled, filename='zchecker.log'):
    import sys
    import logging
    from astropy.time import Time

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

            logfile = logging.FileHandler(filename)
            logfile.setLevel(logging.INFO)
            logfile.setFormatter(formatter)
            logger.addHandler(logfile)

    logger.info('#' * 70)
    logger.info(Time.now().iso + 'Z')
    logger.info('Command line: ' + ' '.join(sys.argv[1:]))
    if enabled:
        for handler in logger.handlers:
            if hasattr(handler, 'baseFilename'):
                logger.info('Logging to ' + handler.baseFilename)

    return logger

class ProgressBar:
    """Progress bar widget for logging.

    Parameters
    ----------
    n : int
      Total number of steps.
    logger : logging.Logger
      The `Logger` object to which to report progress.

    Examples
    --------
    with ProgressBar(1000, logger) as bar:
      for i in range(1000):
        bar.update()

    """

    def __init__(self, n, logger):
        self.n = n
        self.logger = logger

    def __enter__(self):
        self.i = 0
        self.last_tenths = 0
        self.logger.info('-' * 10)
        return self

    def __exit__(self, *args):
        print()
        self.logger.info('#' * 10)

    def update(self):
        self.i += 1
        tenths = int(self.i / self.n * 10)
        if tenths != self.last_tenths:
            self.last_tenths = tenths
            self.logger.info('#' * tenths + '-' * (10 - tenths))
