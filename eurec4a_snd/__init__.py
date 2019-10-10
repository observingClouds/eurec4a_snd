# Create __version__ attribute from setup.py information
from pkg_resources import get_distribution, DistributionNotFound
import os.path

__import__('pkg_resources').declare_namespace(__name__)

try:
    _dist = get_distribution('eurec4a_snd')
    # Normalize case for Windows systems
    dist_loc = os.path.normcase(_dist.location)
    here = os.path.normcase(__file__)
    if not here.startswith(os.path.join(dist_loc, 'eurec4a_snd')):
        # not installed, but there is another version that *is*
        raise DistributionNotFound
except DistributionNotFound:
    __version__ = 'Please install this project with setup.py'
else:
    __version__ = _dist.version
