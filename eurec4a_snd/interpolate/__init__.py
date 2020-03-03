# Create __version__ attribute from setup.py information
from pkg_resources import get_distribution, DistributionNotFound
import os.path
from ._version import get_versions

__import__('pkg_resources').declare_namespace(__name__)

try:
    __version__ = get_versions()['version']
except:
    __version__ = '--'
del get_versions
