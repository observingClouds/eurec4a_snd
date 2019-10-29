import setuptools
import versioneer

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="eurec4a_snd",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="eurec4a folks",
    author_email="",
    description="Common EUREC4A sounding standard",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/observingClouds/eurec4a_snd",
    packages=setuptools.find_packages(),
    setup_requires=['setuptools-git-version'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=['pillow>=6.0.0', 'matplotlib>=3.1.0', 'basemap>=1.2.0',
                      'numpy>=1.14.0', 'netCDF4>=1.4.0'],
    entry_points={'console_scripts':
                    ['sounding_converter=eurec4a_snd.L1_bufr:main',
                     'sounding_visualize=eurec4a_snd.make_quicklooks_rs41:main']},
    package_data={"eurec4a_snd": ["examples/data/*", "config/meta_information_template.ini"]}
)
