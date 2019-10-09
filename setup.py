import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="eurec4a_snd",
    version_format='{tag}.dev{commitcount}+{gitsha}',
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
    entry_points={'console_scripts': ['sounding_converter:eurec4a_snd.L1-rs41.py:main']
)
