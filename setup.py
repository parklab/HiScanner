from setuptools import setup, find_packages

# Read requirements
with open('requirements.txt') as f:
    requirements = [line.strip() for line in f if line.strip() and not line.startswith('#')]

# Read version
exec(open('hiscanner/__version__.py').read())

setup(
    name="hiscanner",
    version=__version__,
    description="Single-cell Copy Number Variation Analysis Pipeline",
    author="Yifan Zhao",
    author_email="yifnzhao@gmail.com",
    packages=find_packages(),
    install_requires=requirements,
    include_package_data=True,
    package_data={
        'hiscanner': [
            'resources/*',
            'resources/default_config.yaml',
            'resources/cluster.yaml',
            'resources/cluster_config.yaml',
            'resources/Snakefile',
        ],
    },
    entry_points={
        'console_scripts': [
            'hiscanner=hiscanner.cli:main',
        ],
    },
    python_requires='>=3.8',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
)