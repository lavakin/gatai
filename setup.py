from setuptools import setup

setup(
    name='trapga',
    version='0.1',
    py_modules=['trapga'],
    install_requires=['numpy', 'scipy', 'pandas', 'argparse', 'scikit-learn', 'tqdm',"setminga"],
    entry_points={
        'console_scripts': [
            'trapga = trapga:cli'
        ]
    }
)

