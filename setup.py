from setuptools import setup

setup(
    name='trapga',
    version='0.3',
    py_modules=['trapga'],
    packages = ['trapga'],
    license='MIT',   
    description = 'Library designed for extracting genes that play a significant role in development. It utilizes transcriptomic data of genes, spanning multiple developmental stages and their respective gene ages',   # Give a short description about your library
    author = 'Nikola Kalábová',              
    author_email = 'nikola@kalabova.eu',     
    url = 'https://github.com/lavakin/trapga',  
    download_url = 'https://github.com/lavakin/trapga/archive/refs/tags/v0.3.tar.gz',    
    keywords = ['Genetic algorithms', 'minimal subset', 'multi-objective', "optimization"],   
    install_requires=['numpy', 'scipy', 'pandas', 'argparse', 'scikit-learn', 'tqdm',"setminga"],
    entry_points={
        'console_scripts': [
            'trapga = trapga.trapga:cli'
        ]
        },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',  
    ] 
)

