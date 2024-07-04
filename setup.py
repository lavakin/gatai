from distutils.core import setup

setup(
    name='tai_chi',
    py_modules=['tai_chi'],
    packages = ['tai_chi'],
    license='MIT',   
    description = 'Library designed for extracting genes that play a significant role in development. It utilizes transcriptomic data of genes, spanning multiple developmental stages and their respective gene ages',  # Give a short description about your library
    long_description = 'Library designed for extracting genes that play a significant role in development. It utilizes transcriptomic data of genes, spanning multiple developmental stages and their respective gene ages', 
    project_urls = {"Documentation" :"https://trapga.readthedocs.io/en/latest/genindex.html"},
    author = 'Nikola Kalábová',              
    author_email = 'nikola@kalabova.eu',     
    url = 'https://github.com/lavakin/tai_chi', 
    version = "{{VERSION_PLACEHOLDER}}",      
    #version = 1.1,
    keywords = ['Genetic algorithms', 'minimal subset', 'multi-objective', "optimization"],   
    install_requires=['numpy', 'scipy', 'pandas', 'argparse', 'scikit-learn', 'tqdm',"setga"],
    entry_points={
        'console_scripts': [
            'tai_chi = tai_chi.tai_chi:cli'
        ]
        },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',  
    ] 
)
