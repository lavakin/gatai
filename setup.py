from distutils.core import setup

setup(
    name='trapga',
    py_modules=['trapga'],
    packages = ['trapga'],
    license='MIT',   
    description = 'Library designed for extracting genes that play a significant role in development. It utilizes transcriptomic data of genes, spanning multiple developmental stages and their respective gene ages',  # Give a short description about your library
    long_description = 'Library designed for extracting genes that play a significant role in development. It utilizes transcriptomic data of genes, spanning multiple developmental stages and their respective gene ages', 
    project_urls = {"Documentation" :"https://trapga.readthedocs.io/en/latest/genindex.html"},
    author = 'Nikola Kalábová',              
    author_email = 'nikola@kalabova.eu',     
    url = 'https://github.com/lavakin/trapga', 
    version = "{{VERSION_PLACEHOLDER}}",      
    keywords = ['Genetic algorithms', 'minimal subset', 'multi-objective', "optimization"],   
    install_requires=['numpy', 'scipy', 'pandas', 'argparse', 'scikit-learn', 'tqdm',"setga"],
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
