from distutils.core import setup

setup(
    name='gatai',
    py_modules=['gatai'],
    packages = ['gatai'],
    license='MIT',   
    description = 'Library designed for extracting genes that play a significant role in development. It utilizes transcriptomic data of genes, spanning multiple developmental stages and their respective gene ages',  # Give a short description about your library
    long_description = 'Library designed for extracting genes that play a significant role in development. It utilizes transcriptomic data of genes, spanning multiple developmental stages and their respective gene ages', 
    project_urls = {"Documentation" :"https://trapga.readthedocs.io/en/latest/genindex.html"},
    author = 'Nikola Kalábová',              
    author_email = 'nikola@kalabova.eu',     
<<<<<<< HEAD
    url = 'https://github.com/lavakin/gatai', 
    version = "{{VERSION_PLACEHOLDER}}",      
    #version = 1.1,
=======
    url = 'https://github.com/lavakin/tai_chi', 
    #version = "{{VERSION_PLACEHOLDER}}",      
    version = 1.1,
>>>>>>> parent of 0682a12 (release)
    keywords = ['Genetic algorithms', 'minimal subset', 'multi-objective', "optimization"],   
    install_requires=['numpy', 'scipy', 'pandas', 'argparse', 'scikit-learn', 'tqdm',"setga"],
    entry_points={
        'console_scripts': [
            'gatai = gatai.gatai:cli'
        ]
        },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',  
    ] 
)
