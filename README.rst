.. image:: https://readthedocs.org/projects/trapga/badge/?version=latest
    :target: https://trapga.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

TraP-GA 
=============

Overview
-------------------

The TraP-GA is a Python library built upon SetMiG, designed for extracting genes that play a significant role in development. It utilizes transcriptomic data of genes, spanning multiple developmental stages and their respective gene ages (for more information on how to get the gene ages, see [GeneEra](https://github.com/josuebarrera/GenEra)).

The project features an algorithm designed to identify genes contributing to the observed hourglass pattern in developmental gene expression data. The hourglass pattern refers to a characteristic developmental pattern observed in various organisms, where the morphological and genetic similarities between individuals are most pronounced at a specific stage of development.

This algorithm aims to identify a subset of genes that, if removed from the dataset, would significantly reduce the presence of the hourglass pattern. By employing a multi-objective genetic algorithm, it maximizes the removal of the hourglass pattern while minimizing the number of removed genes.

The algorithm utilizes the DEAP (Distributed Evolutionary Algorithms in Python) library, which provides a flexible framework for implementing genetic algorithms. It offers various selection methods, mutation operators, and genetic operators to evolve populations of candidate solutions.

Additionally, to enhance its search capability and avoid being trapped in local optima, the algorithm employs an island model. This approach involves maintaining multiple subpopulations, or "islands," that evolve independently. Periodic migration of individuals between islands allows for sharing of genetic information and prevents the algorithm from converging prematurely to suboptimal solutions. This utilization of the island model enhances the algorithm's ability to explore the solution space and discover more globally optimal solutions.

https://github.com/lavakin/ga_for_hourglass/assets/76702901/5435d04b-151b-46da-b179-d48ca9a7e5ce

Features
-------------------

- **Gene Extraction:** Automatically identifies and extracts genes significant to development from transcriptomic data.
- **Built on SetMiG:** Leveraging the capabilities of SetMiG, a library for extracting a minimal subset which optimizes a given function.

Installation
-------------------

.. code-block:: bash

  git clone https://github.com/lavakin/TraP-GA

.. code-block:: bash

  chmod +x DEAP_island.py


Usage
-------------------

.. code-block:: bash

   ./DEAP_island.py data.csv output_folder --save_plot




Contributing
-------------------

Contributions to this project are welcome. If you have any ideas for improvements, new features, or bug fixes, please submit a pull request. For major changes, please open an issue to discuss the proposed modifications.


License
-------------------

This project is licensed under the MIT License. Feel free to use and modify the code according to the terms of this license.
