.. image:: https://readthedocs.org/projects/trapga/badge/?version=latest
    :target: https://trapga.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

GATAI
=============

Overview
-------------------

GATAI  is a Python library built upon setga, designed for extracting genes that play a significant role in development. It utilizes transcriptomic data of genes, spanning multiple developmental stages and their respective gene ages (for more information on how to get the gene ages, see GeneEra_.

.. _GeneEra: https://github.com/josuebarrera/GenEra

The project features an algorithm designed to identify genes contributing to the observed TAI pattern in developmental gene expression data. GATAI aims to identify a subset of genes that, if removed from the dataset, would significantly reduce the presence of the pattern. By employing a multi-objective genetic algorithm, it maximizes the removal of the TAI pattern while minimizing the number of removed genes.

To determine the significance of the TAI pattern during the optimization, GATAI uses the variance sampling introduced in myTAIs flat-line-test_ by Hajk-Georg Drost

.. _flat-line-test: https://drostlab.github.io/myTAI/reference/FlatLineTest.html

The algorithm utilizes the DEAP (Distributed Evolutionary Algorithms in Python) library, which provides a flexible framework for implementing genetic algorithms. It offers various selection methods, mutation operators, and genetic operators to evolve populations of candidate solutions.

Additionally, to enhance its search capability and avoid being trapped in local optima, the algorithm employs an island model. This approach involves maintaining multiple subpopulations, or "islands," that evolve independently. Periodic migration of individuals between islands allows for sharing of genetic information and prevents the algorithm from converging prematurely to suboptimal solutions. This utilization of the island model enhances the algorithm's ability to explore the solution space and discover more globally optimal solutions.


https://github.com/lavakin/ga_for_hourglass/assets/76702901/5435d04b-151b-46da-b179-d48ca9a7e5ce

Features
-------------------

- **Gene Extraction:** Automatically identifies and extracts genes significant to development from transcriptomic data.
- **Built on setga:** Leveraging the capabilities of SetMiG, a library for extracting a minimal subset which optimizes a given function.

Installation
-------------------

.. code-block:: bash

  pip3 install gatai


Citation
-------------------

Please cite myTAI as well

> Drost et al. __myTAI: evolutionary transcriptomics with R__. _Bioinformatics_ 2018, 34 (9), 1589-1590. [doi:10.1093](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btx835/4772684)





Usage
-------------------
The main function of this tool is to identify a small set of genes driving the TAI pattern, for this having a tsv file, where the first column named "Phylostratum" are the gene ages, second column named "GeneID" are the gene ids and the following columns are the expressions for the respective developmental stages. Those columns need to be sorted time-wise and the replicates collapsed.

To identify the genes driving the TAI pattern run:

.. code-block:: bash

   gatai run_minimizer input_data output_folder

In the output folder a text file with identified genes will be stored

If the run statistics should be stored, run

.. code-block:: bash

   gatai run_minimizer input_data output_folder --save_stats

This will save the summary of the run with the elapsed time, number of generations, number of extracted genes etc. The pickled logbook, best solutions for every generation and the final population will be stored as well.

In case you have sampled variances precomputed and saves in a file with the values separated by newline, you can run

.. code-block:: bash

    gatai run_minimizer input_data output_folder --variances variances_file

In case your dataset is single cell, you can run 

.. code-block:: bash

    gatai run_minimizer input_data output_folder --single_cell

to run the version working eith the expression matrix as a sparse matrix. However, due to TAI not being tested for single cell we do not recommend to draw any conclusions from the identified genes.


Contributing
-------------------

Contributions to this project are welcome. If you have any ideas for improvements, new features, or bug fixes, please submit a pull request. For major changes, please open an issue to discuss the proposed modifications.


License
-------------------

This project is licensed under the MIT License. Feel free to use and modify the code according to the terms of this license.
