=======================
Scipion OpenMM plugin
=======================

**Documentation under development, sorry for the inconvenience**

In order to use this plug-in, you need to have Scipion3 installed
(https://scipion-em.github.io/docs/docs/scipion-modes/how-to-install.html).

To install the plugin,  you have to follow the following steps:

1. **Clone this repository:**

.. code-block::

    git clone https://github.com/scipion-chem/scipion-chem-openmm.git


2. **Switch to the desired branch** (master or devel):

Scipion-chem-openmm is constantly under development and including new features.
If you want a relatively older an more stable version, use master branch (default).
If you want the latest changes and developments, user devel branch.

.. code-block::

            cd scipion-chem-openmm
            git checkout devel

3. **Install**:

.. code-block::

    scipion3 installp -p path_to_scipion-chem-openmm --devel -j <numberOfProcessors>
    
OR
    
**Install the plugin in user mode** (not available yet)

.. code-block::

    scipion3 installp -p path_to_scipion-chem-openmm
