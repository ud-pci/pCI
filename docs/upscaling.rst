Upscaling basis sets
====================

When running larger and larger CI calculations, you may run into problems due to the lack of available computational resources. In this section, we will discuss a method of generating effective configuration spaces by building on top of a subspace of important configurations, i.e. upscaling the basis set. Applications of this method includes expanding a basis set to include configurations with higher principle quantum numbers, higher partial waves, extra configurations and more.

Theory
------

We begin with with an initial configuration space :math:`H_0` from which we were able to obtain the corresponding low-lying eigenvalues :math:`E_0` and eigenvectors :math:`X_0`. Our goal is then to obtain the low-lying eigenvalues :math:`E` and eigenvectors :math:`X` of a larger configuration space :math:`H=H_0+H_1`, where :math:`H_1` is a configuration space generated from an expansion of the initial configuration space :math:`H_0`. Here, the CI calculation for :math:`H` is not possible due to prohibitive computational resources, but we can very well approximate :math:`H` by stripping :math:`H_0` of unimportant configurations, thus reducing the size of the overall configuration space. The goal then is to construct an effective configuration space 

.. math::

    \begin{aligned}
    H^*=H_0^*+H_1 \\
    \approx H_0+H_1
    \end{aligned},

where :math:`H_0^*` is a subspace of :math:`H_0` with important configurations. It is important to check that the energies do not shift significantly when stripping :math:`H_0` of unimportant configurations by running the CI calculation for the :math:`H_0^*` configuration space, e.g. make sure :math:`\Delta E_0=E_0^*-E_0 \approx 0`. 

We can then run the CI procedure to obtain the energies :math:`E^*` and eigenvectors :math:`X^*`` in the effective configuration space :math:`H^*`. Assuming :math:`\Delta E_0 \approx 0` as mentioned above, we can approximate the energies of the full larger configuration space :math:`H` as :math:`E\approx E^*-\Delta E_0`.

Method
------

We can now apply the upscaling method described above using the pCI code package. This method can be repeated to expand the basis set as many times as necessary. 

Run 1
*****

We start from a base run corresponding to :math:`H_0`, where we were successful in obtaining the desired energies and eigenvectors via ``conf``. From this run, we can use the ``con_cut`` program to obtain ``CON_CUT.RES``, which is a ``CONF.INP`` file corresponding to a smaller configuration space corresponding to :math:`H_0^*`, with the important configurations. Note that when running ``con_cut``, you will be given the option to choose a cutoff value based on the weights of configurations. This subspace should have differences of energies within 5 :math:`\text{cm}^{-1}` from the full configuration space, e.g. :math:`E_0^*-E_0 \leq 5` :math:`\text{cm}^{-1}`. Note this 5 :math:`\text{cm}^{-1}` is an arbitrary value to ensure an optimal configuration space, where nearly all important configurations are still included. 

Run 2
*****

Next, we will need to generate the configuration subspace :math:`H_1` corresponding to the desired expansion. This is done by running ``concmp``, which takes in two ``CONF.INP`` files named ``C_A.INP`` and ``C_B.INP``. Here, ``C_A.INP`` should correspond to the configuration space of the initial full run :math:`H_0` and ``C_B.INP`` should correspond to the configuration space of the full desired expansion :math:`H`. Inputting arguments ``0 0 1``, ``concmp`` compares the two configuration lists and outputs a new configuration list ``C_B_new.INP``, which includes all configurations that are in ``C_B.INP``, but not in ``C_A.INP``. We can then use the ``merge_ci`` program to merge ``CON_CUT.RES``, corresponding to :math:`H_0^*`, and ``C_B_new.INP``, corresponding to :math:`H_1`, to obtain ``C_M.INP``, which corresponds to :math:`H^*=H^*_0+H_1`. Note that that program ``merge_ci`` requires two specific input files, so one must rename ``CON_CUT.RES`` to ``C_A.INP`` and ``C_B_new.INP`` to ``CONF.INP``. Note that this specific choice of renaming ``CON_CUT.RES`` to ``C_A.INP`` allows the initial matrix for diagonalization to be constructed in the :math:`H_0^*` subspace. We can then rename the merged ``C_M.INP`` to ``CONF.INP`` and run ``conf`` to obtain the desired energies and eigenvectors. 

The final energies from the expanded basis set can then be calculated by subtracting the energies resulting from ``conf`` using ``CON_CUT.RES``. 

Example
-------

Here is an example of upscaling a basis set from 17g to 18g:

1. Run ``conf`` for full 17g
2. Run ``con_cut`` on 17g to obtain 17g_cut configurations in ``CON_CUT.RES`` 

    * find a ``log_cutoff`` threshold value where the energies obtained from the resulting configuration list is near-identical to those of the full 17g run

3. Run ``concmp`` to obtain (18g-17g) configurations, i.e. configurations only in 18g

    * 17g in ``C_A.INP``, 18g in ``C_B.INP``
    * run ``concmp`` with input ``0 0 1`` to obtain output ``C_B_new.INP``

4. Run ``merge_ci`` to combine 17g_cut with (18g-17g) configurations

    * rename ``C_B_new.INP`` to ``CONF.INP``
    * rename ``CON_CUT.RES`` from 17g_cut to ``C_A.INP``
    
5. Run ``conf`` to obtain energies for 17g_cut + (18g-17g)
6. Obtain energy difference of 18g-17g by subtracting out 17g_cut energies