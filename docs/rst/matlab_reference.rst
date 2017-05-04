Function reference
==================

Throughout this document the sum :math:`\sum_\tau` implies the periodic sum :math:`\sum_{p
\in \mathbb Z^3}` with :math:`\tau = (p_1L_1, p_2L_2, p_3L_3)`.

SE
--

Ewald summation for the 3-periodic electrostatic potential

.. math::

   u(x) = \sum_{\tau} \sum_{n=1}^N \frac{q_n}{|x - x_n + \tau|} .

.. automodule:: SE
   :members:

SE2P
----

Ewald summation for the 2-periodic electrostatic potential.

.. automodule:: SE2P
   :members:

SE1P
----

Ewald summation for the 1-periodic electrostatic potential.

.. automodule:: SE1P
   :members:      

SE2P_Stokes
-----------
Ewald summation for the 2-periodic stokeslet potential.

.. automodule:: SE2P_Stokes
   :members:

SE_fast_gridding
----------------

Routines for fast gaussian gridding (FGG), which is central to Spectral Ewald.

.. automodule:: SE_fast_gridding
   :members:

SE_Rotlet
----------------


Ewald summation for the 3-periodic rotlet potential

.. math::

   u(x) = \sum_{\tau} \sum_{n=1}^N R(x - x_n + \tau) t_n

where

.. math::
   R_{ij}(r) = \epsilon_{ijk} \frac{r_k}{|r|^3} .


.. automodule:: SE_Rotlet
   :members:

SE_Stokes
----------------

Ewald summation for the 3-periodic stokeslet potential

.. math::

   u(x) = \sum_{\tau} \sum_{n=1}^N S(x - x_n + \tau) f_n,

where

.. math::

   S(r) = \frac{I}{|r|} + \frac{r \otimes r}{|r|^3}.


.. automodule:: SE_Stokes
   :members:

SE_Stresslet
----------------


Ewald summation for the 3-periodic stresslet potential

.. math::

   u(x) = \sum_{\tau} \sum_{m=1}^N q_m T(x - x_m + \tau) \hat n_m,

where

.. math::

   T_{ijk}(r) = -6\frac{r_i r_j r_k}{|r|^5} .


.. automodule:: SE_Stresslet
   :members:





