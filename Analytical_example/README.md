# (Adaptive) Sparse Grid - analytic example 

To illustrate how (adaptive) sparse grids can be used, we provide an example in plain python and a Jupyter notebook, which approximates the following function:


![\begin{align*}
f_d(\vec{x}) = \max(0,1-e^{\frac{1}{2}-(\prod_{i=1}^d(x_i+\frac{1}{5}))^{\frac{1}{d}}}),\;\,\text{with}\,\vec{x}=\{x_1,\dots,x_d\}\in[0,1]^d
\end{align*}
](https://render.githubusercontent.com/render/math?math=%5CLarge+%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0Af_d%28%5Cvec%7Bx%7D%29+%3D+%5Cmax%280%2C1-e%5E%7B%5Cfrac%7B1%7D%7B2%7D-%28%5Cprod_%7Bi%3D1%7D%5Ed%28x_i%2B%5Cfrac%7B1%7D%7B5%7D%29%29%5E%7B%5Cfrac%7B1%7D%7Bd%7D%7D%7D%29%2C%5C%3B%5C%2C%5Ctext%7Bwith%7D%5C%2C%5Cvec%7Bx%7D%3D%5C%7Bx_1%2C%5Cdots%2Cx_d%5C%7D%5Cin%5B0%2C1%5D%5Ed%0A%5Cend%7Balign%2A%7D%0A)

as described in [Section 2.4](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3979412). We generate sparse grid as well as adaptive sparse grid approximations to the function and measure their accuracy as follows: We randomly generate 10,000 uniformly distributed test points, and compute the maximum error and the L2-error.

The code is organized as follows:

    1. Computation of errors for (non-adaptive) sparse grids from level 6 to 16 (including the respective output files).
    
    2. Computation of errors for adaptive sparse grids, starting from level 3 and allowing for 
        refinement steps from level 6 to level 16 (including the respective output files).
    3. Visualization of grid points from a (non-adaptive) sparse grid of level 6 compared to 
        the grid points from an adaptive sparse grid that is refined up to level 16.
        
    4. Visualization of the performance in steps 1. and 2.

Alternative test functions to measure the performance of function approximators were discribed for instance by [Genz (1984)](https://dl.acm.org/doi/10.5555/2837.2842), and available 
[here](https://www.sfu.ca/~ssurjano/integration.html).


### Prerequisites / Installation

To run the code for either the analytical example or the international real business cycle (IRBC) model, follow the instructions below.

**TASMANIAN 7.7**

Tasmanian is included in the Python Pip index: https://pypi.org/project/Tasmanian/

```shell
$ python3 -m pip install scikit-build packaging numpy --user (required dependencies)
$ python3 -m pip install Tasmanian --user                    (user installation)
$ python3 -m pip install Tasmanian                           (virtual env installation)
```

## Usage
Launch the python code with
```shell
    $ cd <PATH to the repository>/SparseGridsForDynamicEcon/SparseGrids_in_econ_handbook/Analytical_example
    $ python 2D_example.py
```
or alternatively (Jupyter notebook)

```shell
    $ cd <PATH to the repository>/SparseGridsForDynamicEcon/SparseGrids_in_econ_handbook/Analytical_example
    $ jupyter-notebook 2D_example.ipynb
```
