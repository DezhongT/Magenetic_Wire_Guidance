## [Real-time Model-Based Control of MSCR]

***
Numerical simulation and model-based control of a Magnetic Soft Continuum Robot(MSCR) navigating in confined lumens. Uses [Discrete Elastic Rod (DER)](http://www.cs.columbia.edu/cg/pdfs/143-rods.pdf) framework and incorporates hard magnetic theory and contact. By altering the confined lumen's geometry, this repo can compute the optimal actuation magnetic field to maintain the tool tip along the centerline of the lumen so that the contact can be avoided. Simulation examples can be seen below in Figure 1.


<p align="center">
<img src="images/knot_tying.png" alt>
<br>
<em> Figure 1.  </em>
</p>

***


## How to Use

### Dependencies
Install the following C++ dependencies:
- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
  - Eigen is used for various linear algebra operations.
  - IMC is built with Eigen version 3.4.0 which can be downloaded [here](https://gitlab.com/libeigen/eigen/-/releases/3.4.0). After downloading the source code, install through cmake as follows.
    ```bash
    cd eigen-3.4.0 && mkdir build && cd build
    cmake ..
    sudo make install
    ```

- [OpenGL / GLUT](https://www.opengl.org/)
  - OpenGL / GLUT is used for rendering the knot through a simple graphic.
  - Simply install through apt package manager:
      ```bash
    sudo apt-get install libglu1-mesa-dev freeglut3-dev mesa-common-dev
    ```
- Lapack (*usually preinstalled on your computer*)

***
### Compiling
After completing all the necessary above steps, clone the source repository of IMC and then build the project through cmake.
```bash
mkdir build && cd build
cmake ..
make -j4
```
Two complied executors are generated. ```optimizeDER``` is used to compute the optimal actuation magnetic field. ```simDER``` is used to generate the visualization of the MSCR's navigation controlled by the computed actuation magnetic field.

***

### Setting Parameters

All simulation parameters for ```optimizeDER``` are in the parameter file ```option_optimize.txt```; the parameters for ```simDER``` are in the parameter file ```option_sim.txt```.

Specifiable parameters are as follows (we use SI units):
- ```rodRadius``` - Cross-sectional radius of the rod.
- ```RodLength``` - Contour length of the rod.
- ```youngM``` - Young's modulus.
- ```numVertices``` - Number of nodes on the rod.
- ```Poisson``` - Poisson ratio.
- ```deltaTime``` - Time step size.
- ```totalTime``` - Total simulation time
- ```tol``` and ```stol``` - Small numbers used in solving the linear system. A fraction of a percent, e.g. 1.0e-3, is often a good choice.
- ```maxIter``` - Maximum number of iterations allowed before the solver quits.
- ```caseOpt``` - We pre-defined a series of (e.g., 1 to 5) cases with different confined lumen geometry. Here, by inputting different integers from 1 to 5 to visualize different examples.
- ```density``` - Mass per unit volume.
- ```gVector``` - 3x1 vector specifying acceleration due to gravity.
- ```baVector``` - 3x1 vector specifying external magnetic acutation. Set all zeros here since this value will be computed.
- ```brVector``` - 3x1 magnetized direction vector of the MSCR's tip
- ```scaleRender``` - Scale value for rendering the simulation
- ```deltaBaStep``` - the time difference to update the actuation magnetic field
- ```stopTime``` -
- ```timeWait``` -
- ```speed``` - the pushing speed of MSCR
- ```muZero``` -
- ```saveData (0 or 1)``` - Flag indicating whether pull forces and rod end positions should be reocrded.
- ```dBar``` -
- ```stiffness``` - contact stiffness between the MSCR and the lumen's wall
- ```mu``` -
- ```epsilonV``` - hyperparameters for contact (friction-relevant)
- ```viscosity``` - viscosity of damping effect
- ```render (0 or 1) ```- Flag indicating whether OpenGL visualization should be rendered.

***
### Running the Simulation
Once parameters are set to your liking, the simulation can be ran from the terminal by running the provided script:
```bash
./optimizeDER option_optimize.txt
```
to generate the optimal external actuation magnetic field.
Then, run:
```bash
./simDER option_sim.txt
```
to visualize the MSCR's navigation.

***

### Citation
If our work has helped your research, please cite the following paper.
```
@article{tong2025real,
  title={Real-time simulation enabled navigation control of magnetic soft continuum robots in confined lumens},
  author={Tong, Dezhong and Hao, Zhuonan and Li, Jiyu and Sun, Boxi and Liu, Mingchao and Wang, Liu and Huang, Weicheng},
  journal={arXiv preprint arXiv:2503.08864},
  year={2025}
}


```



