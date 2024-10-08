# About The Project
This project code can generate random porous media with Quartet Structure Generation Set (QSGS) method. QSGS was first proposed by Wang and his collaborators[<sup>1</sup>](#refer-1). Then he further enriched this method at 2007[<sup>2</sup>](#refer-2) and 2009[<sup>3</sup>](#refer-3).

## References
<div id="refer-1"></div>
[1] Wang, M., Wang, J., Pan, N., & Chen, S. (2007). Mesoscopic predictions of the effective thermal conductivity for microscale random porous media.  Physical Review E - Statistical, Nonlinear, and Soft Matter Physics, 75(3), 036702. (https://doi.org/10.1103/PhysRevE.75.036702)

<div id="refer-2"></div>
[2] Wang, M., Meng, F., & Pan, N. (2007). Transport properties of functionally graded materials. Journal of Applied Physics, 102(3), 1–7. (https://doi.org/10.1063/1.2767629)

<div id="refer-3"></div>
[3] Wang, M., & Kang, Q. (2009). Electrokinetic transport in microchannels with random roughness. Analytical Chemistry, 81(8), 2953–2961. (https://doi.org/10.1021/ac802569n)

# Run the code
The most of this project is based on the C++ standard library. 

The method **output2png** of Class **Porous2D** and **Porous3D** can generate PNG file directly, which is base on the **GnuPlot** and **gnuplot-iostream**. If you need this function, make sure to install **GnuPlot** and **Boost C++ library** on your system. If you don't need, just comment the related code.

The method **output2tecplot** of Class **Porous2D** and **Porous3D** can generate plt file for [Tecplot](https://www.tecplot.com/) and __GenerateSTL.py__ program.

## How to run
1. Clone the repo
```shell
    git clone https://github.com/Jeff-Hugh/GenPorMed.git
```
2. Compile

    Make sure the boost path in CMakeLists.txt is correct in your system.
```shell
    cmake ./
    make
```
    or just use command
```shell
    g++ Porous2D.cpp Porous3D.cpp run.cpp -lboost_system -lboost_thread -lboost_iostreams -std=c++11 -fopenmp -o run
```

3. Run
```shell
    ./run
```

4. If you want to generate the STL file, you should use the function __output2tecplot__ of Porous2D to output a file like __lattice.dat__. Python program, numpy and pandas Python packages should also be installed. Then run the py file.
```shell
    pip install numpy pandas
    python GenerateSTL.py
```
# License
Distributed under the GPL v3.0 License.

# Examples
2D porous figure
![porosity = 0.5](fig/porous_2D_0.5.png)


3D porous figure
![porosity = 0.5](fig/porous_3D_0.5.png)

# Please cite this repo
```tex
@software{Li_GenPorMed_2024,
author = {Li, Jianhui},
month = aug,
title = {{GenPorMed}},
url = {https://github.com/Jeff-Hugh/GenPorMed},
version = {1.1},
year = {2024}
}
```