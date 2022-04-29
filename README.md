# Mitochodria Retrograde signalling ODE model

> Codes, data and notebooks for the manuscript titled "Mathematical Modeling and Analysis of Mitochondrial Retrograde Signaling Dynamics"

Forked from Steven's results: <https://github.com/ntumitolab/MitoRetroDynamics>

![Graphical abstract](https://user-images.githubusercontent.com/29009898/130342513-081f4592-3cc6-4468-ba3b-868416f3be6b.png)

## Requirements

1. Jupyter lab
2. Julia (version 1.7+)
3. IJulia

### Pip3 installation

jill (https://github.com/abelsiqueira/jill) is a light installer of Julia on Linux. To install jill, pip3 (https://pip.pypa.io/en/stable/) is required.

**pip3 installation**

```
sudo apt-get install python3-pip
```

### Install Scipy

```
pip3 install scipy
```

### Install matplotlib

```
pip3 install matplotlib
```

### JupyterLab installation

The simulations are executed and displayed via JupyterLab. It is an IDE for literal programming. In this project, Julia kernel is used in Jupyterlab, and extra installation is needed to use Julia with Jupyterlab.

```
pip3 install jupyterlab
```

See https://jupyterlab.readthedocs.io/en/stable/getting_started/installation.html for further instruction.

### Installing Julia

Julia programming language runtime can be installed in the official website (https://julialang.org/downloads/)


### IJulia

IJulia (https://julialang.github.io/IJulia.jl/stable/manual/installation/) is a Julia package for kernel installation for Jupyter.

The installation requires Julia Package manager (Pkg, https://docs.julialang.org/en/v1/stdlib/Pkg/), which is a built-in package in Julia. To install IJulia, one needs to open Julia REPL (Read–eval–print loop) using the following command line.

```
julia
```
![Julia REPL](https://user-images.githubusercontent.com/29009898/130343508-7d8e5e18-7ca8-46f8-b3de-3e4910c42ff3.png)

The Julia REPL begins in the terminal after submitting `julia` command. The next step is to install IJulia with package manager (Pkg). Copy-paste the following code block to Julia REPL

```julia
using Pkg
Pkg.add("IJulia")
Pkg.build("IJulia")
```
![IJulia installation with Pkg](https://user-images.githubusercontent.com/29009898/130343609-a997c935-a209-4364-8c5c-e2e62a3e42b0.png)

### Enabling a multi-threaded kernel

Some notebooks use multithread algorithms to speed up. The multithread process can be set by IJulia in Julia REPL.

```julia
using IJulia
IJulia.installkernel("Julia 8 Threads", env=Dict(
    "JULIA_NUM_THREADS" => "8",
))
```

where `JULIA_NUM_THREADS` represents number of threads that can be used in the Julia process, and this can not be changed after a Julia process is initiated. Though the number of threads is set as `8`, it can be changed and usually depends on the hardware setup.

## Set up of simulation environments

Open the project folder with terminal, and then in the Julia REPL prompt

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Run notebooks

IN the Julia REPL prompt

```julia
using IJulia
IJulia.jupyterlab(dir=pwd())
```

Select the kernel with the name `Julia 8 Threads` to initiate Julia with multiple threads. After successfully initiation, one can see `Julia 8 Threads 1.6` located at the top right corner.

![Jupyter screen shot](https://user-images.githubusercontent.com/29009898/130344079-b98a76f1-13c5-4b8c-a195-dd6f45bd5527.png)

