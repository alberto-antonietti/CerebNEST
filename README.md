# CerebNEST

Tested with:

*branch nest_3_8*: Rocky 8 server, Python 3.13.3 and NEST Release 3.8


### Installation instructions

0. Install NEST with `mamba` (see full instructions on NEST documentation). E.g.:

```
mamba create --name mamba_nest_3_8 -c conda-forge nest-simulator jupyterlab seaborn
[...]
mamba activate mamba_nest_3_8
```   

2. Export an Environment Variable containing the installation directory of NEST. E.g.:
```
export NEST_BIN=$HOME/miniforge3/envs/mamba_nest_3_8/bin/nest-config
```

2. Clone this GitHub Repository in a directory outside NEST source and build directories. E.g.:
```
cd $HOME
git clone https://github.com/alberto-antonietti/CerebNEST/
```

3. Move to CerebNEST directory and create a new folder where you will build the extension module
```
mkdir $HOME/cerebellar_module_build
cd $HOME/cerebelar_module_build

```
4. Run the following CMake command (Tested with CMake 3.2.2)
```
cmake -Dwith-nest="${NEST_BIN}" $HOME/CerebNEST
```

The resulting output should be something similar to:
> [...]
>-------------------------------------------------------
>You can now build and install 'cerebellar_module' using
> 
>  make
> 
>  make install
>
>The library file libcerebellar_module.so will be installed to
> 
>  /home/aantonie/miniforge3/envs/mamba_nest_3_8/lib/nest
> 
>The module can be loaded into NEST using
> 
>  (cerebellar_module) Install       (in SLI)
> 
>  nest.Install(cerebellar_module)   (in PyNEST)
> 
>
>-- Configuring done (0.8s)
> 
>-- Generating done (0.0s)
> 
>-- Build files have been written to: /home/aantonie/cerebellar_module_build
> 

5. Make and install the module
```
make -j8
make install
```

7. Every time you need the module, you can install it in this way:
```
python

```

```
import nest
nest.Install("cerebellar_module")

```

8. You can now use all the nodes and synapses contained in this NEST Module.
