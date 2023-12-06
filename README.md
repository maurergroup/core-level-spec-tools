# core-level-spec-tools

## Installation 
Usage, installation, and dependency management is automated using poetry. To setup the virtual environment with dependencies, execute 

``` sh
poetry install
```

after cloning the repository. Then to build the project, 

``` sh
poetry build
```

Finally, to make scripts available system-wide, from the root directory of the project, 

``` sh
pip install .
```

This is currently not available on PyPi, however this may change in the future

## Usage 
### NEXAFS 
Currently, the only tools that have been interfaced with poetry are `nexafs` and `molpdos`, however other scripts in the repository are still able to be used. To view usage options for `nexafs`, 

``` sh
nexafs --help
```

`molpdos` currently does not take any command line arguments, so any variables need to be changed within the script itself, however this is likely to change in the future. More extensive documentation will also be written at a later date on usage of these tools. 
