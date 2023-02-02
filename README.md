# Evaluating and boosting a single-cell proteomics study


This repository is dedicated to the evaluation of a single-cell proteomics study and boosting of the existing approach with an uncomplicated solution, using good code practices and following the FAIR principles (making the data and code **F**indable, **A**ccessible, **I**nteroperable and **R**eusable).

The existing code can be found in a Jupyter notebook in the following repository: https://github.com/theislab/singlecell_proteomics. The initial observations were that the code is not reusable or  scalable, it is not organised, there is plenty of repeated code, no clear and concise comments, no tests and the data is not fitting the FAIR principles. The code is lacking most of the best coding practices, which is what was the goal of this code boosting.
</br>

## Installation and run

### Poetry

The programming language of the existing repository is Python. Dependency management was done with Poetry. Poetry is a tool for dependency management, packaging and publishing in Python. It allows us to declare the libraries a project depends on and it will manage them for us (including installation and updating). Poetry offers a lockfile to ensure repeatable installs (poetry.lock), and can build a project for distribution.

Once the repository is cloned, next step would be to access the folder where the project is placed and then running the _poetry_insall_ command. After tthe installation, run the _poetry_shell_ command which spawns a shell, or equivalently, activates the poetry virtual environment. This environment will have all the packages and dependencies from the poetry.lock predefined file. 

### Tox

In addition to poetry, tox was used for generic virtual environment management. It is a test command line tool which can be used for code checking and testing. In this implementation, tox is used to assure code best practices such as PEP8 code style guide (flake8), formatting (black), sorting of libraries (isort) and static type checking (mypy). It can be easily extended to support unit-testing through testing frameworks like pytest, however due to time constrains, unit-tests were not implemented. 

- Use  _tox -e fix_styling_ command for automatic code formatting
- Use _tox -e linters_ command for automatic code checks (eg. if all the datatypes match)

### Data

The data provided (access over https://drive.google.com/file/d/1_JMv0-9q40cVLLTft9hzHh2gDzP7DEyr/view?usp=share_link) should be placed in the bioexperiment/data/experiment_data path. 

This concludes the setup necessary for the project to work. 

- To run the project, run the command _python bioexperiment/experiments/proteomics_experiment.py_.

## Implemented best practices

The repository has been organised and documented in a way to follow the Coding best practices. Here are some of the principles that were followed thorough the project development:

- Code reusability and scalability - Firstly, the existing code was evaluated, and based on the analysis it has been concluded that on multiple places, the code was being repeated. To solve this issue, the code had been transformed with the use of objects, classes and functions through Object-oriented and Functional Programming and the principles of inheritance and polymorphism. There are many types of plots in the code that can be seen more than one time, however there are minor differences between the 'almost-duplicates', and this was included in the boosted code. 
In order to run the code with different data, one needs to redo the loader_configuration. The load methods can be reused since the file formats are the same and the variables need to be rewired inside the code. Since the notebook / analysis is structured like this, the code cannot be 'plug and play'. 
- Interpretability - the code has been structured in such a way to be easily readable and interpretable.
- Portability - after following the steps explained in the 'Installation' section, the code is ready to be used. This code works in different environments (processors, operating systems, different versions of libraries etc.). All the requirements can be found in the poetry.lock file. 
- Versioning - GitHub was used to keep track of code modifications, dependencies, and also to enable an easy and transparent way of tracking the changes of the project code. 
-  Formatting - the code was formatted according to standards using black (including the indentation ). Moreover, the functions are not lengthy, and the lines of code are not long. 
- Naming conventions - the mentioned classes, functions as well as variables had been named in a self-explanatory way and follow a consistent theme throughout the code following the naming best practices. 
- Clear and concise comments - the updated code functions and classes are described with DocStrings in a concise and easily readable manner. Also, there are additional comments added in parts of the code where thought necessary. 
- Typing - for easier understanding of the inputs and outputs of the methods (with tox)
- Testing - the code quality can be tested with tox (see commands above). 

## Code organisation

The code is structured in the following way (directories):

- data/
    - experimental_data/
    - loader.py
    - loader_definitions.py
    - preprocess.py
- experiments/
    - base_experiment/py
    - proteomics_experiment.py
- visualisations/
    - config.py
    - plots.py
- (results/) 

And that is mainly the order that should be followed. Firstly, the data (paths specified in the loader_definitions) are loaded in the loader.py folder, and there are functions that read all different types of data that belong to this project. 
After the data is loaded, it is preprocessed for further analysis and visualisations. 
Which brings us to the second folder group dedicated to the experiments, and the main experiment file (proteomics_experiment) is split in three analytical parts with its 'sub-functions' (just like the Jupyter notebook is split initially), that is: 

- True Single-Cell Proteomics - Cell Cycle Analysis
- True Single-Cell Proteomics - Comparison with Transcriptomics Data
- Supplementary analyses (for the sake of time, this group was omitted; the idea is completely the same as in the first two groups)

The third folder is formed for the sake of repository organisation, so it is separated from the other two folders. This folder holds all the plotting definitions. 

Lastly, the folder mentioned in brackets above is where the plots get saved, for further analysis.

## Next steps

Here are a few next steps, or additional good coding practices that I would implement besides the practices mentioned above:

- Unit tests - for the sake of isolating the code to test it thoroughly and determine if every part of it works as intended. It is a very important step in the development process, because it can help detect issues in the code very early on the code development process, which otherwise would not be detected (or would be very hard to detect) later on. In addition, with unit tests there is more control over the code, and errors can be avoided. 
- Containerisation - a faster and more secure way of project deployment. In the traditional way (Container-free), code is being developed in a specific coding environment which often can result in errors and bugs when transferred to a new location. For eg. Docker Containers offer fast deployment, fast migration, ease of maintaining and moving the application, better security and fewer software dependencies. 
