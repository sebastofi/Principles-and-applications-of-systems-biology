# Project in Principles-and-applications-of-systems-biology

Goal: recreate figures 3-6 produduced by the authors of the article 'Improving the phenotype predictions of a yeast genome-scale metabolic model by incorporating enzymatic constraints'. (https://www.embopress.org/doi/full/10.15252/msb.20167411)

Method: Use GECKO, a model developed by the authors, and apply the content learned in the exercise sessions

## Installation of GECKO (https://github.com/SysBioChalmers/GECKO)

__Github of the GECKO mentioned in the paper:__
- Gecko v1.0 (https://github.com/SysBioChalmers/GECKO/releases/tag/v1.0)
  - README file of the installed source code indicates to install **SBML** (https://sbml.org/software/libsbml/libsbml-docs/installation/)
    - run installSBML.m in matlab to install SBML
  - **cobratoolbox** (https://github.com/opencobra/cobratoolbox)
    - make sure to add the IBM path with subfolders, initCobraToolbox.m will need it to run (IBM Cplex is one of the solver that is required)
- not sure how to run it. I followed the steps of cloning the repo on their GECKO page, probably running v3? 
   
__Installing steps of the required packages for GECKO:__
  1) Gurobi (this is another solver, not sure if it is redundant with IBM Cplex, or if step 2) needs necessarily Gurobi)
     - Register @ Gurobi
     - install and run https://www.gurobi.com/downloads/gurobi-software/, restart computer
     - a new tab opens (https://portal.gurobi.com/iam/licenses/request/), generate a license key (one should be connected to the EPFL VPN, and ask the academic license) then, you request the license and it generates you a specific key with the format: grggetkey
     and run in your terminal (tutorial link: https://support.gurobi.com/hc/en-us/articles/4534161999889), save license at default location.
     - For running it, you should run in the command lines of matlab the following path:
     
     '>>' cd /Library/gurobi1003/macos_universal2/matlab
     '>>' gurobi_setup
       tutorial-link (https://support.gurobi.com/hc/en-us/articles/4533938303505-How-do-I-install-Gurobi-for-Matlab-)

  2) RAVEN
     - clone repository (git clone --depth=1 https://github.com/SysBioChalmers/RAVEN.git)
     - change to folder, run 'checkInstallation'
  3) Docker

## Working with GECKO

check out protocols .. https://github.com/SysBioChalmers/GECKO/tree/main/tutorials/full_ecModel
