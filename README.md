# Agent-based model implementation for CWD dynamics 

## Overview
The following repo contains code for implementing grid-based ABM modeling for Chronic Wasting Disease (CWD). We utilize `Agents.jl`, a Julia-based framework for efficient ABM simulation. We implement [Modeling Routes of Chronic Wasting Disease
Transmission: Environmental Prion Persistence Promotes
Deer Population Decline and Extinction](https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0019896&type=printable), which provides an aspatial model for long-term CWD transmission in deer.

## Setup
### Installation
#### macOS/Unix
- Install Julia:
    - Run `curl -fsSL https://install.julialang.org | sh` in a terminal
    - Restart terminal, type in `julia` to make sure it does something

#### Windows:
- Install Julia:
    - Run `winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore`
    - Restart Powershell, type in `julia` to make sure it does something

### Environment
VSCode is the environment of choice. Add the Julia extension on VSCode for linting and a language server. Navigate to the `intro.jl` file and press F5 (or go to the debug tab and press play), and the necessary packages should be loaded and you should see some output