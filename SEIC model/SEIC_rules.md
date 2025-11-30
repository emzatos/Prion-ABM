# Model Overview:
The agent-based model simulates CWD transmission in a population of deer living on a 2D grid. Disease can spread through direct transmission from infected deer agents (S, E, I, C) to susceptible agents and indirect transmission from environmental sources (disease reservoir V) to susceptible agents.
- Grid Setup: Agents are arranged on a 2D grid with a Chebyshev metric (i.e. each agent interacts with all eight adjacent cells).
- Disease Progression: Agents belong to one of four groups marking disease progression:
    - Susceptible (S): No exposure, vulnerable to CWD.
    - Exposed (E): Exposed, incapable of transmitting.
    - Infected (I): Infectious, capable of transmitting.
    - Clinical (C): Infectious and showing symptoms, capable of transmitting.
- Environmental Disease Reservoir (V): Prions present in  environment, can infect susceptible agents through indirect transmission. Treat as spatially distributed disease pool on the grid.

During each time step (per week):
- Transmit CWD to suceptible agents via direct and indirect transmission
- For Exposed, Infected, Clinical agents, transition disease state
- Simulate births and (non-prion) deaths within populations
- Infected and Clinical agents shed prions
- Envionmental prions decay

Initial conditions: 
- Populations (set N = 1500, vary distribution)
    - Infected deer in middle of grid (I = n)
    - No infected deer, nonzero environmental prion load (V = n)
- Refer to Table 1 for parameter values
- To consider: two populations of deer

- Remove infected deer from population
    - Try at different states of simulation (1 year, 5 years, 10 years, etc)
    - Analogous to culling program early vs late
- Introduce physical barriers on grid, analogous to isolation efforts

Data to Collect:
- Population over time, sub aggregate per SEIR state
- Prion load: environmental vs direct (visualize using heatmap)
- Death rate per year over 50-100 year simulation
- Fraction of original population size
- Proportion of deer infected by direct vs indirect transmission
- Compare birth rate to death rate to calclate mean growth rate
- Calculate reproductive number R0

Agent Behavior:
- Agent Movement: Deer move along grid via a random walk (To consider: movement could be altered to avoid clinical deer, i.e. deer with visible symptoms)
- Diease Transmission: 
    - Direct: Infected (I) and Clinical (C) agents directly transmit the disease to neighboring Susceptible (S) agents with a certain probability.
        - Rate determined by the strength of infection $\beta_d$ or dependent on $\lambda$.
        - To consider: Distance between infected and suceptible agents: i.e. Agents on diagonals are less likely to be infected compared to adjacent deer?
    - Indirect: Susceptible agents become infected by coming into contact with enviromental proins (V) present in residing grid.
        - Prion concentration to determine probability for transmission
- Agent Classes: SEIC classes divided into subclasses
    - E: 35 sub classes, transition probability equal to 1 for $E_i\to E_{i+1}$ $i = 1,\dots,15$, and $1-\rho_i$ for $E_{i+15}\to E_{i+16}$ $i = 1,\dots,20$.
    - I: 44 sub classes, transition probability equal to 1 for $I_i\to I_{i+1}$ $i = 1,\dots,23$, and $1-\sigma_i$ for $I_{i+23}\to I_{i+24}$ $i = 1,\dots,20$.
    - C: 35 sub classes, transition probability $1-\mu_i$ for $C_{i}\to C_{i+1}$ $i = 1,\dots,34$.
- Calculate prion load on Infected and Clinical deer

Environmental Prions:
- Shedding: Infected (I) and Clinical (C) deer shed disease into the environment at a rate $\tau$ per agent per time step.
    - Spatially distributed: prions spread across the grid based on agent's position and movement
- Decay: Prions decay probabilistically over time, with a weekly decay rate $\phi$. To consider: more complex decay, e.g. concentration/time/spatially dependent.
- Environmental Concentration: Prion concentration in a given grid cell is updated at each time step according to shedding rate of nearby agents and decay rate.

Parameters:
- $\tau$: weekly shedding rate.
- $\phi$: weekly decay rate.
- $\beta_d$ and $\beta_i$: transimission probabilities.
- $rho_i$, $\sigma_i$, $\mu_i$: transition probabilities.