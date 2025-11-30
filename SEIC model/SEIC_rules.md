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

Agent Behavior:
- Agent Movement: Deer move along grid via a random walk (To consider: movement could be altered to avoid clinical deer, i.e. deer with visible symptoms)
- Diease Transmission: 
    - Direct: Infected (I) and Clinical (C) agents directly transmit the disease to neighboring Susceptible (S) agents with a certain probability.
        - Rate determined by the strength of infection $\beta_d$ or dependent on $\lambda$.
        - To consider: Distance between infected and suceptible agents: i.e. Agents on diagonals are less likely to be infected compared to adjacent deer?
    - Indirect: Susceptible agents become infected by coming into contact with enviromental proins (V) present in residing grid.
        - Prion concentration to determine probability for transmission
- Agent Classes: SEIC classes divided into subclasses
    - E: 35 sub classes, transition probabilities 

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