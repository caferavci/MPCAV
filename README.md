# MPCAV
Quadratic Programming and Model Predictive Control Based Trajectory Optimization for platoons of equipped Cooperative Adaptive Cruise Control (CACC) Vehicles

## Authors

1. Cafer Avci - Transportation Engineering, Aalto University, Espoo, Finland

2. Konstantinos Mattas - European Commission Joint Research Centre. Ispra, Italy

3. Giovanni Albano - European Commission Joint Research Centre. Ispra, Italy

4. Claudio Roncoli - Transportation Engineering, Aalto University, Espoo, Finland

5. Biagio Ciuffo - European Commission Joint Research Centre. Ispra, Italy


## Motivation

The objective of this research is to develop and test a trajectory optimisation strategy for platoons of CACC-equipped vehicles that is:

* able to optimise trajectories cooperatively for all vehicles in a platoon, by considering vehicle dynamics based on state-of-the-art CACC models, including also desired system-level objectives;

* robust not only in fully automated environment but also in mixed vehicle traffic environment;

* able to work on-line in different traffic conditions

* able to adapt headway to optimal macroscobic headway obtained from upper-level optimization layer.

## Methodology

A platoon is described by a sequence of vehicles where each vehicle feature sensors capable of measuring distance and speed from the preceding (and, in some cases, also the following) one, while being also capable of communicating its control actions, e.g., accelerations, via V2V (vehicle-to-vehicle) communication. In particular, we assume that all vehicles in a platoon, with exception of the leading vehicle, i.e., the first vehicle in the platoon, implement CACC as the longitudinal controller. This study adopts the model by Van Arem et al., based on previous simulation studies focusing on deterministic acceleration, to calculate the acceleration of CAVs. This model has been chosen as, with appropriate parameter choice, guarantees smooth responses and string stability. 

Van Arem, B., Van Driel, C. J., & Visser, R. 2006. The impact of cooperative adaptive cruise control on traffic-flow characteristics. IEEE Transactions on intelligent transportation systems, 7(4), 429-436.

![alt text](https://github.com/caferavci/MPCAV/blob/main/Media/CACC_Model.jpg)

As previously stated, the QP problem is implemented in an MPC framework, in order to reject any past inaccuracies and to maintain the difference between the model predictions, e.g., the HDV trajectory, and the real process outcome at low levels. 

![alt text](https://github.com/caferavci/MPCAV/blob/main/Media/MPC_Fig.jpg)

The idea behind MPC is sketched in Figure 1. The procedure for implementing the controller follows the following steps:
1. The platoon leader detects the distance and speed of its HDV predecessor at time k; platoon-level information is exchanged via V2V among vehicle in the platoon.
2. A prediction of the HDV trajectory is made for horizon [k,k+t_hor), by using current and historical information on its movement.
3. The optimisation problem presented above is solved for time horizon [k,k+t_hor).
4. The optimal control action is implemented for the duration of the control step, e.g., [k,k+1).

We utilise the optimisation software Gurobi (27) for solving the optimisation problem. Simulation codes written are written in C++, and experiments are performed with Aimsun traffic modelling software. Optimization framework dynamically link to Aimsun API so that MPC can be executed in every control time step. After each control time step, predicted trajectory HDV is calculated for time horizon by averaging  using current and past real control time step movements which are exchanged with Aimsun API and optimization framework. Proposed algorithm and simulation environment produce results within very short time because of the optimisation approach structure and selected programming language. Computation times can be change initial settings of MPC and simulation environment. 

## Experiments

### Behaviours of CAVs in platoon level

In this section, we perform two experiments to demonstrate how to optimize vehicle trajectories to show the effectiveness of the concept. To compare effectiveness of developed model, benchmark has been performed with ACC model (Milanés and Shladover, 2014).

Milanés, V., Shladover, S.E., 2014. Modeling cooperative and autonomous adaptive cruise control dynamic responses using experimental data. Transportation Research Part C: Emerging Technologies 48. doi:10.391016/j.trc.2014.09.001.

We simulate a platoon of 5 CAVs following an HDV. To simulate behaviors of CAVs, wedesigned two different scenarios as follows:

* In the first scenario, the HDV brakes suddenly with minimum deceleration rate in the1middle of simulation from 108 km/h to 36 km/h.
* In the second scenario, the HDV accelerates with maximum rate suddenly in the middle3of simulation from 36 km/h to 108 km/h.


