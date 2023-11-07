# GyB_TDK_2023
Scripts used for the paper "Utilizing the memory effect of elastic tires in the steering of a self-balancing e-scooter", which was presented in BME Students' Scientific Conference, 2023.

## Abstract
The main inspiration of this paper was the self-balancing motorcycle developed by Honda [1], which was first unveiled in 2017. The motor is able to stabilize itself when standing still (with no longitudinal motion) using the Honda Riding Assist technology. The motorcycle is balanced by the steering assist system, instead of using gyroscopes, as other companies did previously. This industrial development area inspired the topic of this work, because there is an ongoing research at BME Department of Applied Mechanics to propose similar control strategies for a self-balancing vehicle.

To test the theoretical control laws, first, an experimental 2-wheeled vehicle needs to be developed. In this paper, we design a power steering system for a Xiaomi Mi Electric Scooter 3 [2], by which the fork-handlebar assembly of the scooter can be steered using a brush-less DC (BLDC) motor. Therefore, it can be used for self-balancing tasks.

In our former work [3], we focused on estimating the coefficient of friction at the steered wheel. Knowledge about the frictional relation between the road surface and the elastic tire would provide important information for the steering control, which is necessary for balancing the electric scooter. The above-specified experimental set-up made it possible to design measurements and experiments to further enhance the knowledge about the effect of friction forces on the accuracy of steering control.

Based on measurement data, we aim to enhance the accuracy of our delayed tire model, which describes the relation between the ground and the elastic tire. The simulation results are compared to the acquired experimental data. Hence, we lay out the foundation of a further improved friction estimation method.
### Reference
[1] Honda, “Honda riding assist,” https://global.honda/innovation/CES/2017/002.html, 2017, [Online; accessed: 2023-09-09].

[2] Xiaomi, “Mi electric scooter 3,” https://xiaomiofficial.hu/okoseszkoz/elektromos-rollerek/mi-electric-scooter-3/, [Online; accessed: 2023-09-09].

[3] B. M. Györök, “Estimating the coefficient of friction based on the memory effect of an elastic tire,” Budapest University of Technology and Economics, Students’ Scientific Conference, 2022.

# User guide
- In the *SimulationEnv* directory, the closed-loop and open-loop simulation files can be found. Running the `MAIN_sim.m` file various steering control algorithms can be tested while running the `OpenLoop.m` file results in an open-loop simulation based on measurement inputs.
- The *ODrive_experiments* directory contains the jupyter notebooks for different experiments with the self-balancing electric scooter.
- The *ANN_ErrorPrediction* directory contains the training and test files and the script for the ANN-based conditional estimation.


