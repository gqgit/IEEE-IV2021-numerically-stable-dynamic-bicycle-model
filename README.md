# IEEE-IV2021-numerically-stable-dynamic-bicycle-model
Numerically Stable Dynamic Bicycle Model for Discrete-time Control, accepted by IEEE IV2021 workshop  
Authors: Qiang Ge, Qi Sun, Shengbo Eben Li, Sifa Zheng*, Wei Wu and Xi Chen  
Affiliation: School of Vehicle and Mobility, Tsinghua University  
Contact: gq17@mails.tsinghua.edu.cn

Abstract: Dynamic/Kinematic model is of great significance in decision and control of intelligent vehicles. However, due to the singularity of dynamic models at low speed, kinematic models have been the only choice under such driving scenarios. Inspired by the concept of backward Euler method, this paper presents a discrete dynamic bicycle model feasible at any low speed. We further give a sufficient condition, based on which the numerical stability is proved. Simulation verifies that (1) the proposed model is numerically stable while the forward-Euler discretized dynamic model diverges; (2) the model reduces forecast error by up to 65% compared to the kinematic model. As far as we know, it is the first time that a dynamic bicycle model is qualified for urban driving scenarios involving stop-and-go tasks.  

For readers/users: Thank you for your attention! You may find the discrete model named 'dynamic_linear_backward(x0,u0,T)' in dynamic_vs_kinematic.py or forward_vs_backward.py. You may run dynamic_vs_kinematic.py to check the accuracy, or run forward_vs_backward.py to check the numerical stability. Variable 'speed' and steplength 'T' could be adjusted to observe the trends.

