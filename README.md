# Project-Thesis
Code for project thesis

**BENCH**
A bench is the frame that is simulated at each timestep.
The user should specify:
- Screen position
- Number of rays
- Surface position



**Using fminunc**
- Requires optimization toolbox, see Rays.m 424

% The objective function and its derivatives must be of type double
% because optimization algorithms typically rely on numerical stability
% and precision that double-precision floating-point numbers provide.