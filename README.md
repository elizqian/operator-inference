# Operator Inference

This code implements the operator learning approach described in

1. Peherstorfer, B. and Willcox, K. [Data-driven operator inference for non-intrusive projection-based model reduction.](https://cims.nyu.edu/~pehersto/preprints/Non-intrusive-model-reduction-Peherstorfer-Willcox.pdf) Computer Methods in Applied Mechanics and Engineering, 306:196-215, 2016.
2. Qian, E., Kramer, B., Marques, A., and Willcox, K. 
[Transform & Learn: A data-driven approach to nonlinear model reduction](https://arc.aiaa.org/doi/10.2514/6.2019-3707).
In the AIAA Aviation 2019 Forum, June 17-21, Dallas, TX. ([Download](https://www.dropbox.com/s/5znea6z1vntby3d/QKMW_aviation19.pdf?dl=0))

The `examples` folder contains scripts that set-up and run several examples:
* The heat equation example from [1]
* The Burgers' equation from [1]
* The Euler equation example from [2]. This example uses MATLAB's Curve Fitting Toolbox to generate the random initial conditions.

This code was developed in MATLAB R2019a and is in a beta state. Please contact Elizabeth Qian (elizqian@mit.edu) with any questions or comments.