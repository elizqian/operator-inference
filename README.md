<style>

.myDiv:target{

display:block !important;

}

</style>

# Operator Inference

This code implements the operator learning approach described in

1. Peherstorfer, B. and Willcox, K. 
[Data-driven operator inference for non-intrusive projection-based model reduction.](https://www.sciencedirect.com/science/article/pii/S0045782516301104)
Computer Methods in Applied Mechanics and Engineering, 306:196-215, 2016.
([Download](https://cims.nyu.edu/~pehersto/preprints/Non-intrusive-model-reduction-Peherstorfer-Willcox.pdf))<details><summary>BibTeX</summary><pre>
@article{Peherstorfer16DataDriven,
    title   = {Data-driven operator inference for nonintrusive projection-based model reduction},
    author  = {Peherstorfer, B. and Willcox, K.},
    journal = {Computer Methods in Applied Mechanics and Engineering},
    volume  = {306},
    pages   = {196-215},
    year    = {2016},
}</pre></details>

2. Qian, E., Kramer, B., Marques, A., and Willcox, K. 
[Transform & Learn: A data-driven approach to nonlinear model reduction](https://arc.aiaa.org/doi/10.2514/6.2019-3707).
In the AIAA Aviation 2019 Forum, June 17-21, Dallas, TX. ([Download](https://www.dropbox.com/s/5znea6z1vntby3d/QKMW_aviation19.pdf?dl=0))<details><summary>BibTeX</summary><pre>
@inbook{QKMW2019aviation,
author = {Elizabeth Qian and Boris Kramer and Alexandre N. Marques and Karen E. Willcox},
title = {Transform \&amp; Learn: A data-driven approach to nonlinear model reduction},
booktitle = {AIAA Aviation 2019 Forum},
doi = {10.2514/6.2019-3707},
URL = {https://arc.aiaa.org/doi/abs/10.2514/6.2019-3707},
eprint = {https://arc.aiaa.org/doi/pdf/10.2514/6.2019-3707}
}</pre></details>

The `examples` folder contains scripts that set-up and run several examples:
* The heat equation example from [1].
* The Burgers' equation from [1].
* The Euler equation example from [2]. This example uses MATLAB's Curve Fitting Toolbox to generate the random initial conditions.

This code was developed in MATLAB R2019a and is in a beta state. Please contact Elizabeth Qian (elizqian@mit.edu) with any questions or comments.