# Theory-of-Mechanisms
This is the repository for the programs used for the course ME6221 (Theory of Mechanisms) to design, analyse and synthesize the [Jansen Linkage](https://en.wikipedia.org/wiki/Jansen%27s_linkage)

## File structure
```python
.
├── data
├── outputs                                 # All the final plots are saved here
│   ├── dimension_synthesis
│   ├── dynamic_analysis
│   └── kinematic_analysis
├── README.md
└── scripts
    ├── dimension_synthesis_jansen.ipynb    # Program for dimension synthesis 
    ├── dynamic_analysis_jansen.ipynb       # Program for dynamic analysis
    ├── jansen_linkage.py
    ├── kinematic_analysis_jansen.ipynb     # Program for kinematic analysis
    └── tom_utils.py

```
## Dimension synthesis
Dimension synthesis is the process of determining the dimensions of the components of a mechanism. The following plots show parts of the mechanism assembled in 4 key positions.
<p float="left">
  <img src="outputs/dimension_synthesis/dyad_assembly.png" width="300"/>
  <img src="outputs/dimension_synthesis/four_bar.png" width="300"/>
</p>


## Kinematic analysis
The objective of kinematic analysis is to determine the position, velocity and acceleration of the links of a mechanism.
<p float="left">
  <img src="outputs/jansen_linkage.gif" width="500" title="Kinematic animation of the Jansen Linkage"/>
</p>
<p float="left">
  <img src="outputs/kinematic_analysis/position_angle_plots.png" width="250" title="Position vs time plot of each joint"/>
  <img src="outputs/kinematic_analysis/velocity_angle_plots.png" width="250" title="Velocity vs time plot of each joint"/>
    <img src="outputs/kinematic_analysis/acceleration_angle_plots.png" width="250" title="Acceleration vs time plot of each joint"/>
</p>

## Dynamic analysis
The objective of kinematic analysis is to determine the reaction forces at the joints and input torques when a mechanism is subject to external forces.
*(Naming conventions are presented in data/jansen_dynamics.jpg)*
<p float="left">
  <img src="outputs/dynamic_analysis/center_of_gravity.png" width="500" title="Plot of the center of gravity of the system and the links"/>
</p>
<p float="left">
  <img src="outputs/dynamic_analysis/walking_force_plots.png" width="500" title="The reaction forces at the joints of the mechanism, when subjected to a ground reaction force"/>
</p>


