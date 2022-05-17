# Interactive 1D basin model

Return orbits of stars with data from the Gaia mission.
The query gets sky positions, distances (parallaxes), proper motions, and radial velocities 
for a set of stars that are close to the Sun. 
The observed, heliocentric kinematic measurements are transforme dto Galactocentric Cartesian coordinates 
and we use the positions and velocities as initial conditions to compute the orbits of these stars in the galaxy using the gala Python package.

```python
import pydaisi as pyd
orbits_query = pyd.Daisi("Gaia 2 Orbits Query")
result = orbits_query.get_orbits(nb = 4096, 
                                parallax_over_error = 8, 
                                parallax = 10, 
                                BP_RP_range = [0.5, 0.7], 
                                M_G_range = [2.0, 3.75], 
                                time_step = 1., 
                                time = 500.).value
```