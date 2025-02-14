# TPCFillLevel
Utility for calculating the mass filled and fill level in a time projection chamber of a given geometry.

## Volume calculation
The fluid volume is the sum of all (positive) void volumes and (negative) solid volumes in the detector. For each of these subvolumes we can define the cumulative volume $V(y)$ up to a height $y$:

$$V(y) = \int_0^{y}A(y^\prime)\,dy^\prime$$

where $A(y)$ is the *horizontal* cross-sectional area of the subvolume at height $y$. Rather than computing $A(y)$ by doing a two-dimensional integral at each $y$, we can instead define the *width profile*, $w(y)$, as the width of a rectangle with a length $\ell$ and an equal area:

$$w(y)\equiv \frac{A(y)}{\ell}$$

The cumulative volume is calculated by integrating over $w(y)$ only:

$$V(y) = \ell\int_0^{y}w(y^\prime)\,dy^\prime$$

This simplifies the calculation significantly, as many different objects can be described using only a few common width profiles. For example, a circular width profile describes horizontal cylinders, cones, or spheres, and a rectangular width profile describes vertical cylinders or prisms of any kind. With the width profile specified, the user only has to provide the total volume (calculated using Solidworks or from drawings) so that $\ell$ can be determined.

## Defining the geometry
A detector is defined in its own YAML file in [geometries](https://github.com/clarkehardy/TPCFillLevel/tree/main/geometries). The YAML file contains a key for each of the components in the detector. Each component has a set of attributes that must be filled in:
- `volume`: `float`, the total volume of the object or void
- `height`: `float`, the maximum height of the object or void
- `y_position`: `float`, the vertical coordinate defining the lowest point on the object
- `number`: `int`, the number of instances for objects that are repeated multiple times
- `spacing`: `float`, the vertical spacing between the lowest point of one instance and the lowest point of the next (zero if they are all at the same level)
- `profile`: `string`, which width profile to use for the object. Currently, `rectangular` and `circular` are supported, but other profiles can be added by defining a new functions in [FillLevel.py](https://github.com/clarkehardy/TPCFillLevel/blob/main/FillLevel.py)
- `solid`: `bool`, whether the object is solid (`True`) or empty space (`False`)

## Usage
The [example notebook](https://github.com/clarkehardy/TPCFillLevel/blob/main/example.ipynb) should be sufficient to get started using the code. After initializing a `FillLevel` object with a YAML file, the `build_fill_funcs()` method will build the functions which allow for the mass and fill levels to be calculated. Methods which return the mass filled, fill level, or make plots for visualization can then be called.
