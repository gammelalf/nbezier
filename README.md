# nbezier
nbezier aims to be a general purpose library for working with bezier curves of any degree.

It is backed by [nalgebra](https://nalgebra.org/) (hence the name) and optimised for linear, quadratic and cubic curves.

## Current Features
- store curve as list of control points
- [de Castlejau's Algorithm](https://en.wikipedia.org/wiki/De_Casteljau%27s_algorithm)
  - evaluate a point
  - split a curve
- a curve's polynomial and its derivative
  - normal and tangental vectors
  - minimal axis aligned bounding box (currently only for cubic or less)
- a curve's control points' axis aligned bounding box and convex hull

## Experimental Features
- find a point on a curve
- find all intersection points between two curves

## Planned Features
- degree elevation and reduction
- any suggestions?

## Optimization
As this project was originally intended to power a wasm web app, it is optimised for cubic and lower curves
which are the only ones directly supported by svg or canvas.

**But what exactly is optimised?**
1. [smallvec](https://github.com/servo/rust-smallvec) is used to avoid heap allocations for said degrees
2. A lot of magic constants are used which are precalculated from the general formulas
