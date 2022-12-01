# nbezier
<a href="https://github.com/gammelalf/nbezier/blob/master/LICENSE">
  <img src="https://img.shields.io/github/license/gammelalf/nbezier" alt="license">
</a>
<a href="https://crates.io/crates/nbezier">
  <img src="https://img.shields.io/crates/v/nbezier" alt="crates.io">
</a>
<a href="https://docs.rs/nbezier/0.2.1/nbezier/">
  <img src="https://img.shields.io/docsrs/nbezier" alt="docs">
</a>

nbezier aims to be a general purpose library for working with bezier curves of any degree.

It uses [nalgebra](https://nalgebra.org/) (hence the name) to implement a generic `BezierCurve`.

This library also provides a non-generic type `SimpleCurve` hiding nalgebra's complexity.
`SimpleCurve` is optimised for cubic bezier curves while supporting arbitrary degree.

## Current Features
- store curve as list of control points (a matrix' column vector to be precise)
- [de Castlejau's Algorithm](https://en.wikipedia.org/wiki/De_Casteljau%27s_algorithm)
  - evaluate a point
  - split a curve
- compute a curve's polynomial and its derivative
  - normal and tangental vectors
- a curve's control points' axis aligned bounding box and convex hull
- raise or reduce a curve's degree

## Experimental Features
- find a point on a curve
- find all intersection points between two curves

## Planned Features
- any suggestions?

## How is `SimpleCurve` optimised?

Using nalgebra `BezierCurve` is generic over its degree.
`SimpleCurve` is an enum storing curves of degree 1, 2, 3 and anything above in its 4 variants.
Since these low degrees are their own variant with dedicated type, rust can monomorphize these
computing a lot of "magical constants" at compile time.
Also these degree's variants are stored exclusively on the stack
giving them an enormous performance boost.
