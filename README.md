# libtwobody - two body orbital mechanics

[![Build Status](https://travis-ci.org/rikusalminen/twobody.svg?branch=master)](https://travis-ci.org/rikusalminen/twobody)
[![Coverage Status](https://coveralls.io/repos/rikusalminen/twobody/badge.svg?branch=master&service=github)](https://coveralls.io/github/rikusalminen/twobody?branch=master)
[![Coverity Status](https://scan.coverity.com/projects/6810/badge.svg)](https://scan.coverity.com/projects/rikusalminen-twobody)

libtwobody is a software library for two body orbital mechanics.
It can be used to determine orbits and predict the motion of satellites,
planets and moons.
libtwobody accurately and realistically solves the two body problem of two
orbiting massive particles in an unperturbed orbital motion with no external
forces.

## Features
* Elliptic, hyperbolic and parabolic trajectories
* Classical and universal variable formulations
* Determine orbital elements from position and velocity vectors
* Advanced solver for time of flight equations (Kepler's equation,
    Gudermannian, Barker's equation, Universal Kepler's equation)
* Predict position and velocity vectors (and other quantities) at any point in
    time using true anomaly, eccentric/hyperbolic/parabolic anomaly or
    universal variables
* Find closest approaches between satellites or sphere of influence
    transitions on lunar/interplanetary trajectories

## Tests

libtwobody is extensively tested with a purpose-built test framework (called
numtest).
The test framework generates test cases from a "seed" (64 bit integer) and
verifies physical properties (e.g. conservation of energy and angular
momentum) and mathematical invariants (e.g. focus-directrix property of conic
sections) apply.
Test cases are repeatable and reproducible by using the same seed value.
Particular attention is paid to floating point issues (NaN, inf).
Half of libtwobody exists to verify that the other half works correctly.
Running the tests takes tens of minutes of CPU time.

## Bibliography

* Bate, Mueller, White: Fundamentals of Astrodynamics
* Conway, Prussing: Orbital Mechanics
* Battin, R.H.: Introduction to the Mathematics and Methods of Astrodynamics
* Danby, J.M.A.: Fundamentals of Celestial Mechanics
* Conway, B.: An improved algorithm due to Laguerre for the
    solution of Kepler's equation
* Danby, J.M.A.: The solution of Kepler's Equation, Parts I, II, III
* Hoots, Crawford, Roehrich: An analytic method to determine future close
    approaches between satellites
* Rodriguez, Fadrique, Klinkrad: Collision risk assessment with a "smart
    sieve" method
