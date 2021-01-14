![GitHub](https://img.shields.io/github/license/fneum/benchmark-teplopf)

![GitHub repo size](https://img.shields.io/github/repo-size/fneum/benchmark-teplopf)

# Benchmark of PyPSA Transmission Expansion Planning (TEP)

This repository contains the snakemake workflow to reproduce the benchmarks from the paper

**Fabian Neumann and Tom Brown. 2020. Transmission Expansion Planning Using Cycle Flows. In Proceedings of the Eleventh ACM International Conference on Future Energy Systems (e-Energy ’20).**

> The common linear optimal power flow (LOPF) formulation that underlies most
transmission expansion planning (TEP) formulations uses bus voltage angles as
auxiliary optimization variables to describe Kirchhoff's voltage law. As well
as introducing a large number of auxiliary variables, the angle-based
formulation has the disadvantage that it is not well-suited to considering the
connection of multiple disconnected networks, It is, however, possible to
circumvent these auxiliary variables and reduce the required number of
constraints by expressing Kirchhoff's voltage law directly in terms of the
power flows, based on a cycle decomposition of the network graph. In
computationally challenging benchmarks such as generation capacity expansion
with multi-period LOPF, this equivalent reformulation was shown in previous
work to reduce solving times for LOPF problems by an order of magnitude.
Allowing line capacity to be co-optimized in a discrete TEP problem makes it a
non-convex mixed-integer problem. This paper develops a novel cycle-based
reformulation for the TEP problem with LOPF and compares it to the standard
angle-based formulation. The combinatorics of the connection of multiple
disconnected networks is formalized for both formulations, a topic which has
not received attention in the literature. The cycle-based formulation is shown
to conveniently accommodate synchronization options. Since both formulations
use the big-$M$ disjunctive relaxation, useful derivations for suitable big-$M$
values are provided. The competing formulations are benchmarked on a realistic
generation and transmission expansion model of the European transmission system
at varying spatial and temporal resolutions. The cycle-based formulation solves
up to 31 times faster for particular cases, while averaging at a speed-up of
factor 4.

For a benchmark of PyPSA Linear Optimal Power Flow (LOPF) formulations, see https://github.com/FRESNA/benchmark-lopf.

## Cite As

```
@inproceedings{NeumannTEP2020,
    author = {Neumann, Fabian and Brown, Tom},
    title = {Transmission Expansion Planning Using Cycle Flows},
    year = {2020},
    booktitle = {Proceedings of the Eleventh ACM International Conference on Future Energy Systems},
    series = {e-Energy ’20}
}
```
