title: Adding New Models

# Developer Guide

The library is highly extensible thanks to the use of abstract base classes.

## How to add a new Hardening Law
1. Create a new module in `src/plasticity/hardening_laws/`.
2. Define a type that extends `Base_hardening_laws`.
3. Implement the deferred procedures:
    - `stress(ep)`
    - `dstress_dep(ep)`
    - `ddstress_ddep(ep)`

## How to add a new Yield Criterion
1. Extend `Base_yield_critera`.
2. Implement `stress_eq(stress)`.
3. (Optional) Implement analytical derivatives for higher performance, or use the default numerical ones.