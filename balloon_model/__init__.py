""" Physics-based balloon forward model.

    All physics is written as vectorised NumPy functions with no Parcels
    dependency, so the same code serves the standalone column integrator and
    the Parcels kernel wrapper. Equation labels in comments (e.g. [E-DYN-1])
    name the equations the code implements.
"""
