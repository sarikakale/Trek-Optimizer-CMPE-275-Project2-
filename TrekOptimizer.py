import TrekOptimizer
# In TrekOptimizer.cpp we expose hello() function, and it now exists in the TrekOptimizer module.
assert 'trekPlanner' in dir(TrekOptimizer)
# TrekOptimizer.hello is a callable.
assert callable(TrekOptimizer.trekPlanner)
# Call the C++ hello() function from Python.
print TrekOptimizer.trekPlanner()

