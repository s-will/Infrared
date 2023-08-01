##
# @namespace infrared.libinfrared
# @brief The libinfrared module exposing C++ infrared functionality to Python
#
# The module libinfrared is the Python interface to the low level C++
# library. It is specified using pybind11 in libpython.cpp.


# NOTE: here, we define dummy classes/functions only for documentation.
# These entities are defined in libinfrared.cpp as exports from C++
#

def seed():
    """
     @brief seed the C++-side random number generator
     Args:
        x: integer seed value
    """
    pass

class FiniteDomain:
    """
    @brief A finite domain of a variable

    Defines lower and upper bound of a contiguous domain for a finite domain
    variable.

    Interfaces the C++ class ired::FiniteDomain. 
    Supports (deep) copying.
    """
    
    def __init__(self,*args):
        """
        Construct finite domain
        
        Args:
            args: size or (lb,ub) construct as 0..size-1 or lb..ub, resp.
        """
        pass

    def lb(self):
        """
        Lower bound

        Returns: lower bound
        """
        pass

    def ub(self):
        """
        Upper bound

        Returns: upper bound
        """
        pass

    def size(self):
        """
        Domain size

        Returns: size ub-lb+1 or zero if empty
        """
        pass

    def empty(self):
        """
        Test for empty domain
        
        Returns: true, if empty
        """
        pass

    def contains(self,val):
        """
        Test for membership

        Note: `in' works as alias
        
        Returns: true, if val is in the domain
        """
        pass

    def undet(self,val):
        """
        Value to flag 'undetermined'
        
        for internal use

        Returns: undet value specific for this domain
        """
        pass


class Function:
    """
    @brief A real valued feature function

    Typically not used direclty, but through wrappers.

    Interfaces the C++ class ired::Function &lt;double&gt;. Exposes constructor as
     __init__, operator() as __value__, and vars
    """
    pass

class IntFunction:
    """
    @brief An integer valued feature function

    Typically not used direclty, but through wrappers.
    Used by the system to hold a scaled and rounded version of a feature function for optimization.

    Interfaces the C++ class ired::Function &lt;int&gt;. Exposes constructor as
    _init__, operator() as __value__, and vars
    """
    pass

class Constraint:
    """
    @brief A constraint

    Interfaces the C++ class ired::Function &lt;bool&gt;
    Exposes constructor as __init__, operator() as __value__, and vars
    """
    pass

class Assignment:
    """
    @brief An assignment (of values to variables)

    Holds solutions as returned by solvers.

    Interfaces the C++ class ired::Assignment.
    """

    def values():
        """Access values of the assignment"""
        pass
    
    pass
