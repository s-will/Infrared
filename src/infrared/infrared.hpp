#ifndef INFRARED_HPP
#define INFRARED_HPP

/**
 * @mainpage
 * InfraRed -- A generic framework for Dynamic Programming, Optimization and Boltzmann 
 * sampling over constraint models
 *
 * (C) Sebastian Will, 2018--2021
 *
 * InfraRed provides a generic framework for tree decomposition-based
 * Dynamic Programming, including optimization and Boltzmann sampling
 * based on constraint models 
 *
 * Declarative specification of constraint models which are automatically
 * transformed to function/constraint annotated tree decompositions in
 * Python. Specific functions and constraints can be encoded in Python. C++
 * engine optimizes; evaluates partition functions; and samples from
 * Boltzmann distributions. Python code extends such functionality to
 * Multi-dimensional Boltzmann sampling.
 *
 * Keep it simple, but flexible on the C++ side
 *
 * - variables are indexed from 0..n-1; using index type int
 *
 * - variable domains are consecutive lb...ub; value type is int
 *
 * - function value type templated
 *
 * - function values are combined according to policy, which supports evaluation
 *   algebras
 *
 * - support sparse, memory-optimized data structures for intermediary
 *   results to enable larger problem instances
 *
 * Provide a flexible interface and added functionality on the Python side
 *
 * - support constraint modeling syntax to declarativly and compositionally
 *   describe constraint problems
 *
 * - support definintion of 'custom' constraints and functions in Python
 * 
 * - support efficient evaluation of models for optimization and Boltzmann sampling
 */

/**
 * @file
 * @brief Includes all headers of the (header-only) Infrared core library
 */

// package parts
#include "assignment.hpp"
#include "functions.hpp"
#include "cluster.hpp"
#include "constraint_network.hpp"
#include "cluster_tree.hpp"

#endif
