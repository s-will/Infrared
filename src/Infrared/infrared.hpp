#ifndef INFRARED_HPP
#define INFRARED_HPP

/**
 * @mainpage
 * InfraRed -- A generic engine for Boltzmann sampling over constraint networks
 *
 * (C) Sebastian Will, 2018
 *
 * InfraRed provides a generic framework for tree decomposition-based
 * Boltzmann sampling over constraint networks
 *
 * Design goals: engine allows specifying function/constraint
 * annotated tree decomposition in Python (with functions/constraints
 * encoded in either C++ or Python); evaluate partition functions; and
 * sample from Boltzmann distribution
 *
 * Keep it simple:
 *   -- variables are indexed from 0..n-1; using index type int
 *   -- variable domains are always 0..k-1 (where domains can have different size); value type is int
 *
 * Keep it flexible:
 *   -- function value type templated
 *   -- function values are combined according to policy (supports evaluation algebras)
 *   -- functions / and constraints are polymorph
 *
 * Issue: Since functions and constraints are polymorph, we need to store pointers, references.
 * Policy: we guarantee persistence of the objects passed to the framework. The CN owns
 * functions and constraints. All other classes hold references.
 * The user can decided whether to transfer ownership (unique pointers) or require cloning.
 * Cloning must be implemented for the polymorph function and constraint classes.
 *
 * The ClusterTreeDecomposition holds a constraint network; when
 * manually constructed, it autosyncs with its CN
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
