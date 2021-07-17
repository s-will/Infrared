#ifndef INFRARED_FINITEDOMAIN_HPP
#define INFRARED_FINITEDOMAIN_HPP

/*
 * InfraRed ---  A generic engine for Boltzmann sampling over constraint networks
 * (C) Sebastian Will, 2018
 *
 * This file is part of the InfraRed source code.
 *
 * InfraRed provides a generic framework for tree decomposition-based
 * Boltzmann sampling over constraint networks
 */

/**
 * @file
 *
 * @brief Defines finite domains 
 */

namespace ired {

    /**
     * Specification of a finite domain
     *
     * @note Even if currently only contiguous domains are supported,
     * users shouldn't assume this, but use the provided methods to
     * manipulate domain values.
     */
    class FiniteDomain {
        public:
            /**
             * @brief construct as domain 0..(domain size - 1)
             * @param domsize size of the domain
             */
            explicit
            FiniteDomain(int domsize = 0)
                : lb_(0),
                  ub_(domsize-1) {
            }
            
            /**
             * @brief construct as contiguous domain lb..ub
             * @param bounds lower and upper bound
             *
             * @note both, lower and upper bound, are included in the domain
             */
            explicit
            FiniteDomain(std::pair<int,int> bounds)
                : lb_(bounds.first),
                  ub_(bounds.second) {
                assert(lb_ <= ub_);
            }

            /**
             * @brief construct as contiguous domain lb..ub
             * @param lb lower bound
             * @param ub upper bound
             *
             * @note both, lower and upper bound, are included in the domain
             */
            FiniteDomain(int lb, int ub)
                : lb_(lb),
                  ub_(ub) {
                assert(lb_ <= ub_);
            }

            /**
             * @brief size of the domain
             * @return size
             */
            int
            size() const {
                return ub_ - lb_ + 1;
            }
            
            /**
             * @brief lower bound
             * @return bound
             */
            int lb() const {
                return lb_;
            }

            /**
             * @brief upper bound
             * @return bound
             */
            int ub() const {
                return ub_;
            }

            /**
             * @brief increment domain value
             * @param v a domain value
             *
             * Sets v to the next larger domain value after v
             * or the lowest domain value if v was undetermined.
             *
             * @note undefined behavior if v not in the domain
             * @note if v was the last domain variable, then
             * afterwards, in(v) is false and v exceeds the upper bound.
             */
            void
            inc(int & v) const {
                assert( lb_<= v && v < ub_ );
                v++;
            }

            /**
             * @brief undetermined
             * @return domain-specific undetermined value 
             */
            int
            undet() const {
                return lb() - 1;
            }

            /**
             * @brief test domain membership
             * @param value a domain value
             * @return whether v is in the domain 
             */
            bool
            in(int v) const {
                return lb_ <= v and v <= ub_;    
            }

        private:
             int lb_; //!< lower bound
             int ub_; //!< upper bound
    };

    using FiniteDomains = std::vector<FiniteDomain>;

}
#endif

