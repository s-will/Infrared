#ifndef BSCF_SAMPLE_HPP
#define BSCF_SAMPLE_HPP


/**
 * Constraint-Network Boltzmann sampling engine
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
 * The ClusterTreeDecomposition holds a constraint network; when manually constructed, it autosyncs with its CN
 *
 *
 * @todo [core][interface] write python bindings
 * @todo [interface] write demo python program
 * @todo [core] factorize code, distribute to files
 * @todo [core] distinguish function and constraints by tags and get rid of redundant methods
 * @todo [RNA design] write sequence-sampling spefific convenience classes/methods
 * @todo [RNA design] write C++ function classes for parametrizable base-pair and stacking model, Turner model
 */

#include <iostream>

#include<limits>
#include <vector>
#include <memory>

#include <type_traits>


#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>


//namespace Boltzmann Sampling Constraint Framework
namespace bscf {

    // (partial) assignment assigns cn variables to values
    class Assignment {
    public:

        using assignment_t = Assignment;

        using var_idx_t = int;
        using var_value_t = int;
        static const int Undetermined = std::numeric_limits<int>::max();

        //! @brief construct 'empty' with size
        Assignment(int size)
            :values_(size, Undetermined)
        {
        }

        decltype(auto)
        operator [](var_idx_t i) const {return values_[i];}

        decltype(auto)
        operator [](var_idx_t i) {return values_[i];}

        const auto &
        values() const { return values_; }

        //! @brief is a variable determined by the assignment?
        auto
        is_det(var_idx_t i) {return values_[i] != Undetermined;}

        void
        reset(std::vector<var_idx_t> &vars) {
            for (auto var : vars) {
                values_[var] = var_value_t();
            }
        }


        // Iterate over possible assignments
        // assignment iterator works destructively on parent assignment!
        //
        // Special case: empty variable set: exactly one valid assignment
        template<class Constraint>
        class AssignmentIterator {
        public:
            using constraint_t = Constraint;

            /**
             * @see AssignmentIterator::reset()
             */
            AssignmentIterator(assignment_t &a,
                               const std::vector<var_idx_t> &vars,
                               const std::vector<int> &domsizes,
                               const std::vector<const Constraint *> &constraints)
                : a_(a),
                  vars_(vars),
                  domsizes_(domsizes),
                  constraint_board_(create_constraint_board(constraints))
            {
                reset();
            }

            // auto &
            // operator * () {
            //     return a_;
            // }


            /**
             * @brief Reset iterator
             *
             * Resets the assignment such that the iteration starts again
             *
             * After reset, the valuation is not necessarily valid (use ++ operator!)
             */
            auto
            reset() {
                // the entries for vars_ in a are (ab)used as a stack
                // top_ points to the next to be determined variable
                // var=vars_[top_]; a_[var] is the next candidate
                // value
                for (auto var : vars_) {
                    a_[var] = var_value_t();
                }
                top_ = 0;
                return *this;
            }

            //! check constraints for current top_ entry
            bool
            is_valid() const {
                //check all constraints at top_
                auto constraints = constraint_board_[top_];

                return all_of( constraints.begin(), constraints.end(),
                               [&](auto c) { return (*c)(a_); } );
            }

            /** @brief Next valid valuation
             *
             * Sets the assignemnt to the next valid valuation of
             * vars_
             */
            auto
            operator ++ () {

                if (top_ == vars_.size()) {
                    top_--;
                    if (vars_.size()==0) {
                        return *this;
                    }
                    a_[vars_[top_]]++;
                }

                while ( top_ < vars_.size() ) {
                    while( a_[vars_[top_]] < domsizes_[vars_[top_]] && !is_valid() ) {
                        a_[vars_[top_]] ++;
                    }
                    if ( a_[vars_[top_]] < domsizes_[vars_[top_]] ) { // i.e. the current value is valid!
                        top_++;
                    } else {
                        a_[vars_[top_]] = var_value_t();
                        top_--;
                        if (top_ < 0) break; // stack is empty -> termination
                        a_[vars_[top_]]++;
                    }
                }

                return *this;
            }


            /** @brief Check for termination
             */
            bool finished() {
                if ( vars_.size()==0 && !empty_case_flag_) {
                    empty_case_flag_=true;
                    return false;
                }

                return top_ < 0;
            }

        private:

            using constraint_board_t = std::vector<std::vector<const constraint_t *>>;

            assignment_t &a_;
            const std::vector<var_idx_t> &vars_;
            const std::vector<int> &domsizes_;
            constraint_board_t constraint_board_;
            int top_;
            bool empty_case_flag_ = false;


            constraint_board_t
            create_constraint_board(const std::vector<const constraint_t *> &constraints) {

                // the constraint board cb holds a collection cb[i] of to-be-checked constraints
                // for each variable i.
                constraint_board_t cb(vars_.size());

                for( decltype(auto) c : constraints ) {

                    // find 'last' variable of constraint c
                    int last_idx = -1;
                    for ( auto var : c->vars() ) {
                        // unless variable is assigned by a
                        if ( a_[var] != Undetermined ) continue;

                        // determine maximum index
                        auto idx = std::distance(vars_.begin(), find(vars_.begin(),vars_.end(),var));
                        last_idx = std::max(last_idx, (int)idx);
                    }

                    if ( 0 <= last_idx && last_idx < vars_.size() ) {
                        cb[last_idx].push_back( c );
                    }
                }

                return cb;
            }
        };

        /**
         * @see AssignmentIterator::reset()
         */
        template<class Constraint>
        auto
        make_iterator(const std::vector<var_idx_t> &vars,
                      const std::vector<int> &domsizes,
                      const std::vector<const Constraint *> &constraints) {
            return AssignmentIterator<Constraint>(*this, vars, domsizes, constraints);
        }

    private:
        std::vector<var_value_t> values_;
    };

    const int Assignment::Undetermined;

    // for debugging
    std::ostream &
    operator <<(std::ostream &out, const std::vector<int> &x) {
        std::copy(x.begin(),x.end(),std::ostream_iterator<int>(out,", "));
        return out;
    }


    /**
     * Dependencies specify a dependency between variables
     */
    class Dependency {
    public:
        using var_idx_t = int;

        Dependency(const std::vector<var_idx_t> &vars) : vars_(vars) {}

        const std::vector<var_idx_t> &
        vars() const {return vars_;}

    private:
        const std::vector<var_idx_t> vars_;
    };

    /**
     * Functions evaluate assignments of a subset of variables
     */
    template<class FunValue>
    class Function : public Dependency {
    public:
        using self_t = Function<FunValue>;
        using base_t = self_t;
        using assignment_t = Assignment;
        using fun_value_t = FunValue;

        Function(const std::vector<var_idx_t> &vars) : Dependency(vars) {}

        Function(const Function &fun) : Dependency(fun) {}

        virtual
        fun_value_t
        operator () (const assignment_t &) const {
            return fun_value_t();
        }

        virtual
        std::unique_ptr<base_t>
        clone() const {
            return std::make_unique<self_t>(*this);
        }

        virtual
        ~Function() {}
    };

    template<class FunValue, class UnderlyingFunction>
    class DynamicFunction : public Function<FunValue> {
    public:
        using self_t = DynamicFunction<FunValue,UnderlyingFunction>;
        using parent_t = Function<FunValue>;
        using base_t = typename parent_t::base_t;
        using var_idx_t = typename parent_t::var_idx_t;
        using assignment_t = typename parent_t::assignment_t;
        using fun_value_t = FunValue;


        template<class F>
        DynamicFunction(const std::vector<var_idx_t> &vars, const F &f) : parent_t(vars), f_(f) {}

        fun_value_t
        operator () (const assignment_t & a) const override { return f_(a); }

        std::unique_ptr<base_t>
        clone() const override {
            return std::make_unique<self_t>(*this);
        }

    private:
        using function_t = UnderlyingFunction;
        const function_t f_;
    };

    //! MaterializedFunction is used for the computed messages (DP
    //! matrices)
    //!
    //! @note via template parameter Container, one can
    //! choose to use a sparse data structure like
    //! Container=std::unordered_map<FunValue>
    //
    template<class FunValue, class Container=std::vector<FunValue> >
    class MaterializedFunction : public Function<FunValue> {
    public:
        using self_t = MaterializedFunction<FunValue,Container>;
        using parent_t = Function<FunValue>;
        using base_t = typename parent_t::base_t;

        using var_idx_t = typename parent_t::var_idx_t;
        using assignment_t = typename parent_t::assignment_t;
        using fun_value_t = FunValue;
        using data_t = Container;

        MaterializedFunction(const std::vector<var_idx_t> &vars, const std::vector<int> &domsizes) :
            parent_t(vars),
            domsizes_(domsizes),
            data_( calc_size() )
        {}

        fun_value_t
        operator () ( const assignment_t & a ) const override { return data_[ index_(a) ]; }

        std::unique_ptr<base_t>
        clone() const override {
            return std::make_unique<self_t>(*this);
        }

        void
        set( const assignment_t & a, const fun_value_t &val) {
            data_[ index_(a) ] = val;
        }

    private:

        data_t data_;

        const std::vector<int> domsizes_;

        int
        index_( const assignment_t & a ) const {
            int x = 0;

            for ( auto var : this->vars()) {
                x *= domsizes_[var];
                x += a[var];
            }

            return x;
        }

        auto
        calc_size() const {
            int x=1;

            for ( auto var : this->vars()) {
                x *= domsizes_[var];
            }
            return x;
        }
    };

    using Constraint = Function<bool>;

    template<class FunValue, class UnderlyingFunction>
    using DynamicConstraint = DynamicFunction<bool, UnderlyingFunction>;


    /**
     * Cluster in the cluster tree decomposition
     *
     * A cluster consists of collections
     *  - variables
     *  - constraints
     *  - functions
     */

    template<class FunValue=double>
    class Cluster {
    public:
        using var_idx_t = int;
        using function_t = Function<FunValue>;
        using constraint_t = Constraint;
        using assignment_t = Assignment;


        Cluster()
        {
        }

        Cluster(const std::vector<var_idx_t> &vars)
            : vars_(vars)
        {
        }

        Cluster(const std::vector<var_idx_t> &vars,
                const std::vector<constraint_t *> &constraints,
                const std::vector<function_t *> &functions
                )
            : vars_(vars),
              constrs_(constraints),
              funs_(functions)
        {
        }

        const std::vector<var_idx_t> &
        vars() const {return vars_;}

        const auto &
        constraints() const {
            return constrs_;
        }

        const auto &
        functions() const {return funs_;}

        void
        add_variable(const var_idx_t &x) {
            vars_.push_back(x);
        }

        void
        add_constraint(const constraint_t *c) {
            constrs_.push_back(c);
        }

        void
        add_function(const function_t *f) {
            funs_.push_back(f);
        }

        auto
        sep_vars(const Cluster &parent) {
            auto vars = std::vector<var_idx_t>();
            const auto &pvars = parent.vars();
            for ( auto &var : vars_ ) {
                if (std::find(pvars.begin(), pvars.end(),var)
                    != pvars.end()) {
                    vars.push_back(var);
                }
            }
            return vars;
        }

        auto
        diff_vars(const Cluster &parent) {
            auto vars = std::vector<var_idx_t>();
            const auto &pvars = parent.vars();
            for ( auto &var : vars_ ) {
                if (std::find(pvars.begin(), pvars.end(),var)
                    == pvars.end()) {
                    vars.push_back(var);
                }
            }
            return vars;
        }

    private:
        std::vector<var_idx_t> vars_;
        std::vector<const constraint_t *> constrs_;
        std::vector<const function_t *> funs_;
    };

    /** A constraint network consists of sets of
     *  - variables
     *  - constraints
     *  - functions
     *
     * @note The variables in a contraint network are indexed (0..n-1). Each variable has a domain of values 0..domsize
     */
    template<class FunValue=double>
    class ConstraintNetwork {
    public:
        using var_idx_t = int;
        using var_value_t = int;

        using fun_value_t = FunValue;

        using assignment_t = Assignment;
        using function_t = Function<FunValue>;
        using constraint_t = Constraint;

        //! @brief suitable cluster type for the constraint network
        using cluster_t = Cluster<fun_value_t>;

        /**
         * @brief Construct empty
         */
        ConstraintNetwork() {};

        /**
         * @brief Construct with domains
         */
        ConstraintNetwork(std::vector<int> &domsizes)
            : domsizes_(domsizes) {
        };

        /**
         * @brief Construct with uniform domains
         */
        ConstraintNetwork(int num_vars, int domsize)
            : domsizes_(num_vars,domsize) {
        };

        /**
         * @brief add variable with domain size
         *
         * @returns 0-based index of the new variable
         */
        var_idx_t
        add_var(int domainsize) {
            domsizes_.push_back(domainsize);
            return domsizes_.size()-1;
        }

        /**
         * @brief add constraint (moving ownership)
         */
        auto
        add_constraint(std::unique_ptr<constraint_t> &&x) {
            constraints_.push_back( std::move(x) );
            return constraints_.back().get();
        }

        /**
         * @brief add constraint (cloning)
         * @returns pointer to the new constraint
         */
        auto
        add_constraint(const constraint_t &x) {
            return add_constraint( x.clone() );
        }

        /**
         * @brief add function
         */
        auto
        add_function(std::unique_ptr<function_t> &&x) {
            functions_.push_back( std::move(x) );
            return functions_.back().get();
        }

        /**
         * @brief add function
         */
        auto
        add_function(const function_t &x)  {
            return add_function( x.clone() );
        }

        int
        num_vars() const {return (signed)domsizes_.size();}

        const auto & domsizes() const {
            return domsizes_;
        }

    private:
        std::vector<int> domsizes_;

        std::vector<std::unique_ptr<function_t>> functions_;
        std::vector<std::unique_ptr<constraint_t>> constraints_;
    };


    //NOTE: it should be possible to directly (manually) construct a tree decomposition
    // such that the cn is filled implicitely

    template<class FunValue>
    class StdEvaluationPolicy {
    public:
        using fun_value_t = FunValue;

        static
        fun_value_t
        plus(const fun_value_t &x, const fun_value_t &y) {
            return x+y;
        }

        static
        fun_value_t
        multiplies(const fun_value_t &x, const fun_value_t &y) {
            return x*y;
        }

        static
        fun_value_t
        one() {
            return fun_value_t(1);
        }

        static
        fun_value_t
        zero() {
            return fun_value_t();
        }

    };

    /** @brief A tree decomposition is a tree of Clusters
     */
    template<class FunValue=double, class EvaluationPolicy=StdEvaluationPolicy<FunValue>>
    class ClusterTreeDecomposition {

    public:
        using constraint_network_t = ConstraintNetwork<FunValue>;

        using cluster_t = typename constraint_network_t::cluster_t;
        using assignment_t = typename constraint_network_t::assignment_t;
        using fun_value_t = typename constraint_network_t::fun_value_t;
        using constraint_t = typename constraint_network_t::constraint_t;
        using function_t = typename constraint_network_t::function_t;


        using message_t = MaterializedFunction<fun_value_t>;


        struct vertex_info_t {
            cluster_t cluster;
        };

        struct edge_info_t {
            Function<fun_value_t> *message;
        };

        using tree_t = boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS, vertex_info_t, edge_info_t >;

        using vertex_descriptor = typename boost::graph_traits<tree_t>::vertex_descriptor;

        /**
         * @brief construct empty
         */
        ClusterTreeDecomposition() {
        }

        /**
         * @brief Construct with domains
         */
        ClusterTreeDecomposition(std::vector<int> &domsizes)
            : cn_(domsizes) {
        };

        /**
         * @brief Construct with uniform domains
         */
        ClusterTreeDecomposition(int num_vars, int domsize)
            : cn_(num_vars,domsize) {
        };


        //! @brief read access to constraint network
        const auto & constraint_network() const {
            return cn_;
        }

        //
        // enum class TDMethod {
        //     tdlib
        // };
        //
        //  /** @brief construct from constraint network by some tree decomposition method
        //  */
        // ClusterTreeDecomposition(constraint_network_t &n, TDMethod td_method): cn_(cn) {
        // }


        // for 'manual' construction
        auto
        make_root_cluster(const std::vector<int> &vars) {
            auto node = boost::add_vertex(tree_);
            tree_[node].cluster = cluster_t(vars);
            return node;
        }

        auto
        make_child_cluster( vertex_descriptor parent, const std::vector<int> &vars) {
            auto node = boost::add_vertex(tree_);
            boost::add_edge(tree_, parent, node);

            tree_[node].cluster = cluster_t(vars);
            return node;
        }

        //! @brief add constraint to cluster
        //!
        //! synchronized between cluster and cn
        auto
        add_constraint( vertex_descriptor node, const constraint_t &c ) {
            tree_[node].cluster.add_constraint( cn_.add_constraint(c) );
        }

        //! @brief add function to cluster
        //!
        //! synchronized between cluster and cn
        auto
        add_function( vertex_descriptor node, const function_t &c ) {
            tree_[node].cluster.add_function( cn_.add_function(c) );
        }

        // run the dynamic programming evaluation
        void
        evaluate();

        // generate one sample, using the DP tables
        //
        // this is a real 'sample' only if partition functions are
        // computed due to the evaluation policy!
        auto
        sample();

    private:
        constraint_network_t cn_;
        tree_t tree_;

        bool evaluated_ = false;
        bool single_rooted_ = false;

        vertex_descriptor root_;

        // insert pseudo root to connect trees of the forest (if necessary)
        // @returns (potentially new) root
        auto
        single_root();

        struct evaluate_finish_edge {
            using event_filter = boost::on_finish_edge;

            evaluate_finish_edge(constraint_network_t &cn, tree_t &tree)
                : cn_(cn), tree_(tree) {
            }

            template <class Edge, class Graph>
            void
            operator ()(Edge e, Graph &graph) {

                //std::cout << "evaluate_finish_edge: visit "<<e<<std::endl;

                auto parent = graph[ source(e, graph) ].cluster;
                auto child =  graph[ target(e, graph) ].cluster;


                // compute the message from cluster child to cluster parent

                auto sep  = child.sep_vars(parent);
                auto diff = child.diff_vars(parent);

                //std::cout << "Sep: " << sep << std::endl;
                //std::cout << "Diff: " << diff << std::endl;

                auto message = message_t(sep,cn_.domsizes());

                auto a = Assignment(cn_.num_vars());

                auto it = a.make_iterator(sep, cn_.domsizes(), child.constraints());

                for( ++it ; ! it.finished() ; ++it ) {
                    fun_value_t x = EvaluationPolicy::zero();

                    auto it2 = a.make_iterator(diff, cn_.domsizes(), child.constraints());

                    for( ++it2; ! it2.finished(); ++it2) {

                        fun_value_t p = EvaluationPolicy::one();
                        for ( auto f : child.functions() ) {
                            p = EvaluationPolicy::multiplies(p, (*f)(a) );
                        }
                        x = EvaluationPolicy::plus(x, p);
                    }

                    message.set(a, x);

                    a.reset(diff);
                }

                // register message in cn, such that it persists!
                auto msg = cn_.add_function(message);
                // then, register in cluster parent
                parent.add_function(msg);
                // ... and as edge property
                tree_[e].message = msg;
            }

        private:
            constraint_network_t &cn_;
            tree_t &tree_;
        };

        struct sample_examine_edge {
            using event_filter = boost::on_examine_edge;

            sample_examine_edge(constraint_network_t &cn,assignment_t &a) : cn_(cn), a_(a) {}

            template <class Edge, class Graph>
            void
            operator ()(Edge e, Graph &graph) {
                //std::cout << "sample_examine_edge: visit "<<e<<std::endl;

                auto parent = graph[ source(e, graph) ].cluster;
                auto child =  graph[ target(e, graph) ].cluster;

                auto diff = child.diff_vars(parent);

                auto message = graph[e].message;

                double r = rand()/(RAND_MAX+1.0) * (*message)(a_);
                
                //std::cout << r << " " << (*message)(a_) << std::endl;

                auto x = EvaluationPolicy::zero();
                for( auto it = ++ a_.make_iterator(diff, cn_.domsizes(), child.constraints()); ! it.finished(); ++it ) {
                    fun_value_t p = EvaluationPolicy::one();
                    for ( auto f : child.functions() ) {
                        p = EvaluationPolicy::multiplies(p, (*f)(a_) );
                    }

                    x = EvaluationPolicy::plus(x, p);

                    if ( x > r ) {
                        break;
                    }
                }
            }
        private:
            constraint_network_t &cn_;
            assignment_t &a_;
        };

    };

    template<class ConstraintNetwork, class EvaluationPolicy>
    auto
    ClusterTreeDecomposition<ConstraintNetwork, EvaluationPolicy>::single_root() {
        if (single_rooted_) {return root_;}

        // find all root nodes of the tree
        auto index = boost::get(boost::vertex_index, tree_);

        std::set<int> no_roots;
        std::vector<vertex_descriptor> old_roots;

        for (auto it = edges(tree_).first; it != edges(tree_).second; ++it) {
            no_roots.insert(index[target(*it,tree_)]);
        }

        for (auto it = vertices(tree_).first; it != vertices(tree_).second; ++it) {
            if ( no_roots.find(index[*it]) == no_roots.end() ) {
                old_roots.push_back(*it);
            }
        }

        // insert new root and point to all old roots
        auto new_root = boost::add_vertex(tree_);

        for (auto old_root : old_roots) {
            boost::add_edge(new_root, old_root , tree_);
        }

        single_rooted_ = true;
        root_ = new_root;
        return root_;
    }

    template<class ConstraintNetwork, class EvaluationPolicy>
    void
    ClusterTreeDecomposition<ConstraintNetwork,EvaluationPolicy>
    ::evaluate() {

        auto root = single_root();

        auto cte_visitor = boost::make_dfs_visitor(evaluate_finish_edge(cn_,tree_));

        boost::depth_first_search(tree_, visitor(cte_visitor).root_vertex(root));

        evaluated_ = true;
    }


    template<class ConstraintNetwork, class EvaluationPolicy>
    auto
    ClusterTreeDecomposition<ConstraintNetwork,EvaluationPolicy>::sample() {

        if (!evaluated_) {
            evaluate();
        }

        auto a = assignment_t(cn_.num_vars());

        auto sample_visitor = boost::make_dfs_visitor(sample_examine_edge(cn_,a));

        boost::depth_first_search(tree_, visitor(sample_visitor).root_vertex(single_root()));

        return a;
    }

}

#endif
