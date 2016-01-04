#pragma once
#ifndef ANDRES_GRAPH_LIFTED_GRAPH_HXX
#define ANDRES_GRAPH_LIFTED_GRAPH_HXX

#include <vector>
#include <random>
#include <iostream> 

#include "kernighan-lin.hxx"
#include "andres/graph/runtime_check.hxx"
#include "andres/graph/components.hxx"
#include "greedy-additive.hxx"

namespace andres {
namespace graph {

    template<class ORGINAL_GRAPH>
    class LiftedGraph{
    public:
        class LiftedGraphAdditionalEdge{
        };

    private:
        const ORGINAL_GRAPH & orginalGraph_;
        
    };
}
}

#endif ANDRES_GRAPH_LIFTED_GRAPH_HXX
