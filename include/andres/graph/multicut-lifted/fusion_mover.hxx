#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_LIFTED_FUSION_MOVER_HXX
#define ANDRES_GRAPH_MULTICUT_LIFTED_FUSION_MOVER_HXX

#include <vector>
#include <random>
#include <iostream> 

#include "kernighan-lin.hxx"
#include "andres/graph/runtime_check.hxx"
#include "andres/graph/components.hxx"
#include "greedy-additive.hxx"

namespace andres {
namespace graph {
namespace multicut_lifted {





template<class ORIGINAL_GRAPH, class LIFTED_GRAPH, class E_COSTS>
class FusionMover{
public:
    FusionMover(const ORIGINAL_GRAPH & originalGraph, const LIFTED_GRAPH & liftedGraph, const E_COSTS & edgeCosts)
    :   originalGraph_(originalGraph),
        liftedGraph_(liftedGraph),
        edgeCosts_(edgeCosts),
        lAB_(liftedGraph.numberOfEdges()),
        edgeInLiftedGraph_(originalGraph.numberOfEdges()){

        for (std::size_t i = 0; i < originalGraph.numberOfEdges(); ++i)
        {
            auto v0 = originalGraph.vertexOfEdge(i, 0);
            auto v1 = originalGraph.vertexOfEdge(i, 1);
            edgeInLiftedGraph_[i] = liftedGraph.findEdge(v0, v1).second;
        }  
    }

    template<class L>
    double getEnergy(const L & l){
        auto totalCost = 0.0;
        for(std::size_t edge = 0; edge < liftedGraph_.numberOfEdges(); ++edge){
            if(l[edge]){
                totalCost += edgeCosts_[edge];
            }
        }
        return totalCost;
    }


    template<class L_A, class L_B, class L_RESULT>
    bool fuse(
        const L_A & lA, 
        const L_B & lB, 
        L_RESULT & lResult
    ){
        struct SubgraphWithCut { // a subgraph with cut mask
            SubgraphWithCut(const std::vector<char> & labels, std::vector<std::size_t> const& edge_in_lifted_graph)
                : labels_(labels), edge_in_lifted_graph_(edge_in_lifted_graph)
                {}
            bool vertex(const std::size_t v) const
                { return true; }
            bool edge(const std::size_t e) const
                { return labels_[edge_in_lifted_graph_[e]] == 0; }

            std::vector<std::size_t> const& edge_in_lifted_graph_;
            const std::vector<char> & labels_;
        };

        const double eA = getEnergy(lA);
        const double eB = getEnergy(lB);
        //std::cout<<"eA "<<eA<<"\n";
        //std::cout<<"eB "<<eB<<"\n";
        //const double eA = getEnergy(lA);


        auto anyCut = false;
        for(std::size_t edge = 0; edge < liftedGraph_.numberOfEdges(); ++edge){
            auto isCut = (lA[edge]==1 || lB[edge]==1 ) ? 1 : 0;
            lAB_[edge] = isCut;
            if(isCut)
                anyCut = true;
        }
        if(!anyCut){
            //std::cout << "not cut!\n";
            std::copy(lA.begin(), lA.end(), lResult.begin());
            return false;
        }









        ComponentsBySearch<ORIGINAL_GRAPH> components;
        components.build(originalGraph_, SubgraphWithCut(lAB_, edgeInLiftedGraph_));

        // setup lifted mc problem for the subgraph
        
        const auto numberOfSubgraphNodes = *std::max_element(components.labels_.begin(), components.labels_.end()) +1;
        //std::cout<<"    subgraph nodes "<<numberOfSubgraphNodes<<"\n";

        Graph<> subgraphOriginal(numberOfSubgraphNodes);
        Graph<> subgraphLifted(numberOfSubgraphNodes);

        // fill subgraphOriginal
        for(std::size_t edge = 0; edge < originalGraph_.numberOfEdges(); ++edge){
            const auto v0 = originalGraph_.vertexOfEdge(edge, 0);
            const auto v1 = originalGraph_.vertexOfEdge(edge, 1);
            const auto l0 = components.labels_[v0];
            const auto l1 = components.labels_[v1];
            if(l0 != l1){
                GRAPH_TEST_OP(l0,<,subgraphOriginal.numberOfVertices())
                GRAPH_TEST_OP(l1,<,subgraphOriginal.numberOfVertices())
                if(!subgraphOriginal.findEdge(l0, l1).first){
                    subgraphOriginal.insertEdge(l0, l1);
                }
            }
        }
        //std::cout<<"    subgraph original edges "<<subgraphOriginal.numberOfEdges()<<"\n";

        // fill lifted graph
        for(std::size_t edge = 0; edge < liftedGraph_.numberOfEdges(); ++edge){

            const auto v0 = liftedGraph_.vertexOfEdge(edge, 0);
            const auto v1 = liftedGraph_.vertexOfEdge(edge, 1);
            const auto l0 = components.labels_[v0];
            const auto l1 = components.labels_[v1];
            if(l0 != l1){
                GRAPH_TEST_OP(l0,<,subgraphLifted.numberOfVertices())
                GRAPH_TEST_OP(l1,<,subgraphLifted.numberOfVertices())
                if(!subgraphLifted.findEdge(l0, l1).first){
                    subgraphLifted.insertEdge(l0, l1);
                }
            }
        }
        //std::cout<<"    subgraph lifted edges "<<subgraphLifted.numberOfEdges()<<"\n";
        std::vector<double> subgraphEdgeCosts(subgraphLifted.numberOfEdges(),0.0);

        // fill lifted graph weights
        //std::cout<<"    lifted cost setup "<<numberOfSubgraphNodes<<"\n";
        for(std::size_t edge = 0; edge < liftedGraph_.numberOfEdges(); ++edge){

            const auto v0 = liftedGraph_.vertexOfEdge(edge, 0);
            const auto v1 = liftedGraph_.vertexOfEdge(edge, 1);
            const auto l0 = components.labels_[v0];
            const auto l1 = components.labels_[v1];
            if(l0 != l1){
                const auto findEdge = subgraphLifted.findEdge(l0,l1);
                GRAPH_TEST(findEdge.first);
                GRAPH_TEST_OP(findEdge.second,<,subgraphEdgeCosts.size())
                subgraphEdgeCosts[findEdge.second] += edgeCosts_[edge];
            }
        }
        //std::cout<<"    lifted cost setup done "<<numberOfSubgraphNodes<<"\n";

        std::vector<int> subgraphInput(subgraphLifted.numberOfEdges(),1);
        std::vector<int> subgraphRes(subgraphLifted.numberOfEdges());
        KernighanLinSettings s;
        s.verbose = false;
        kernighanLin(subgraphOriginal,subgraphLifted,subgraphEdgeCosts,subgraphInput, subgraphRes,s);
        //std::cout<<"    inference done "<<numberOfSubgraphNodes<<"\n";
        
        for(std::size_t edge = 0; edge < liftedGraph_.numberOfEdges(); ++edge){
            const auto v0 = liftedGraph_.vertexOfEdge(edge, 0);
            const auto v1 = liftedGraph_.vertexOfEdge(edge, 1);
            const auto l0 = components.labels_[v0];
            const auto l1 = components.labels_[v1];
            if(l0 != l1){
                const auto subgraphLiftedFindEdge = subgraphLifted.findEdge(l0,l1);
                GRAPH_TEST(subgraphLiftedFindEdge.first);
                lAB_[edge] = subgraphRes[subgraphLiftedFindEdge.second];
            }
            else{
                lAB_[edge] = 0;
            }
        }
        
        const double eFS = getEnergy(lAB_);
        if(eFS<std::min(eA,eB)){
            std::copy(lAB_.begin(), lAB_.end(), lResult.begin());
            return true;
        }
        else if(eA<eB){
            std::copy(lA.begin(), lA.end(), lResult.begin());
            return false;
        }
        else{
            std::copy(lB.begin(), lB.end(), lResult.begin());
            return true;
        }
    }
    
private:
    const ORIGINAL_GRAPH & originalGraph_;
    const LIFTED_GRAPH   & liftedGraph_;
    const E_COSTS & edgeCosts_;
    std::vector<char> lAB_;
    std::vector<std::size_t> edgeInLiftedGraph_;
    
};


}
}
}

#endif
