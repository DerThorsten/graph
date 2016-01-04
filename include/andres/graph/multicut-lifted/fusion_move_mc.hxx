#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_LIFTED_FUSION_MC_HXX
#define ANDRES_GRAPH_MULTICUT_LIFTED_FUSION_MC_HXX

#include <vector>
#include <random>
#include <iostream> 

#include "kernighan-lin.hxx"
#include "andres/graph/runtime_check.hxx"
#include "andres/graph/components.hxx"
#include "fusion_mover.hxx"
#include "greedy-additive.hxx"

namespace andres {
namespace graph {
namespace multicut_lifted {


struct FusionMoveMcSettings {
    std::size_t maxNumberOfIterations { 1000 };
    std::size_t maxNumberOfIterationsWithoutImprovement { 200 };
    double nodeLimit = {0.1};
    double sigma{ 1.0};
    int verbose { 1 };
};





template<typename ORIGINAL_GRAPH, typename LIFTED_GRAPH, typename ECA, typename ELA, typename VIS>
inline
void fusionMoveMc(
    const ORIGINAL_GRAPH& originalGraph,
    const LIFTED_GRAPH& liftedGraph,
    const ECA& edgeCosts,
    const ELA& inputLabels,
    ELA& outputLabels,
    VIS& visitor,
    const FusionMoveMcSettings & settings
){

    FusionMover<ORIGINAL_GRAPH, LIFTED_GRAPH, ECA> fusionMover(originalGraph, liftedGraph, edgeCosts);

    // lambda to compute energy 
    auto getEnergy = [&] (const std::vector<char> & edgeLabels_) {
        auto totalCost = 0.0;
        for(std::size_t edge = 0; edge < liftedGraph.numberOfEdges(); ++edge){
            if(edgeLabels_[edge]){
                totalCost += edgeCosts[edge];
            }
        }
        return totalCost;
    };
    
    std::vector<char> lBest(inputLabels.begin(),inputLabels.end()),
                      lNoisy(inputLabels.size()),
                      lFusionRes(inputLabels.size()),
                      lOnes(inputLabels.size(),1);
   


    auto cNonImprovement = 0; 
    auto anyValid = false;
    std::mt19937 gen(42); // TODO externalize this? use seed?
    auto dist = std::normal_distribution<>(0.0,settings.sigma);
    auto rg = [&](){return dist(gen);}; 
    for(auto i=0; i< settings.maxNumberOfIterations; ++i){


        // ( 1 ) build a noisy objective
        //std::cout<<"     build noisy obj\n";
        auto edgeValuesNoisy = std::vector<double>(edgeCosts.begin(),edgeCosts.end()); // not type save TODO fixeme
        for(auto & ec : edgeValuesNoisy)
            ec += rg();

        // ( 2 ) use greedy additive alg. on noise objective, and
        //       contract edges until only N^ are left and
        //       save labeling  as lNoisy
        //std::cout<<"     do  greedy additive obj\n";
        int stopCond = -1;
        if(settings.nodeLimit<1.0){
            stopCond = double(originalGraph.numberOfVertices())*settings.nodeLimit + 0.5;
        }
        else{
            stopCond = settings.nodeLimit + 0.5;
        }
        //std::count<<"stop cond\n";
        stopCond = std::max(2,stopCond);
        greedyAdditiveEdgeContraction(originalGraph,liftedGraph, edgeValuesNoisy, lNoisy, stopCond);
        //KernighanLinSettings s;
        //s.verbose = false;
        //kernighanLin(originalGraph,liftedGraph,edgeValuesNoisy,lOnes, lNoisy,s);

        // check that proposal has at least one cut
        auto anyCut = false;
        for(std::size_t edge = 0; edge < liftedGraph.numberOfEdges(); ++edge){
            auto isCut = (lBest[edge]==1 || lNoisy[edge]==1 ) ? 1 : 0;
            if(isCut)
                anyCut = true;
        }
        if(!anyCut)
            continue;
        else{
            // accept the first proposal with cut (!invalidate starting point)
            if(!anyValid){
                lBest = lNoisy;
                anyValid = true;
                continue;
            }
        }

        // ( 3  ) fuse lBest with lNoisy: 
        // ( 3.1) by setting up "fused objective" O_FS as in beier_2015 
        //
        //      working in the node label domain would be cheaper
        const auto ePre  = getEnergy(lBest);
        fusionMover.fuse(lBest, lNoisy, lFusionRes);
        std::copy(lFusionRes.begin(), lFusionRes.end(), lBest.begin());
        const auto ePost = getEnergy(lBest);
        ++cNonImprovement;
        if(ePost<ePre)
            cNonImprovement = 0;
        if(settings.verbose>= 1)
            std::cout<<"iter "<<i<<" Energy "<<ePost<< " since "<<cNonImprovement<<"\n";

        if(cNonImprovement >= settings.maxNumberOfIterationsWithoutImprovement){
            break;
        }
    }
    std::copy(lBest.begin(), lBest.end(), outputLabels.begin());

}

template<typename ORIGINAL_GRAPH, typename LIFTED_GRAPH, typename ECA, typename ELA>
void fusionMoveMc(
    const ORIGINAL_GRAPH& originalGraph,
    const LIFTED_GRAPH& liftedGraph,
    const ECA& edgeCosts,
    const ELA& inputLabels,
    ELA& outputLabels,
    const FusionMoveMcSettings settings = FusionMoveMcSettings())
{
    struct Visitor {
        constexpr bool operator()(std::vector<std::size_t>& vertex_labels)
            { return true; }
        constexpr bool time_limit_exceeded()
            { return false; }
    } visitor;

    fusionMoveMc(originalGraph, liftedGraph, edgeCosts, inputLabels, outputLabels, visitor, settings);
}

}
}
}

#endif
