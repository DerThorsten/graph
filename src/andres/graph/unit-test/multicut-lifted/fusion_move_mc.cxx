#include <iostream>
#include <random>
#include "andres/graph/runtime_check.hxx"
#include "andres/graph/graph.hxx"
#include "andres/graph/complete-graph.hxx"
#include "andres/graph/multicut-lifted/fusion_move_mc.hxx"


using namespace andres::graph;

void testMulticutLifted()
{
    Graph<> original_graph(5);
    original_graph.insertEdge(0, 1); // 0
    original_graph.insertEdge(0, 3); // 1
    original_graph.insertEdge(1, 2); // 2
    original_graph.insertEdge(1, 4); // 3
    original_graph.insertEdge(3, 4); // 4

    CompleteGraph<> lifted_graph(5);
    
    std::vector<double> weights(10);
    weights[lifted_graph.findEdge(0, 1).second] = 10;
    weights[lifted_graph.findEdge(0, 2).second] = -1;
    weights[lifted_graph.findEdge(0, 3).second] = -1;
    weights[lifted_graph.findEdge(0, 4).second] = -1;
    weights[lifted_graph.findEdge(1, 2).second] = 10;
    weights[lifted_graph.findEdge(1, 3).second] = -1;
    weights[lifted_graph.findEdge(1, 4).second] = 4;
    weights[lifted_graph.findEdge(2, 3).second] = -1;
    weights[lifted_graph.findEdge(2, 4).second] = -1;
    weights[lifted_graph.findEdge(3, 4).second] = 10;

    std::vector<int> edge_labels(lifted_graph.numberOfEdges(), 0);
    multicut_lifted::fusionMoveMc(original_graph, lifted_graph, weights, edge_labels, edge_labels);

    GRAPH_TEST_OP(int(edge_labels[lifted_graph.findEdge(0, 1).second]), ==, 0);
    GRAPH_TEST_OP(int(edge_labels[lifted_graph.findEdge(0, 2).second]), ==, 0);
    GRAPH_TEST_OP(int(edge_labels[lifted_graph.findEdge(0, 3).second]), ==, 1);
    GRAPH_TEST_OP(int(edge_labels[lifted_graph.findEdge(0, 4).second]), ==, 1);
    GRAPH_TEST_OP(int(edge_labels[lifted_graph.findEdge(1, 2).second]), ==, 0);
    GRAPH_TEST_OP(int(edge_labels[lifted_graph.findEdge(1, 3).second]), ==, 1);
    GRAPH_TEST_OP(int(edge_labels[lifted_graph.findEdge(1, 4).second]), ==, 1);
    GRAPH_TEST_OP(int(edge_labels[lifted_graph.findEdge(2, 3).second]), ==, 1);
    GRAPH_TEST_OP(int(edge_labels[lifted_graph.findEdge(2, 4).second]), ==, 1);
    GRAPH_TEST_OP(int(edge_labels[lifted_graph.findEdge(3, 4).second]), ==, 0);
}


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

void testMulticutLifted2()
{



    auto sX = 20;
    auto sY = 20;
    double beta = 20.0;
    std::random_device rd;
    std::mt19937 gen(42); // TODO externalize this? use seed?
    auto dist = std::normal_distribution<>(0.0,1.0);
    auto rg = [&](){return 4.0;}; 


    Graph<> original_graph(sX*sY);
    Graph<> lifted_graph(sX*sY);

    auto getVi = [&](int x, int y){
        return x +y*sX;
    };

    for(auto y=0; y<sY; ++y)
    for(auto x=0; x<sX; ++x){
        if(x+1<sX){
            original_graph.insertEdge(getVi(x,y),getVi(x+1,y)); 
            lifted_graph.insertEdge(getVi(x,y),getVi(x+1,y)); 
        }
        if(y+1<sY){
            original_graph.insertEdge(getVi(x,y),getVi(x,y+1)); 
            lifted_graph.insertEdge(getVi(x,y),getVi(x,y+1)); 
        }
    }

    for(auto y=0; y<sY; ++y){
        // edge from very left to very right
        // very repuslive
        lifted_graph.insertEdge(getVi(0,y),getVi(sX-1,y)); 
    }

    for(auto x=0; x<sX; ++x){
        // edge from top to below 
        // very attractive
        lifted_graph.insertEdge(getVi(x,0),getVi(x,sY-1)); 
    }


    std::vector<double> weights(lifted_graph.numberOfEdges());

    //fill weights
    for(auto y=0; y<sY; ++y)
    for(auto x=0; x<sX; ++x){
        if(x+1<sX){
            auto e = lifted_graph.findEdge(getVi(x,y),getVi(x+1,y)).second;
            weights[e] = rg();
        }
        if(y+1<sY){
            auto e = lifted_graph.findEdge(getVi(x,y),getVi(x,y+1)).second;
            weights[e] = rg();
        }
    }

    for(auto y=0; y<sY; ++y){
        // edge from very left to very right
        // very repuslive
        auto e  = lifted_graph.findEdge(getVi(0,y),getVi(sX-1,y)).second; 
        weights[e] = -1.0*beta;
    }

    for(auto x=0; x<sX; ++x){
        // edge from top to below 
        // very attractive
        auto e = lifted_graph.findEdge(getVi(x,0),getVi(x,sY-1)).second; 
        weights[e] = beta;
    }

    // in the top middle and bottom middle great repulsiveness
    
    auto tm = lifted_graph.findEdge(getVi(sX/2-1,0),getVi(sX/2,0)).second; 
    auto bm = lifted_graph.findEdge(getVi(sX/2-1,sY-1),getVi(sX/2,sY-1)).second; 
    weights[tm] = -1.0 * beta;
    weights[bm] = -1.0 * beta;



    // lambda to compute energy 
    auto getEnergy = [&] (const std::vector<char>  & edgeLabels_) {
        auto totalCost = 0.0;
        for(std::size_t edge = 0; edge < lifted_graph.numberOfEdges(); ++edge){
            if(edgeLabels_[edge]){
                totalCost += weights[edge];
            }
        }
        return totalCost;
    };



    std::vector<char> edge_labels(lifted_graph.numberOfEdges(), 0);
    multicut_lifted::fusionMoveMc(original_graph, lifted_graph, weights, edge_labels, edge_labels);
    GRAPH_TEST_OP(getEnergy(edge_labels),==,-368);


    std::vector<size_t> edgeInLiftedGraph(original_graph.numberOfEdges());
    for (std::size_t i = 0; i < original_graph.numberOfEdges(); ++i)
    {
        auto v0 = original_graph.vertexOfEdge(i, 0);
        auto v1 = original_graph.vertexOfEdge(i, 1);
        edgeInLiftedGraph[i] = lifted_graph.findEdge(v0, v1).second;
    }  

    ComponentsBySearch<Graph<> > components;
    components.build(original_graph, SubgraphWithCut(edge_labels, edgeInLiftedGraph));

    for(auto y=0; y<sY; ++y){
        for(auto x=0; x<sX; ++x){
             std::cout.width(4);
            std::cout<<components.labels_[getVi(x,y)]<<" ";
        }
        std::cout<<"\n";
    }

    //GRAPH_TEST_OP(int(edge_labels[lifted_graph.findEdge(0, 1).second]), ==, 0);
   // GRAPH_TEST_OP(int(edge_labels[lifted_graph.findEdge(0, 2).second]), ==, 0);
   // GRAPH_TEST_OP(int(edge_labels[lifted_graph.findEdge(0, 3).second]), ==, 1);
   // GRAPH_TEST_OP(int(edge_labels[lifted_graph.findEdge(0, 4).second]), ==, 1);
   // GRAPH_TEST_OP(int(edge_labels[lifted_graph.findEdge(1, 2).second]), ==, 0);
   // GRAPH_TEST_OP(int(edge_labels[lifted_graph.findEdge(1, 3).second]), ==, 1);
   // GRAPH_TEST_OP(int(edge_labels[lifted_graph.findEdge(1, 4).second]), ==, 1);
   // GRAPH_TEST_OP(int(edge_labels[lifted_graph.findEdge(2, 3).second]), ==, 1);
   // GRAPH_TEST_OP(int(edge_labels[lifted_graph.findEdge(2, 4).second]), ==, 1);
   // GRAPH_TEST_OP(int(edge_labels[lifted_graph.findEdge(3, 4).second]), ==, 0);
}

int main()
{
    testMulticutLifted();
    testMulticutLifted2();
    return 0;
}
