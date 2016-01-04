// this must define the same symbol as the main module file (numpy requirement)
#define PY_ARRAY_UNIQUE_SYMBOL andres_graph_PyArray_API
#define NO_IMPORT_ARRAY

#include <Python.h>
#include <boost/python.hpp>
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>

#include "andres/graph/grid-graph.hxx"
#include "andres/graph/graph.hxx"
#include "andres/graph/multicut-lifted/fusion_move_mc.hxx"

namespace bp = boost::python;
namespace agraph = andres::graph;



// struct FusionMoveMcSettings {
//     std::size_t maxNumberOfIterations { 1000 };
//     std::size_t maxNumberOfIterationsWithoutImprovement { 200 };
//     double nodeLimit = {0.1};
//     double sigma{ 1.0};
//     int verbose { 1 };
// };



template<class OG, class LG>
vigra::NumpyAnyArray pyLiftedMcFusionMoves(
    const OG & og,
    const LG & lg,
    vigra::NumpyArray<1, float>   weights,
    vigra::NumpyArray<1, uint8_t> startingPoint,
    vigra::NumpyArray<1, uint8_t> result,
    const std::size_t maxNumberOfIterations,
    const std::size_t maxNumberOfIterationsWithoutImprovement,
    const double nodeLimit,
    const int verbose
){
    
    vigra::TinyVector<int, 1> shape(lg.numberOfEdges());
    agraph::multicut_lifted::FusionMoveMcSettings s;
    s.maxNumberOfIterations = maxNumberOfIterations;
    s.maxNumberOfIterationsWithoutImprovement = maxNumberOfIterationsWithoutImprovement;
    s.nodeLimit = nodeLimit;
    s.verbose = verbose;

    startingPoint.reshapeIfEmpty(shape);
    result.reshapeIfEmpty(shape);
    {;
        auto weights_ = vigra::MultiArrayView<1, float>(weights);
        auto startingPoint_ = vigra::MultiArrayView<1, uint8_t>(startingPoint);
        auto result_ = vigra::MultiArrayView<1, uint8_t>(result);
        agraph::multicut_lifted::fusionMoveMc(og, lg, weights_, startingPoint_, result_, s);
    }
    return result;
}


template<class OG, class LG>
vigra::NumpyAnyArray pyImageToLocalWeights(
    const OG & og,
    const LG & lg,
    vigra::NumpyArray<2, float>   image,
    vigra::NumpyArray<1, float> out
){
    
    vigra::TinyVector<int, 1> shape(lg.numberOfEdges());
    out.reshapeIfEmpty(shape);
    for(auto e=0; e<og.numberOfEdges(); ++e){
        auto u = og.vertexOfEdge(e, 0);
        auto v = og.vertexOfEdge(e, 1);
        out[e] = (image[u] + image[v])/2.0;
    }
    return out;
}


template<class OG, class LG>
vigra::NumpyAnyArray pyLiftedEdgeLabesToImage(
    const OG & og,
    const LG & lg,
    vigra::NumpyArray<1, uint8_t> edgeLabels,
    vigra::NumpyArray<2, uint64_t>   image
){
    

    struct SubgraphWithCut { // a subgraph with cut mask
        SubgraphWithCut(const vigra::MultiArrayView<1,uint8_t> & labels, 
                        std::vector<std::size_t> const& edge_in_lifted_graph)
            : labels_(labels), edge_in_lifted_graph_(edge_in_lifted_graph)
            {}
        bool vertex(const std::size_t v) const
            { return true; }
        bool edge(const std::size_t e) const
            { return labels_[edge_in_lifted_graph_[e]] == 0; }

        std::vector<std::size_t> const& edge_in_lifted_graph_;
        const vigra::MultiArrayView<1,uint8_t> & labels_;
    };

    vigra::TinyVector<int, 2> shape(og.shape(0),og.shape(1));
    image.reshapeIfEmpty(shape);
    
    std::vector<size_t> edgeInLiftedGraph(og.numberOfEdges());
    for (std::size_t i = 0; i < og.numberOfEdges(); ++i){
        auto v0 = og.vertexOfEdge(i, 0);
        auto v1 = og.vertexOfEdge(i, 1);
        edgeInLiftedGraph[i] = lg.findEdge(v0, v1).second;
    }  
    auto edgeLabels_ = vigra::MultiArrayView<1,uint8_t>(edgeLabels);
    agraph::ComponentsBySearch<OG > components;
    components.build(og, SubgraphWithCut(edgeLabels_, edgeInLiftedGraph));

    for(size_t i=0; i<image.size(); ++i){
        image[i] =components.labels_[i];
    }
    return image;
}



template<class OG, class LG>
void exportLiftedMcFusionMoves(){
    bp::def("liftedMcFusionMoves",
        vigra::registerConverters(&pyLiftedMcFusionMoves<OG,LG>),
        (
            bp::arg("graph"),
            bp::arg("liftedGraph"),
            bp::arg("weights"),
            bp::arg("startingPoint") = bp::object(),
            bp::arg("out") = bp::object(),
            bp::arg("maxNumberOfIterations") = 1000,
            bp::arg("maxNumberOfIterationsWithoutImprovement") = 200,
            bp::arg("nodeLimit") = 0.1,
            bp::arg("verbose") = 1
        )
    )
    ;
}


void exportLiftedMc(){
    typedef agraph::GridGraph<2> GridGraph2D;
    typedef agraph::GridGraph<3> GridGraph3D;
    typedef agraph::Graph<> Graph;

    exportLiftedMcFusionMoves<GridGraph2D, Graph>();
    exportLiftedMcFusionMoves<GridGraph3D, Graph>();
    exportLiftedMcFusionMoves<Graph, Graph>();


    bp::def("imageToLocalWeights",vigra::registerConverters(&pyImageToLocalWeights<GridGraph2D,Graph>),
        (
            bp::arg("graph"),
            bp::arg("liftedGraph"),
            bp::arg("image"),
            bp::arg("out") = bp::object()
        )
    );

    bp::def("liftedEdgeLabesToImage",vigra::registerConverters(&pyLiftedEdgeLabesToImage<GridGraph2D,Graph>),
        (
            bp::arg("graph"),
            bp::arg("liftedGraph"),
            bp::arg("edgeLabels"),
            bp::arg("out") = bp::object()
        )
    );
}
