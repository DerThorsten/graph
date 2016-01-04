// this must define the same symbol as the main module file (numpy requirement)
#define PY_ARRAY_UNIQUE_SYMBOL andres_graph_PyArray_API
#define NO_IMPORT_ARRAY

#include <Python.h>
#include <boost/python.hpp>
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>

#include <vector>


#include "andres/graph/grid-graph.hxx"
#include "andres/graph/graph.hxx"
#include "andres/graph/threadpool.hxx"
#include "andres/graph/multicut-lifted/lifted_mc_model.hxx"
namespace bp = boost::python;
namespace agraph = andres::graph;







template<class MODEL>
void addLongRangeEdges(
    MODEL & model,
    vigra::NumpyArray<2, float> edgePmap, // in 0-1
    const float beta = 0.5f,
    const int minRadius = 2,
    const int maxRadius = 7                                                
){  
    GRAPH_CHECK_OP(model.originalGraph().shape(0), == ,edgePmap.shape(0),"");
    GRAPH_CHECK_OP(model.originalGraph().shape(1), == ,edgePmap.shape(1),"");

    typedef vigra::TinyVector<int,   2> Coord;
    const auto shape = edgePmap.shape();

    auto clipToImg = [&](const Coord & coord){
        auto c = coord;
        for(auto i=0; i<2; ++i){
            c[i] = std::max(0, c[i]);
            c[i] = std::min(int(shape[i]),c[i]);
        }
        return c;
    };


    // move a coordinate to local min
    auto moveToMin = [&](const Coord & coord){
        Coord coords[5] = {
                coord,
                clipToImg(coord+Coord(0,1)),
                clipToImg(coord+Coord(1,0)),
                clipToImg(coord+Coord(-1,-1)),
                clipToImg(coord+Coord(-1, 0))
        };
        auto minVal = std::numeric_limits<float>::infinity();
        auto minCoord = Coord();
        for(size_t i=0; i<5; ++i){
            const auto val = edgePmap[coord];
            if(val<minVal){
                minVal = val;
                minCoord = coord;
            }
        }
        return minCoord;
    };
    
    auto & originalGraph = model.originalGraph();
    auto & liftedGraph = model.liftedGraph();






    int rad = 5;
    
    auto pOpt = agraph::ParallelOptions();
    const auto nThreads = pOpt.numThreads(-1).getActualNumThreads();
    
    typedef vigra::TinyVector<float, 2> FCoord;
    auto node = [&](const Coord & coord){
        return size_t(coord[0] + coord[1]*shape[0]);
    };

    std::mutex graphMutex;
    std::mutex mapMutex;

    std::set<size_t > processed;

    struct ToAdd{
        size_t u,v;
        float w;
    };

    size_t bufferSize = 1000000;
    std::vector<std::vector<ToAdd> > buffers(nThreads);
    for(auto & vec : buffers){
        vec.reserve(bufferSize+1);
    }



    auto addToBuffer = [&](std::vector<ToAdd> & buffer,const size_t u_, const size_t v_, const float w_){

        ToAdd ta;
        ta.u=u_;
        ta.v=v_;
        ta.w=w_;
        buffer.push_back(ta);
        //std::cout<<"buffer size"<<buffer.size()<<"\n";
        if(buffer.size()>=bufferSize){

            std::unique_lock<std::mutex> lock(graphMutex);
            //std::cout<<"clear buffer\n";
            for(const auto & ta  : buffer){
                const auto fe = liftedGraph.findEdge(ta.u,ta.v);
                // not yet  in lifted graph
                // therefore cannot be processed
                if(!fe.first){
                    //std::cout<<"a\n";
                    const auto e = model.setCost(ta.u,ta.v,ta.w,false);
                    processed.insert(e);
                }
                // edge is in lifted graph
                else{
                    //std::cout<<"b\n";
                    auto fm = processed.find(fe.second);
                    // not yet processed
                    if(fm == processed.end()){
                        std::cout<<"b1\n";
                        const auto e = model.setCost(ta.u,ta.v,ta.w,false);
                        processed.insert(e);
                    }
                    else{
                        std::cout<<"b2\n";
                    }
                }
            }
            buffer.resize(0);
        }
    };


    std::cout<<"lifted graph edge num "<<liftedGraph.numberOfEdges()<<"\n";
    agraph::parallel_foreach(
        nThreads,
        shape[1],
        [&](const int threadId, const int y){
            auto & buffer = buffers[threadId];
            for(int x=0; x<shape[0]; ++x){
                const auto p = Coord(x,y);
                const auto u = node(p);
                GRAPH_CHECK_OP(u,<,originalGraph.numberOfVertices(),"");
                GRAPH_CHECK_OP(u,<,liftedGraph.numberOfVertices(),"");
                const auto start = clipToImg(p-maxRadius);
                const auto end = clipToImg(p+maxRadius+1);
                auto q = Coord();
                for(q[0]=start[0]; q[0]<end[0]; ++q[0])
                for(q[1]=start[1]; q[1]<end[1]; ++q[1]){

                    GRAPH_CHECK_OP(q[0],>=,0,"");
                    GRAPH_CHECK_OP(q[1],>=,0,"");
                    GRAPH_CHECK_OP(q[0],<,shape[0],"");
                    GRAPH_CHECK_OP(q[1],<,shape[1],"");
                    const auto v = node(q);

                    size_t e;
                    if( norm(p-q) < float(minRadius))
                        continue;
                    if(p==q || v>u){
                        continue;
                    }

                    
                    GRAPH_CHECK_OP(v,<,originalGraph.numberOfVertices(),"");
                    GRAPH_CHECK_OP(v,<,liftedGraph.numberOfVertices(),"");
                    const auto qf = FCoord(q);
                    const auto pq = q-p;
                    const auto dist = vigra::norm(pq);
                    const auto step =  pq*(1.0f / (dist * 1.3f + 0.5f));
                    auto pOnLine = FCoord(p);
                    auto noMax = true;
                    auto maxVal = -1.0f*std::numeric_limits<float>::infinity();
                    while(Coord(pOnLine)!=q){
                        //std::cout<<"pol "<<pOnLine<<"\n";
                        auto iCord = Coord(pOnLine);
                        if(iCord != p){
                            noMax = false;
                            maxVal = std::max(edgePmap[iCord], maxVal);
                        }
                        pOnLine += step;
                    }
                    const double p1 = std::max(std::min(maxVal,0.999f),0.001f);
                    const double p0 = 1.0 - p1;
                    auto w = std::log(p0/p1) + beta;

                    addToBuffer(buffer, u,v,w);
                }
            }
        }
    );

    // clear whats left in buffers
    for(const auto & buffer : buffers){
        for(const auto & ta : buffer){
            const auto fe = liftedGraph.findEdge(ta.u,ta.v);
            // not yet  in lifted graph
            // therefore cannot be processed
            if(!fe.first){
                const auto e = model.setCost(ta.u,ta.v,ta.w,false);
                processed.insert(e);
            }
            // edge is in lifted graph
            else{
                auto fm = processed.find(fe.second);
                // not yet processed
                if(fm == processed.end()){
                    const auto e = model.setCost(ta.u,ta.v,ta.w,false);
                    processed.insert(e);
                }
            }
        }
    }
    std::cout<<"lifted graph edge num "<<liftedGraph.numberOfEdges()<<"\n";
}










template<class LiftedMcModel>
void setCosts(
    LiftedMcModel & liftedMcModel,
    vigra::NumpyArray<1, vigra::TinyVector<uint64_t, 2> > uv,
    vigra::NumpyArray<1, float> costs,
    const bool overwrite 
){
    GRAPH_CHECK_OP(uv.size(), == , costs.size(), "shape mismatch: uv and costs have different size");

    const auto u = uv.bindElementChannel(0);
    const auto v = uv.bindElementChannel(1);
    const auto c = vigra::MultiArrayView<1, float>(costs);
    liftedMcModel.setCosts(u.begin(), u.end(), v.begin(), c.begin(), overwrite);
}


template<class LiftedMcModel>
vigra::NumpyAnyArray edgeLabelsToNodeLabels(
    const LiftedMcModel & model,
    vigra::NumpyArray<1, uint8_t> edgeLabels,
    vigra::NumpyArray<1, uint64_t> nodeLabels
){
    vigra::TinyVector<int,1> shape(model.originalGraph().numberOfVertices());
    nodeLabels.reshapeIfEmpty(shape);
    model.getNodeLabels(edgeLabels, nodeLabels);
    return nodeLabels;
}


template<class LiftedMcModel>
vigra::NumpyAnyArray nodeLabelsToEdgeLabels(
    const LiftedMcModel & model,
    vigra::NumpyArray<1, uint64_t> nodeLabels,
    vigra::NumpyArray<1, uint8_t> edgeLabels
){
    vigra::TinyVector<int,1> shape(model.liftedGraph().numberOfEdges());
    edgeLabels.reshapeIfEmpty(shape);
    model.getEdgeLabels(nodeLabels, edgeLabels);
    return edgeLabels;
}


template<class LiftedMcModel>
void fuseGtObjective(
    LiftedMcModel & model,
    vigra::NumpyArray<2, uint64_t> nodeLabels,
    const size_t rr,
    const double beta,
    const bool verbose
){
    if(verbose)
        std::cout<<"nodeLabel.shape "<<nodeLabels.shape()<<"\n";

   
    const auto & originalGraph = model.originalGraph();
    const auto & liftedGraph = model.liftedGraph();
    auto nV = originalGraph.numberOfVertices();
    auto nGt = nodeLabels.shape(1);

    std::vector<std::set<size_t> > extendedNh(nV);
    std::vector<std::set<size_t> > extendedNh2(nV);


    for(auto n=0; n<originalGraph.numberOfVertices(); ++n){
        for(auto iter = originalGraph.verticesFromVertexBegin(n); iter!=originalGraph.verticesFromVertexEnd(n); ++iter){
            extendedNh[n].insert(*iter);
            extendedNh2[n].insert(*iter);
        }
    }

    for(auto r=0;r<2;++r){

        for(auto n=0; n<originalGraph.numberOfVertices(); ++n){
            auto & thisNodeNhSet = extendedNh2[n];
            for(auto iter = originalGraph.verticesFromVertexBegin(n); iter!=originalGraph.verticesFromVertexEnd(n); ++iter){
                auto otherNode = *iter;
                const auto & otherNodeNhSet = extendedNh[otherNode];
                thisNodeNhSet.insert(otherNodeNhSet.begin(), otherNodeNhSet.end());
                thisNodeNhSet.erase(n);
            }
        }
        extendedNh = extendedNh2;
    }

    for(auto u=0; u<originalGraph.numberOfVertices(); ++u){
        const auto & enh = extendedNh[u];
        //std::cout<<"node "<<u<<" |enh| = "<<enh.size()<<"\n";

        const auto uGt = nodeLabels.bindInner(u); 
        for(const auto v : enh){
            const auto vGt = nodeLabels.bindInner(v); 
            auto p1 = 0.0;
            for(size_t gtc=0; gtc<nGt; ++gtc){
                p1 += uGt[gtc] != vGt[gtc] ?  1 : 0;
            }
            p1 /= nGt;
            //if (p1>0.1)
            //    std::cout<<"   p1 "<<p1<<"\n";

            p1 = std::max(std::min(0.999, p1),0.001);
            auto p0 = 1.0 - p1;
            auto w = std::log(p0/p1) + std::log((1.0-beta)/beta);
            model.setCost(u,v,w);
        }
    }
}



template<class LiftedMcModel>
void fuseGtObjectiveGrid(
    LiftedMcModel & model,
    vigra::NumpyArray<3, uint64_t> nodeLabels,
    vigra::NumpyArray<1, double>    pExpert,
    vigra::NumpyArray<2, float>     cutRegularizer,
    const size_t rr,
    const double beta,
    //const double cLocal,
    const bool verbose
){
    if(verbose)
        std::cout<<"nodeLabel.shape "<<nodeLabels.shape()<<"\n";

   
    const auto & originalGraph = model.originalGraph();
    const auto & liftedGraph = model.liftedGraph();
    auto nV = originalGraph.numberOfVertices();
    typedef vigra::TinyVector<int,2> Coord;
    Coord shape(originalGraph.shape(0), originalGraph.shape(1));

    auto nGt = nodeLabels.shape(2);

    std::vector<std::set<size_t> > extendedNh(nV);
    std::vector<std::set<size_t> > extendedNh2(nV);

    auto clipToImg = [&](const Coord & coord){
        auto c = coord;
        for(auto i=0; i<2; ++i){
            c[i] = std::max(0, c[i]);
            c[i] = std::min(int(shape[i]),c[i]);
        }
        return c;
    };
    auto node = [&](const Coord & coord){
        return size_t(coord[0] + coord[1]*shape[0]);
    };

    Coord p;

    for(p[1]=0; p[1]<shape[1]; ++p[1])
    for(p[0]=0; p[0]<shape[0]; ++p[0]){

        const auto u = node(p);
        const auto start = clipToImg(p-int(rr));
        const auto end = clipToImg(p+int(rr)+1);

        auto lp = nodeLabels.bindInner(p[0]).bindInner(p[1]);
        auto q = Coord();

        for(q[1]=start[1]; q[1]<end[1]; ++q[1])
        for(q[0]=start[0]; q[0]<end[0]; ++q[0]){
            const auto v = node(q);
            if(p!=q && u < v){
                auto d = vigra::norm(p-q);
                if(d<=float(rr)){
                    if(d<=1.5)
                        d*=0.5;
                    auto lq = nodeLabels.bindInner(q[0]).bindInner(q[1]);
                    auto p0 = 0.0;
                    auto p1 = 0.0;
                    for(size_t gtc=0; gtc<nGt; ++gtc){
                        p0 += (lp[gtc] == lq[gtc]) ?  pExpert(gtc) : 0.0;
                        p1 += (lp[gtc] != lq[gtc]) ?  pExpert(gtc) : 0.0;
                    }
                    auto Z = p0+p1;
                    p0/=Z;
                    p1/=Z;
                    //if (p1>0.1)
                    //    std::cout<<"   p1 "<<p1<<"\n";

                    p1 = std::max(std::min(0.99999, p1),0.00001);
                    p0 = 1.0 - p1;
                    auto w = std::log(p0/p1) + std::log((1.0-beta)/beta);
                    if(d<=1.01){
                        auto c = (cutRegularizer[p] + cutRegularizer[q])/2.0;
                        w += c;
                    }
                    model.setCost(u,v,w*d);
                }
            }
        }
    }
}




template<class OG, class F>
void exportLiftedMcModelT(const std::string & clsName, F && f){
    typedef OG originalGraph;
    typedef agraph::multicut_lifted::LiftedMcModel<originalGraph, float> LiftedMcModel;


    auto cls = bp::class_<LiftedMcModel>
    (
        clsName.c_str(), 
        bp::init<
            const originalGraph&
        >(
            bp::arg("originalGraph")
        )
    )
    .def("_setCosts",vigra::registerConverters(&setCosts<LiftedMcModel>))
    .def("edgeLabelsToNodeLabels",vigra::registerConverters(&edgeLabelsToNodeLabels<LiftedMcModel>),
        (
            bp::arg("edgeLabels"),
            bp::arg("out") = bp::object()
        )
    )
    .def("nodeLabelsToEdgeLabels",vigra::registerConverters(&nodeLabelsToEdgeLabels<LiftedMcModel>),
        (
            bp::arg("nodeLabels"),
            bp::arg("out") = bp::object()
        )
    );
    f(cls);


    bp::def("fuseGtObjective",vigra::registerConverters(&fuseGtObjective<LiftedMcModel>),
        (
            bp::arg("model"),
            bp::arg("nodeLabels"),
            bp::arg("rr") = 2,
            bp::arg("beta") = 0.5,
            bp::arg("verbose") = true
        )
    );


}


template<class MODEL>
vigra::NumpyAnyArray flattenLabels(
    const MODEL & model,
    vigra::NumpyArray<2, uint64_t> labels2d,
    vigra::NumpyArray<1, uint64_t> out
){
    vigra::TinyVector<int, 1> shape(labels2d.size());
    out.reshapeIfEmpty(shape);
    auto c=0;
    for(auto y=0; y<labels2d.shape(1); ++y)
    for(auto x=0; x<labels2d.shape(0); ++x){
        out[c] = labels2d(x,y);
        ++c;
    }
    return out;
}


void exportLiftedMcModel(){
    typedef agraph::GridGraph<2> GridGraph2D;
    typedef agraph::GridGraph<3> GridGraph3D;
    typedef agraph::Graph<> Graph;
    typedef agraph::multicut_lifted::LiftedMcModel<GridGraph2D, float> LiftedMcModelGridGraph2D;

    {
        typedef agraph::multicut_lifted::LiftedMcModel<GridGraph2D, float> LiftedMcModelGridGraph2D; 
        exportLiftedMcModelT<GridGraph2D>("LiftedMcModelGridGraph2D",
            [&](
                bp::class_< LiftedMcModelGridGraph2D > & cls
            ){
                cls
                    .def("flattenLabels",
                        vigra::registerConverters(&flattenLabels< LiftedMcModelGridGraph2D >),
                        (
                            bp::arg("labels"),
                            bp::arg("out") =  bp::object()
                        )
                    )
                ;
            }
        );
    }
    exportLiftedMcModelT<GridGraph3D>("LiftedMcModelGridGraph3D",
        [&](
            bp::class_<agraph::multicut_lifted::LiftedMcModel<GridGraph3D, float> > & cls
        ){
            
        }
    );
    exportLiftedMcModelT<Graph>("LiftedMcModelGraph",
        [&](
            bp::class_<agraph::multicut_lifted::LiftedMcModel<Graph, float> > & cls
        ){
            
        }
    );


    bp::def("addLongRangeEdges",vigra::registerConverters(&addLongRangeEdges<LiftedMcModelGridGraph2D>))
    ;


    bp::def("fuseGtObjectiveGrid",vigra::registerConverters(&fuseGtObjectiveGrid<LiftedMcModelGridGraph2D>),
        (
            bp::arg("model"),
            bp::arg("nodeLabels"),
            bp::arg("pExpert"),
            bp::arg("cutRegularizer"),
            bp::arg("rr") = 2,
            bp::arg("beta") = 0.5,
            //bp::arg("cLocal") = 100,
            bp::arg("verbose") = true
        )
    );
}
