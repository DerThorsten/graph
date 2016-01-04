#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_PARALLEL_LIFTED_MC_HXX
#define ANDRES_GRAPH_MULTICUT_PARALLEL_LIFTED_MC_HXX

#include <iomanip>
#include <stdexcept>
#include <unordered_set>
#include <vector>
#include <stack>
#include <cmath>
#include <random>
#include <mutex>

#include <vigra/resizeimage.hxx>
#include <vigra/convolution.hxx>

#include "andres/graph/grid-graph.hxx"
#include "andres/graph/threadpool.hxx"
#include "andres/graph/multicut-lifted/fusion_mover.hxx"
#include "andres/graph/multicut-lifted/greedy-additive.hxx"
#include "andres/graph/multicut-lifted/kernighan-lin.hxx"

namespace andres {
namespace graph {
namespace multicut_lifted {


/*
    decompose grid into overlapping blocks
 
    solve each block on its on
 
    past beach block solution in a global solution
 
    fuse them via fusion moves
    
*/


class SolverSettings{
    enum Settings{
        KL = 0,
        GA = 1,
        FMR = 2
    };
};

template<class LIFTED_MC_MODEL>
class RandomizedProposalGenerator;

template<class LIFTED_MC_MODEL>
class RescaledProposalGenerator;

template<class LIFTED_MC_MODEL>
class SubgraphProposalGenerator;

template<class LIFTED_MC_MODEL, class PROPOSAL_GEN = SubgraphProposalGenerator<LIFTED_MC_MODEL> >
class ParallelSolver;

template<class LIFTED_MC_MODEL>
class SubgraphProposalGenerator{
public:
    typedef GridGraph<2> GridGraph2D;
    typedef typename GridGraph2D::VertexCoordinate Coord2d;
    typedef std::vector<uint8_t>  LiftedEdgeLabels;
    struct  Settings{
        size_t subgraphRadius = {40};
        int seed = {-1};
    };

    SubgraphProposalGenerator(const LIFTED_MC_MODEL & model, const Settings & settings)
    :   model_(model),
        settings_(settings),
        rd_(),
        gen_(settings.seed == -1 ? std::mt19937(rd_()): std::mt19937(size_t(settings.seed)) ),
        //gen_(),
        nodeDist_(0,model.originalGraph().numberOfVertices()-1)
    {

    }


    void generate(const LiftedEdgeLabels & currentBest, LiftedEdgeLabels & proposal){

        this->getSubgraph(model_.originalGraph(), currentBest, proposal );
    }

private:

    template<class SUB_MODEL>
    void optimizeSubmodel(const SUB_MODEL & subModel,
                          std::vector<uint8_t> & subgraphRes){

        // hard work is here
        std::vector<uint8_t> subgraphInput(subModel.liftedGraph().numberOfEdges(),1);
        subgraphRes.resize(subModel.liftedGraph().numberOfEdges());
        KernighanLinSettings s;
        s.verbose = false;
        kernighanLin(subModel.originalGraph(),subModel.liftedGraph(),
                     subModel.edgeCosts(),subgraphInput, subgraphRes,s);

        // transfer back to node solution
        

    }

    template<class G>
    void getSubgraph(const G & g, const LiftedEdgeLabels & currentBest, LiftedEdgeLabels & proposal){
    }

    void getSubgraph(const GridGraph<2> & g, const LiftedEdgeLabels & currentBest, LiftedEdgeLabels & proposal){
        Coord2d ggStart,subShape;   
        std::tie(ggStart, subShape) = this->getRandSubGridGraph(g);
        GridGraph2D gg(subShape);

        LiftedMcModel<GridGraph2D, float> subModel(gg);

        const auto & liftedGraph = model_.liftedGraph();
        const auto & edgeCosts  =model_.edgeCosts();

        auto uLocal = 0;
        for(auto ly=0; ly<subShape[1]; ++ly)
        for(auto lx=0; lx<subShape[0]; ++lx){
            auto ugx = lx + ggStart[0];
            auto ugy = ly + ggStart[1];
            auto uGlobal = ugx + ugy*g.shape(0);

            // all edges of global node where
            // both, u and v are in subgraph
            auto aIter  = liftedGraph.adjacenciesFromVertexBegin(uGlobal);
            auto aEnd  = liftedGraph.adjacenciesFromVertexEnd(uGlobal);
            while(aIter != aEnd){
                auto a = *aIter;
                auto e = a.edge();
                auto v = a.vertex();
                auto vgy = v/ g.shape(0);
                auto vgx = v - vgy*g.shape(1);
                if(vgx>=ggStart[0] && vgx <ggStart[0] + subShape[0] &&
                   vgy>=ggStart[1] && vgy <ggStart[1] + subShape[1]){
                    auto vlx = vgx - ggStart[0];
                    auto vly = vgy - ggStart[1];
                    auto vLocal  = vlx + vly*subShape[0];
                    subModel.setCost(uLocal,vLocal, edgeCosts[e]);
                }
                ++aIter;
            }
            ++uLocal;
        }
        LiftedEdgeLabels subLabls;
        this->optimizeSubmodel(subModel, subLabls);

        uLocal = 0;
        for(auto ly=0; ly<subShape[1]; ++ly)
        for(auto lx=0; lx<subShape[0]; ++lx){
            auto ugx = lx + ggStart[0];
            auto ugy = ly + ggStart[1];
            auto uGlobal = ugx + ugy*g.shape(0);

            // all edges of global node where
            // both, u and v are in subgraph
            auto aIter  = liftedGraph.adjacenciesFromVertexBegin(uGlobal);
            auto aEnd  = liftedGraph.adjacenciesFromVertexEnd(uGlobal);
            while(aIter != aEnd){
                auto a = *aIter;
                auto e = a.edge();
                auto v = a.vertex();
                auto vgy = v/ g.shape(0);
                auto vgx = v - vgy*g.shape(1);
                if(vgx>=ggStart[0] && vgx <ggStart[0] + subShape[0] &&
                   vgy>=ggStart[1] && vgy <ggStart[1] + subShape[1]){
                    auto vlx = vgx - ggStart[0];
                    auto vly = vgy - ggStart[1];
                    auto vLocal  = vlx + vly*subShape[0];
                    auto eLocal = subModel.liftedGraph().findEdge(vLocal,  uLocal).second;
                    proposal[e] = subLabls[eLocal];
                }
                else{
                    proposal[e] = 1;
                }
                ++aIter;
            }
            ++uLocal;
        }

    }


    std::pair<Coord2d, Coord2d>
    getRandSubGridGraph(const GridGraph<2> & g){

        const size_t randVar = nodeDist_(gen_);
        const int y = randVar / g.shape(0);
        const int x = randVar - y*g.shape(0);
        const int r = settings_.subgraphRadius;
        const int startX = std::max(0, x - r);
        const int startY = std::max(0, y - r);
        const int stopX = std::min(int(g.shape(0)), x+ r + 1);
        const int stopY = std::min(int(g.shape(1)), y+ r + 1);
        Coord2d start = {size_t(startX),size_t(startY)};
        Coord2d shape = {size_t(stopX-startX),size_t(stopY-startY)};
        return std::pair<Coord2d, Coord2d>(start,shape);
    }


    const LIFTED_MC_MODEL & model_;
    Settings settings_;

    // rand gen
    std::random_device rd_;
    std::mt19937 gen_;
    std::uniform_int_distribution<> nodeDist_;
};

template<class LIFTED_MC_MODEL>
class RescaledProposalGenerator{
public:
    typedef GridGraph<2> GridGraph2D;
    typedef typename GridGraph2D::VertexCoordinate Coord2d;
    typedef std::vector<uint8_t>  LiftedEdgeLabels;
    struct  Settings{
        float reducingFactorMean = {3.0};
        float reducingFactorSigma = {2.0};
        int seed  = {-1};
    };

    RescaledProposalGenerator(const LIFTED_MC_MODEL & model, const Settings & settings)
    :   model_(model),
        settings_(settings),
        rd_(),
        gen_(settings.seed == -1 ? std::mt19937(rd_()): std::mt19937(size_t(settings.seed)) ),
        facDist_(settings.reducingFactorMean,settings.reducingFactorSigma)
    {

    }


    void generate(const LiftedEdgeLabels & currentBest, LiftedEdgeLabels & proposal){
        this->rescale(model_.originalGraph(), proposal);
    }

    void rescale(const GridGraph2D & graph,LiftedEdgeLabels & proposal){

        

        const auto & liftedGraph = model_.liftedGraph();
        const auto & originalGraph = model_.originalGraph();
        const auto & edgeCosts = model_.edgeCosts();
        auto shapeX = originalGraph.shape(0);
        auto shapeY = originalGraph.shape(1);
        auto fShapeX = float(shapeX);
        auto fShapeY = float(shapeY);
        auto f = facDist_(gen_);
        f = std::max(1.5,f);
        auto iShapeX = int(fShapeX/f + 0.5);
        auto iShapeY = int(fShapeY/f + 0.5);
        
        vigra::TinyVector<int, 2> shape(shapeX,shapeY);
        vigra::TinyVector<int, 2> ishape(iShapeX,iShapeY);
        vigra::MultiArray<2, float> hval(shape,0.0);
        vigra::MultiArray<2, float> vval(shape,0.0);
        vigra::MultiArray<2, float> ihval(ishape,0.0);
        vigra::MultiArray<2, float> ivval(ishape,0.0);
        auto node = [&](const int x,const int y){
            return x+y*shapeX;
        };
        for(auto y=0; y<shapeY; ++y)
        for(auto x=0; x<shapeX; ++x){
            if(x+1<shapeX){
                const auto e = liftedGraph.findEdge(node(x,y),node(x+1,y)).second;
                hval(x,y) += 0.5f *model_.edgeCosts()[e];
                hval(x+1,y) += 0.5f *model_.edgeCosts()[e];
            }
            if(y+1<shapeY){
                const auto e = liftedGraph.findEdge(node(x,y),node(x,y+1)).second;
                vval(x,y) += 0.5f * model_.edgeCosts()[e];
                vval(x,y+1) += 0.5f * model_.edgeCosts()[e];
            }
        }

        vigra::MultiArray<2, float> hvals(hval.shape());
        vigra::MultiArray<2, float> vvals(vval.shape());

        vigra::gaussianSmoothing(hval, hvals, 1.0);
        vigra::gaussianSmoothing(vval, vvals, 1.0);

        vigra::resizeImageSplineInterpolation(hvals, ihval, vigra::BSpline<2, float>());
        vigra::resizeImageSplineInterpolation(vvals, ivval, vigra::BSpline<2, float>());


        GridGraph2D iGridGraph({std::size_t(iShapeX), std::size_t(iShapeY)});
        LiftedMcModel<GridGraph2D,float> iModel(iGridGraph);


        // add long range costs
        auto inode = [&](const vigra::TinyVector<int,2> & coord_){
            return coord_[0]+coord_[1]*iShapeX;
        };

        auto fac = shape/ishape;

        for(auto iy=0; iy<iShapeY; ++iy)
        for(auto ix=0; ix<iShapeX; ++ix){
            if(ix+1<iShapeX){
                const auto val = 0.5f*(ihval(ix,iy) + ihval(ix+1,iy));
                iModel.setCost(
                    inode(vigra::TinyVector<int,2>(ix,iy) ),
                    inode(vigra::TinyVector<int,2>(ix+1,iy) ),
                    val
                );
            }
            if(iy+1<iShapeY){
                const auto val = 0.5f*(ivval(ix,iy) + ivval(ix,iy+1));
                iModel.setCost(
                    inode(vigra::TinyVector<int,2>(ix,iy) ),
                    inode(vigra::TinyVector<int,2>(ix,iy+1) ),
                    val
                );
            }
        }
        for(auto edge =0; edge < liftedGraph.numberOfEdges(); ++edge){
            auto v0 = liftedGraph.vertexOfEdge(edge, 0);
            auto v1 = liftedGraph.vertexOfEdge(edge, 1);
            if(!originalGraph.findEdge(v0,v1).first){

                auto v0y = v0/shapeX;
                auto v0x = v0 - v0y*shapeX;
                auto v1y = v1/shapeX;
                auto v1x = v1 - v1y*shapeX;

                auto fic0 = vigra::TinyVector<float,2>(v0x, v0y)/fac + 0.5f;
                auto fic1 = vigra::TinyVector<float,2>(v1x, v1y)/fac + 0.5f;
                auto ic0 = vigra::TinyVector<int,2>(fic0);
                auto ic1 = vigra::TinyVector<int,2>(fic1);
                if(ic0 != ic1 && ic0[0]>=0 && ic0[0]<iShapeX && ic0[1]>=0 && ic0[1]<iShapeY && 
                                 ic1[0]>=0 && ic1[0]<iShapeX && ic1[1]>=0 && ic1[1]<iShapeY){
                    //std::cout<<"ic0 "<<ic0<<" ic1 "<<ic1<<"\n";
                    auto iu = inode(ic0);
                    auto iv = inode(ic1);
                    iModel.setCost(iu,iv,edgeCosts[edge]);
                }
            }
        }

        // solve model on this scale
        std::vector<uint8_t> subgraphInput(iModel.liftedGraph().numberOfEdges(),1);
        std::vector<uint8_t> subgraphRes(iModel.liftedGraph().numberOfEdges(),1);
        subgraphRes.resize(iModel.liftedGraph().numberOfEdges());
        KernighanLinSettings s;
        s.verbose = false;
        kernighanLin(iModel.originalGraph(),iModel.liftedGraph(),
                     iModel.edgeCosts(),subgraphInput, subgraphRes,s);



        struct SubgraphWithCut { // a subgraph with cut mask
            SubgraphWithCut(const std::vector<uint8_t> & labels, 
                            std::vector<std::size_t> const& edge_in_lifted_graph)
                : labels_(labels), edge_in_lifted_graph_(edge_in_lifted_graph)
                {}
            bool vertex(const std::size_t v) const
                { return true; }
            bool edge(const std::size_t e) const
                { return labels_[edge_in_lifted_graph_[e]] == 0; }

            std::vector<std::size_t> const& edge_in_lifted_graph_;
            const std::vector<uint8_t> & labels_;
        };

        std::vector<size_t> edgeInLiftedGraph(iModel.originalGraph().numberOfEdges());
        for (std::size_t i = 0; i < iModel.originalGraph().numberOfEdges(); ++i){
            auto v0 = iModel.originalGraph().vertexOfEdge(i, 0);
            auto v1 = iModel.originalGraph().vertexOfEdge(i, 1);
            edgeInLiftedGraph[i] = iModel.liftedGraph().findEdge(v0, v1).second;
        } 

        vigra::MultiArray<2, uint64_t> inodeLabels(ishape);
        vigra::MultiArray<2, uint64_t> nodeLabels(shape);

        ComponentsBySearch<GridGraph2D > components;
        components.build(iModel.originalGraph(), SubgraphWithCut(subgraphRes, edgeInLiftedGraph));
        for(std::size_t n=0; n<iModel.originalGraph().numberOfVertices(); ++n){
            auto iy = n/iShapeX;
            auto ix = n - iy*iShapeX;
            auto l = components.labels_[n];
            inodeLabels(ix,iy) = l;
        }

        vigra::resizeImageNoInterpolation(inodeLabels, nodeLabels);

        for(auto edge =0; edge < liftedGraph.numberOfEdges(); ++edge){
            auto v0 = liftedGraph.vertexOfEdge(edge, 0);
            auto v1 = liftedGraph.vertexOfEdge(edge, 1);
            auto v0y = v0/shapeX;
            auto v0x = v0 - v0y*shapeX;
            auto v1y = v1/shapeX;
            auto v1x = v1 - v1y*shapeX;
            proposal[edge] = nodeLabels(v0x,v0y) != nodeLabels(v1x,v1y)  ? 1 : 0 ;
        }

        
    }
    template<class GRAPH>
    void rescale(const GRAPH & graph,LiftedEdgeLabels & proposal){

    }

private:

   

    const LIFTED_MC_MODEL & model_;
    Settings settings_;

    // rand gen
    std::random_device rd_;
    std::mt19937 gen_;
    std::uniform_real_distribution<> facDist_;
};


template<class LIFTED_MC_MODEL>
class RandomizedProposalGenerator{
public:
    typedef GridGraph<2> GridGraph2D;
    typedef typename GridGraph2D::VertexCoordinate Coord2d;
    typedef std::vector<uint8_t>  LiftedEdgeLabels;
    struct  Settings{
        double sigma = {15.0};
        double nodeLimit = {0.05};
        int seed  = {-1};
        bool useGA = {true};
    };

    RandomizedProposalGenerator(const LIFTED_MC_MODEL & model, const Settings & settings)
    :   model_(model),
        settings_(settings),
        nEdgeCosts_(model.liftedGraph().numberOfEdges()),
        rd_(),
        gen_(settings.seed == -1 ? std::mt19937(rd_()): std::mt19937(size_t(settings.seed)) ),
        nDist_(0,settings.sigma)
    {

    }


    void generate(const LiftedEdgeLabels & currentBest, LiftedEdgeLabels & proposal){

        const auto & originalGraph = model_.originalGraph();
        const auto & liftedGraph = model_.liftedGraph();
        const auto & edgeCosts = model_.edgeCosts();

        for(auto e=0; e<model_.liftedGraph().numberOfEdges(); ++e){
            nEdgeCosts_[e] =edgeCosts[e]+nDist_(gen_);
        }

        if(settings_.useGA){
            int stopCond = -1;
            if(settings_.nodeLimit<1.0)
                stopCond = double(originalGraph.numberOfVertices())*settings_.nodeLimit + 0.5;
            else
                stopCond = settings_.nodeLimit + 0.5;
            greedyAdditiveEdgeContraction(originalGraph,liftedGraph,nEdgeCosts_, proposal, stopCond);
        }
        else{
            KernighanLinSettings s;
            s.verbose = false;
            std::vector<uint8_t> ones(currentBest.size(), 1);
            kernighanLin(originalGraph,liftedGraph,nEdgeCosts_,ones, proposal,s);

        }
    }



private:

   

    const LIFTED_MC_MODEL & model_;
    Settings settings_;

    // rand gen
    std::random_device rd_;
    std::mt19937 gen_;
    std::normal_distribution<> nDist_;

    std::vector<float> nEdgeCosts_;
};



template<class LIFTED_MC_MODEL, class PROPOSAL_GEN>
class ParallelSolver{
public:
    typedef LIFTED_MC_MODEL LiftedMcModel;
    typedef PROPOSAL_GEN ProposalGen;
    typedef FusionMover<
                typename LiftedMcModel::OriginalGraph,
                typename LiftedMcModel::LiftedGraph,
                typename LiftedMcModel::EdgeCosts
    > Fm;

    typedef typename ProposalGen::Settings ProposalGenSettings;

    struct Settings{
        std::size_t maxNumberOfIterations{4};
        std::size_t nParallelProposals { 100 };
        std::size_t reduceIterations {1};
        int seed{-1};
        std::size_t verbose{1};
        ProposalGenSettings proposalsGenSettings {ProposalGenSettings()};

        std::vector< std::vector<uint8_t> > externalProposals_;
    };

    ParallelSolver(const LiftedMcModel & model, const Settings & settings = Settings())
    :   model_(model),
        settings_(settings),
        bestEdgeLabels_(model.liftedGraph().numberOfEdges(),0){
    }
    

    template<class LABELS_IN, class LABELS_OUT>
    void run(
        const LABELS_IN & inputEdgeLabels,
        LABELS_OUT & outputLabels
    ){
        // shortcuts
        const auto & originalGraph = model_.originalGraph();
        const auto & liftedGraph = model_.liftedGraph();
        const auto & edgeCosts = model_.edgeCosts();

        // lambda to compute energy 
        auto getEnergy = [&] (const std::vector<uint8_t> & edgeLabels_) {
            auto totalCost = 0.0;
            for(std::size_t edge = 0; edge < liftedGraph.numberOfEdges(); ++edge){
                if(edgeLabels_[edge]){
                    totalCost += edgeCosts[edge];
                }
            }
            return totalCost;
        };


        // setup threadpool
        auto pOpt = ParallelOptions();
        const auto nThreads = pOpt.numThreads(-1).getActualNumThreads();
        ThreadPool threadpool(pOpt);

        // store all proposals of one round
        auto nExternalProposal = settings_.externalProposals_.size();
        std::vector< std::vector<uint8_t>  > proposals(settings_.nParallelProposals + nExternalProposal,
                std::vector<uint8_t>(liftedGraph.numberOfEdges(),0));

        for(auto i=0; i<nExternalProposal; ++i){
            proposals[i+settings_.nParallelProposals] = settings_.externalProposals_[i];
        }
        std::cout<<"nThreads "<<nThreads<<"\n";
        // get as many proposal generators as there are threads
        std::vector<ProposalGen*> proposalGens(nThreads, nullptr);
        auto pgs = 0;
        for(auto & pg : proposalGens){
            auto s = settings_.proposalsGenSettings;
            if(settings_.seed == -1)
                s.seed = -1;
            else
                s.seed  = settings_.seed + pgs;
            pg = new ProposalGen(model_, s);
        }
        




        std::vector<Fm*> fms(nThreads);
        for(auto & fm : fms){
            fm = new Fm(originalGraph, liftedGraph, model_.edgeCosts());
        }
        
        // outer iterations
        for(auto outerIter=0; outerIter<settings_.maxNumberOfIterations; ++outerIter){

            ////////////////////////////////////////////////
            // Generate proposals in parallel
            ////////////////////////////////////////////////
            {
                auto nParallelProposals = settings_.nParallelProposals;
                std::function<void(int)> cb;
                if(settings_.verbose >=1)
                    std::cout<<"generate proposals..\n";
                ProgressCallback progressCallback(nParallelProposals, cb);
                progressCallback.setVerboseCallback(settings_.verbose);
                parallel_foreach(threadpool, nParallelProposals,[&](const int threadId, const int p){

                    // generate the proposal
                    auto & proposal = proposals[p];
                    std::fill(proposal.begin(), proposal.end(), 0);
                    proposalGens[threadId]->generate(bestEdgeLabels_, proposal);
                    
                    // report progress
                    progressCallback.increment(1);

                });
                if(settings_.verbose >=1);
                    std::cout<<"...done\n";
            }


            std::vector< std::vector<uint8_t> > toFuse;
            std::vector< std::vector<uint8_t> > toFuse2;
            std::mutex toFuseMutex;

            /////////////////////////////////////////////////////
            /// reduce
            /// ///////////////////////////////////////////////////
            std::cout<<"reduce\n";
            toFuse = proposals;
            auto level = 0;
            while(level<settings_.reduceIterations && true){
                //std::cout<<toFuse.size()<<"\n";
                if(toFuse.size() == 1){
                    break;
                }
                if(toFuse.size()%2 == 1){
                    toFuse2.push_back(toFuse.back());
                    toFuse.pop_back();
                }
                auto nJobs = toFuse.size()/2;
                parallel_foreach(threadpool, nJobs  ,[&](const int threadId, const int i){
                    auto fm = fms[threadId];
                    auto & pa = toFuse[i*2];
                    auto & pb = toFuse[i*2+1];
                    for(size_t jj=0; jj<pa.size(); ++jj){
                        pa[jj] = std::max(pa[jj],pb[jj]);
                    }
                    {
                        std::unique_lock<std::mutex> lock(toFuseMutex);
                        toFuse2.push_back(pa);
                    }
                });

                toFuse = toFuse2;
                toFuse2.resize(0);
                ++level;
            }



            // fuse all proposals with best in parallel
            // (but only if we are not in the first iteration)
            if(outerIter>=1)
            {
                if(settings_.verbose >=1)
                    std::cout<<"fuse proposals..\n";
                ProgressCallback progressCallback(toFuse.size());
                progressCallback.setVerboseCallback(settings_.verbose);
                parallel_foreach(threadpool, toFuse.size(),[&](const int threadId, const int p){
                    auto fm = fms[threadId];
                    auto & proposal = toFuse[p];
                    std::vector<uint8_t> best(bestEdgeLabels_.size());
                    std::vector<uint8_t> fuseResults(liftedGraph.numberOfEdges());
                    bool changes = fm->fuse(bestEdgeLabels_, proposal, fuseResults);
                    //std::cout<<"ePost "<<ePost<<"\n";
                    if(std::any_of(fuseResults.cbegin(), fuseResults.cend(), [](uint8_t val_){ return val_==1; }))
                    {
                        {
                            std::unique_lock<std::mutex> lock(toFuseMutex);
                            toFuse2.push_back(fuseResults);
                            const auto ePost = getEnergy(fuseResults);
                            //std::cout<<"ePost "<<ePost<<"\n";
                        }
                    }
                    // report progress
                    progressCallback.increment(1);
                });
                if(settings_.verbose >=1);
                    std::cout<<"...done\n";
                toFuse = toFuse2;
                toFuse2.resize(0);
            }
            else{
                //toFuse = proposals;
            }
            //std::cout<<"hierarchical Fusion "<<toFuse.size()<<"\n";
            level = 0;
            while(true){
                //std::cout<<toFuse.size()<<"\n";
                if(toFuse.size() == 1){
                    break;
                }
                if(toFuse.size()%2 == 1){
                    toFuse2.push_back(toFuse.back());
                    toFuse.pop_back();
                }
                auto nJobs = toFuse.size()/2;
                if(settings_.verbose >=1)
                    std::cout<<"hierarchical..\n";
                ProgressCallback progressCallback(nJobs);
                progressCallback.setVerboseCallback(settings_.verbose);
                parallel_foreach(threadpool, nJobs  ,[&](const int threadId, const int i){
                    auto fm = fms[threadId];
                    auto & pa = toFuse[i*2];
                    auto & pb = toFuse[i*2+1];
                    std::vector<uint8_t> fuseResults(liftedGraph.numberOfEdges());
                    bool changes = fm->fuse(pa, pb, fuseResults);
                    {
                        std::unique_lock<std::mutex> lock(toFuseMutex);
                        toFuse2.push_back(fuseResults);
                    }
                    progressCallback.increment(1);
                });

                toFuse = toFuse2;
                toFuse2.resize(0);
            }
            bestEdgeLabels_ = toFuse[0];
            const auto ePost = getEnergy(bestEdgeLabels_);
            std::cout<<"ePost "<<ePost<<"\n";
            
        }

        for(auto pg : proposalGens)
            delete pg;
        for(auto fm : fms)
            delete fm;
        std::copy(bestEdgeLabels_.begin(), bestEdgeLabels_.end(), outputLabels.begin());

    }
    const LiftedMcModel & getModel()const{
        return model_;
    }
private:


    LiftedMcModel model_;
    Settings settings_;
    std::vector<uint8_t> bestEdgeLabels_;
};







} // end namespace multicut-lifted
} // end namespace graph
} // end namespace andres

#endif /*ANDRES_GRAPH_MULTICUT_PARALLEL_LIFTED_MC_HXX*/
