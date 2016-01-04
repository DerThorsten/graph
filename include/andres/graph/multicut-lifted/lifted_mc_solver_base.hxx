#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_LIFTED_SOLVER_BASE_HXX
#define ANDRES_GRAPH_MULTICUT_LIFTED_SOLVER_BASE_HXX

#include <iomanip>
#include <stdexcept>
#include <unordered_set>
#include <vector>
#include <stack>
#include <cmath>
#include <random>
#include <mutex>


namespace visitors{

    template<class MODEL>
    class VisitorBase{
    public:

    private:

    };

    template<class MODEL>
    class VerboseVisitor{
    public:

    private:

    };

}


template<class MODEL>
class LiftedMulticutSolverBase{
public:
    typedef MODEL Model;
    typedef LiftedMulticutSolverBase<MODEL> VisitorBase;
    typedef std::vector<uinEdgeLabels;


    // to implement
    virtual void run(VisitorBase * visitor) = 0;
    void argmin(uint8_t * result) = 0;

    // with default impl
    virtual void run(){
        return this->run(nullptr)
    }


    virtual void setStartingPoint(uint8_t * result){
        
    }


private:

};


#endif // ANDRES_GRAPH_MULTICUT_LIFTED_SOLVER_BASE_HXX
