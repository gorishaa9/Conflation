#ifndef RELATION_H
#define RELATION_H

#include <vector>
#include <boost/shared_ptr.hpp>

#include "Tuple.hpp"

using namespace std;

class Relation{

private:
    long long int _id;
    vector<boost::shared_ptr<Tuple>> _tuples;

public:
    Relation(long long int id) :
            _id(id) {
            }

    long long int get_id(){
        return _id;
    }

    void insert(boost::shared_ptr<Tuple> tuple){
        _tuples.push_back(tuple);
    }

    vector<boost::shared_ptr<Tuple>>& get_entries(){
        return _tuples;
    }

    size_t size(){
        return _tuples.size();
    }
};

#endif
