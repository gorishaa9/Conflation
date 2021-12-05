#ifndef COUNTHANDLER_H
#define COUNTHANDLER_H

#include <osmium/handler.hpp>
#include <boost/shared_ptr.hpp>
#include <osmium/osm/object.hpp>

#include <map>
#include <iostream>

#include "terminal_node.hpp"

using namespace std;

using key_type = osmium::object_id_type;
using value_type = int;

extern vector<std::map<key_type, value_type>> nodeset;
extern vector<std::map<key_type, boost::shared_ptr<Terminal_node>>> terminal_nodes;

nodeset.push_back(ini1);
nodeset.push_back(ini2);

extern int dataset;

struct CountHandler : public osmium::handler::Handler{
public:
  void way(const osmium::Way& ways){
    int count = 1;                        //count of the number of nodes in ways
    int size = ways.nodes().size();       //total number of nodes in the way
    for(const auto& n : ways.nodes()){
        //std::cout << n.lon() << " ";
        auto it = nodeset[dataset].find(n.ref());
        if(it != nodeset[dataset].end()){
            int val = it->second;
            (count==1 || count==size) ? it->second = val+1 : it->second = val+2;
        }
        else{
            (count==1 || count==size) ? nodeset[dataset].insert({n.ref(), 1}): nodeset[dataset].insert({n.ref(), 2});
        }
        ++count;
    }
  }
};

#endif
