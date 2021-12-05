#ifndef INITIALIZE_HPP
#define INITIALIZE_HPP

#include "CountHandler.hpp"
#include "SegmentHandler.hpp"

std::map<key_type, value_type> nodeset_ini1;
std::map<key_type, value_type> nodeset_ini2;
std::map<key_type, boost::shared_ptr<Terminal_node>> tn_ini1;
std::map<key_type, boost::shared_ptr<Terminal_node>> tn_ini2;
std::map<SKey, SVal> seg_ini1;
std::map<SKey, SVal> seg_ini2;

class Initialize{

public:
    void initialize(){
        nodeset.push_back(nodeset_ini1);
        nodeset.push_back(nodeset_ini2);
        terminal_nodes.push_back(tn_ini1);
        terminal_nodes.push_back(tn_ini2);
        segments.push_back(seg_ini1);
        segments.push_back(seg_ini2);
    }
};

#endif
