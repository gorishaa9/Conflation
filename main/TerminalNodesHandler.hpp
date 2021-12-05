#include <osmium/osm/object.hpp>
#include <boost/shared_ptr.hpp>

#include <fstream>

#include "terminal_node.hpp"
#include "CountHandler.hpp"

ofstream file("terminal_nodes.txt");

struct TerminalNodesHandler : public osmium::handler::Handler {
public:
    void way(const osmium::Way& ways){
      int size = ways.nodes().size();
      int count=1;
      for(const auto& n : ways.nodes()){
        int deg = nodeset[dataset].find(n.ref())->second;

        if(count==1 || count==size || deg==1 || deg>=3){
          boost::shared_ptr<Terminal_node> ptr(new Terminal_node(n.ref(), Coordinates{n.location()}, deg));
          terminal_nodes[dataset].insert({n.ref(), ptr});
          //file << n.ref() << "\t" << deg << "\t" << ptr->get_location().x << "\t" << ptr->get_location().y << "\n";
        }
        ++count;
      }
    }
};
