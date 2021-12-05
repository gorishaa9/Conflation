#ifndef SEGMENTHANDLER_H
#define SEGMENTHANDLER_H

#include "SegmentUtils.hpp"
#include "IntNode.hpp"

#define ce 5
#define PI 3.14159265

using namespace std;
using namespace std::this_thread; // sleep_for, sleep_until
using namespace std::chrono; // nanoseconds, system_clock, seconds

using namespace geos;
using namespace geos::geom;
using namespace osmium::geom;

using SKey = long long int;
using SVal = boost::shared_ptr<Segment_creator>;
using key_type = osmium::object_id_type;
using value_type = int;

ofstream ref_data("ref_data.txt");
ofstream tar_data("tar_data.txt");

int flag, ctr;
double length;

char* data_source;

double minx;
double maxx;
double miny;
double maxy;

CoordinateArraySequence* cl;
vector<IntNode> nodeseq;
SVal previous;

map< SKey, map< SKey, SKey >> small_ways;
vector< double > total_length;

struct RoadLengthHandler : public osmium::handler::Handler {

    // If the way has a "highway" tag, find its length and add it to the
    // overall length.

    void way(const osmium::Way& ways) {
        if( strcmp( data_source, "OSM" ) == 0 && !ways.tags().has_key(key) )
            return;

          double length = osmium::geom::haversine::distance(ways.nodes());
          if( dataset == 0 ){
              // cout << "length ... " << length << "\n";
              if( length <= 1 ){
                  small_ways[ ways.nodes().front().ref() ][ ways.nodes().back().ref() ] = ways.id();
                  small_ways[ ways.nodes().back().ref() ][ ways.nodes().front().ref() ] = ways.id();
              }
          }
          total_length[ dataset ] += length;
    }

};

struct SegmentHandler : public osmium::handler::Handler{

    void way(const osmium::Way& ways){
      if( strcmp( data_source, "OSM" ) == 0 )
          osm_way( ways );
      else
          shp_way( ways );
    }

    void osm_way( const osmium::Way& ways ){
        const osmium::TagList& tags = ways.tags();
        if(!(tags.has_key(key)))
          return;

        if(tags.has_key(junc) && !strcmp(tags.get_value_by_key(junc), "roundabout")){

          SVal seg = createSegmentFromWay(ways);

          if(ways.nodes().front().ref() == ways.nodes().back().ref()){
            segments[dataset].insert({seg->get_id(), seg});
            if(!dataset)
              ref_data << seg->get_id() << "\t" << seg->get_first()->get_id() << "\t" << seg->get_second()->get_id() << "\t" << seg->get_id() << "\t" << seg->get_length() << "\t" << 0 << "\t" << seg->get_buffer()->getArea() << "\t" << "***\n";
            else
              tar_data << seg->get_id() << "\t" << seg->get_first()->get_id() << "\t" << seg->get_second()->get_id() << "\t" << seg->get_id() << "\t" << seg->get_length() << "\t" << 0 << "\t" << seg->get_buffer()->getArea() << "\t" << "***\n";

            return;
          }
          incomplete_roundabouts.push_back(seg);
          return;
        }
        createSegmentsFromWay(ways);
    }

    void shp_way( const osmium::Way& ways ){

      if(ways.nodes().front().ref() == ways.nodes().back().ref()) {

          SVal seg = createSegmentFromWay(ways);

          segments[dataset].insert({seg->get_id(), seg});
          if(!dataset)
            ref_data << seg->get_id() << "\t" << seg->get_first()->get_id() << "\t" << seg->get_second()->get_id() << "\t" << seg->get_id() << "\t" << seg->get_length() << "\t" << 0 << "\t" << seg->get_buffer()->getArea() << "\t" << "***\n";
          else
            tar_data << seg->get_id() << "\t" << seg->get_first()->get_id() << "\t" << seg->get_second()->get_id() << "\t" << seg->get_id() << "\t" << seg->get_length() << "\t" << 0 << "\t" << seg->get_buffer()->getArea() << "\t" << "***\n";

          return;
      }
      createSegmentsFromWay(ways);
    }

    void createSegmentsFromWay(const osmium::Way& ways){

      flag = 0; length = 0; ctr = 0;

      minx = DBL_MAX;
      maxx = -DBL_MAX;
      miny = DBL_MAX;
      maxy = -DBL_MAX;

      cl = new CoordinateArraySequence();

      boost::shared_ptr<Terminal_node> first;
      boost::shared_ptr<Terminal_node> second;

      int size = ways.nodes().size();
      Coordinates c1{ways.nodes().front().location()};

      for(const auto& n : ways.nodes()){

        ctr++;

        Coordinates c_obj{n.location()};
        Coordinates c_mcp = lonlat_to_mercator(c_obj);
        IntNode int_node{n.ref(), n.location()};
        cl->add(Coordinate(pm->makePrecise(c_mcp.x), pm->makePrecise(c_mcp.y)));
        nodeseq.push_back(int_node);

        double x = c_mcp.x;
        double y = c_mcp.y;
        assert(x >= -20036376 && x <= 20050000);
        assert(y >= -20050000 && y <= 20050000);

        minx = std::min(minx, x);
        miny = std::min(miny, y);
        maxx = std::max(maxx, x);
        maxy = std::max(maxy, y);

        osmium::geom::Coordinates c2{n.location()};

        length += osmium::geom::haversine::distance(c1, c2); //haversine distance is in metres.
        c1 = c2;

        auto it = terminal_nodes[dataset].find(n.ref());

        if(ctr == 1){
            first = it->second;
            continue;
        }

        if(ctr == size){
            assert(first);
            second = it->second;

            segmentHelper(ways, first, second);

            first = second;
            minx = maxx = c_mcp.x;
            miny = maxy = c_mcp.y;
            cl->add(Coordinate(pm->makePrecise(c_mcp.x), pm->makePrecise(c_mcp.y)));
            nodeseq.push_back(int_node);
            continue;
        }

        if(it != terminal_nodes[dataset].end() && length > 1){ //"it" is one of the terminal nodes.

            second = it->second;
            segmentHelper(ways, first, second);

            first = second;
            minx = maxx = c_mcp.x;
            miny = maxy = c_mcp.y;
            cl->add(Coordinate(pm->makePrecise(c_mcp.x), pm->makePrecise(c_mcp.y)));
            nodeseq.push_back(int_node);
        }

        else{
            Coordinates c1 = ways.nodes()[ctr-2].location(), c2 = ways.nodes()[ctr-1].location(), c3=c2, c4 = ways.nodes()[ctr].location();

            if(withinAngle(c1, c2, c3, c4) || length <= 1)
                continue;

            boost::shared_ptr<Terminal_node> tnode{new Terminal_node(n.ref(), n.location(), -1)};
            second = tnode;
            segmentHelper(ways, first, second);
            first = second;
            minx = maxx = c_mcp.x;
            miny = maxy = c_mcp.y;
            cl->add(Coordinate(pm->makePrecise(c_mcp.x), pm->makePrecise(c_mcp.y)));
            nodeseq.push_back(int_node);
          }
      }
      delete cl;
      nodeseq.clear();
    }

    void segmentHelper(const osmium::Way& ways, boost::shared_ptr<Terminal_node> first, boost::shared_ptr<Terminal_node> second){
      const Geometry* ls = global_factory->createLineString(cl);
      Geometry* buff = (ls->buffer(ce));

      uint64_t minX = (minx-nnn);
      uint64_t minY = (miny-nnn);
      uint64_t maxX = (maxx-nnn);
      uint64_t maxY = (maxy-nnn);

      std::vector<boost::uint64_t> lower = {minX, minY};
      std::vector<boost::uint64_t> upper = {maxX, maxY};
      boost::shared_ptr<Rectangle> rect{new Rectangle(lower, upper)};

      if( length < 1 )
          length = 1;

      boost::shared_ptr<Segment_creator> ptr;
      if(first->get_location().x < second->get_location().x){
        boost::shared_ptr<Segment_creator> seg(new Segment_creator((const osmium::object_id_type)seg_id, first, second, max(1.0, length), ways.id(), rect, buff, nodeseq, ce, ways.tags()));
        ptr = seg;
      }
      else if(first->get_location().x > second->get_location().x){
        std::reverse(nodeseq.begin(), nodeseq.end());
        boost::shared_ptr<Segment_creator> seg(new Segment_creator((const osmium::object_id_type)seg_id, second, first, max(1.0, length), ways.id(), rect, buff, nodeseq, ce, ways.tags()));
        ptr = seg;
      }
      else if(first->get_location().y < second->get_location().y){
        boost::shared_ptr<Segment_creator> seg(new Segment_creator((const osmium::object_id_type)seg_id, first, second, max(1.0, length), ways.id(), rect, buff, nodeseq, ce, ways.tags()));
        ptr = seg;
      }
      else{
        std::reverse(nodeseq.begin(), nodeseq.end());
        boost::shared_ptr<Segment_creator> seg(new Segment_creator((const osmium::object_id_type)seg_id, second, first, max(1.0, length), ways.id(), rect, buff, nodeseq, ce, ways.tags()));
        ptr = seg;
      }

      segments[dataset].insert({seg_id, ptr});
      if(!dataset)
        ref_data << seg_id << "\t" << ptr->get_first()->get_id() << "\t" << ptr->get_second()->get_id() << "\t" << ways.id() << "\t" << length << "\t" << ls->getLength()<< "\t" << buff->getArea() << "\t" << "\n";
      else
        tar_data << seg_id << "\t" << ptr->get_first()->get_id() << "\t" << ptr->get_second()->get_id() << "\t" << ways.id() << "\t" << length << "\t" << ls->getLength()<< "\t" << buff->getArea() << "\t" << "\n";

      ++seg_id;
      previous = ptr;
      length = 0;

      cl->clear();
      nodeseq.clear();
    }

    bool withinAngle(Coordinates c1, Coordinates c2, Coordinates c3, Coordinates c4){
        double num1 = c2.y - c1.y;
    		double deno1 = c2.x - c1.x;
    		double num2 = c3.y - c4.y;
    		double deno2 = c3.x - c4.x;
        double val1=0, val2=0;
        val1 = num1/deno1;
        val2 = num2/deno2;
        double atan1 = atan(val1);
        double atan2 = atan(val2);
        if(deno1 == 0)
          atan1 = PI/2;
        if(deno2 == 0)
          atan2 = PI/2;

    		if(min(abs(atan1 - atan2), PI-abs(atan1 - atan2)) <= PI/3)
    			return true;
    		return false;
    }

    SVal concatR(SVal base, SVal other){
        osmium::object_id_type id = seg_id;
        ++seg_id;
        boost::shared_ptr<Terminal_node> first = base->get_first();
        boost::shared_ptr<Terminal_node> second = other->get_second();
        double length = base->get_length() + other->get_length();
        osmium::object_id_type way_id = base->get_way_id();
        boost::shared_ptr<Rectangle> box = concatMBRR(base->get_box(), other->get_box());
        vector<IntNode> cl = concatIntermediateNodesR(base->get_intermediate_nodes(), other->get_intermediate_nodes());

        if(first->get_id() == second->get_id()){
            Geometry* buff = concatBuffR(base, other, true);
            boost::shared_ptr<Segment_creator> joined(new Segment_creator(id, first, second, max(1.0, length), way_id, box, buff, cl, ce, true, base->tags()));
            return joined;
        }
        Geometry* buff = concatBuffR(base, other, false);
        boost::shared_ptr<Segment_creator> joined(new Segment_creator(id, first, second, max(1.0, length), way_id, box, buff, cl, ce, false, base->tags()));
        return joined;
    }

    void empty_roundabouts(){
        if(incomplete_roundabouts.size() == 0)
            return;

        // cout << "All ids ... \n";
        // for( auto it : incomplete_roundabouts )
        //     cout << it->get_way_id() << "\n";
        // cout << " ------------------- \n";

        for(auto it1 = incomplete_roundabouts.begin(); it1 != incomplete_roundabouts.end(); it1++){
            // cout << "id of roundabout ... " << ( *it1 )->get_way_id() << "\n";
            for(auto it2 = it1+1; it2 != incomplete_roundabouts.end(); it2++){
                SVal seg1 = *it1, seg2 = *it2;

                if(seg1->get_first()->get_id() == seg2->get_second()->get_id()){
                    SVal joined = concatR(seg2, seg1);
                    incomplete_roundabouts.erase(it2--);
                    incomplete_roundabouts.erase(it1--);
                    if(joined->is_roundabout()){
                      for(auto& node : joined->get_intermediate_nodes()){
                        if(terminal_nodes[dataset].find(node.get_id()) != terminal_nodes[dataset].end())
                          joined->insert_tn(terminal_nodes[dataset].find(node.get_id())->second);
                      }
                      segments[dataset].insert({joined->get_id(), joined});
                      if(!dataset)
                        ref_data << joined->get_id() << "\t" << joined->get_first()->get_id() << "\t" << joined->get_second()->get_id() << "\t" << joined->get_way_id() << "\t" << joined->get_length() << "\t" << 0 << "\t" << joined->get_buffer()->getArea() << "\t" << "***\n";
                      else
                        tar_data << joined->get_id() << "\t" << joined->get_first()->get_id() << "\t" << joined->get_second()->get_id() << "\t" << joined->get_way_id() << "\t" << joined->get_length() << "\t" << 0 << "\t" << joined->get_buffer()->getArea() << "\t" << "***\n";
                    }
                    else
                      incomplete_roundabouts.push_back(joined);
                    break;
                }
                if(seg1->get_second()->get_id() == seg2->get_first()->get_id()){
                    SVal joined = concatR(seg1, seg2);
                    incomplete_roundabouts.erase(it2--);
                    incomplete_roundabouts.erase(it1--);
                    if(joined->is_roundabout()){
                      for(auto& node : joined->get_intermediate_nodes()){
                        if(terminal_nodes[dataset].find(node.get_id()) != terminal_nodes[dataset].end())
                          joined->insert_tn(terminal_nodes[dataset].find(node.get_id())->second);
                      }
                      segments[dataset].insert({joined->get_id(), joined});
                      if(!dataset)
                        ref_data << joined->get_id() << "\t" << joined->get_first()->get_id() << "\t" << joined->get_second()->get_id() << "\t" << joined->get_way_id() << "\t" << joined->get_length() << "\t" << 0 << "\t" << joined->get_buffer()->getArea() << "\t" << "***\n";
                      else
                        tar_data << joined->get_id() << "\t" << joined->get_first()->get_id() << "\t" << joined->get_second()->get_id() << "\t" << joined->get_way_id() << "\t" << joined->get_length() << "\t" << 0 << "\t" << joined->get_buffer()->getArea() << "\t" << "***\n";
                    }
                    else
                      incomplete_roundabouts.push_back(joined);
                    break;
                }
            }
        }
        for( auto it : incomplete_roundabouts ){
          segments[dataset].insert({it->get_id(), it});
          if(!dataset)
            ref_data << it->get_id() << "\t" << it->get_first()->get_id() << "\t" << it->get_second()->get_id() << "\t" << it->get_way_id() << "\t" << it->get_length() << "\t" << 0 << "\t" << it->get_buffer()->getArea() << "\t" << "*\n";
          else
            tar_data << it->get_id() << "\t" << it->get_first()->get_id() << "\t" << it->get_second()->get_id() << "\t" << it->get_way_id() << "\t" << it->get_length() << "\t" << 0 << "\t" << it->get_buffer()->getArea() << "\t" << "*\n";
        }
        // assert(incomplete_roundabouts.size() == 0);
    }
};


#endif
