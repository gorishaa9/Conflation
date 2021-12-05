#ifndef SEGMENT_UTILS_H
#define SEGMENT_UTILS_H

#include <osmium/handler.hpp>
#include <boost/shared_ptr.hpp>
#include <osmium/osm/object.hpp>
#include <osmium/osm/location.hpp>
#include <osmium/osm/tag.hpp>
#include <osmium/geom/haversine.hpp>
#include <osmium/geom/coordinates.hpp>
#include <osmium/geom/mercator_projection.hpp>

#include <geos/geom/PrecisionModel.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Point.h>
#include <geos/geom/LineString.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/GeometryCollection.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/CoordinateArraySequence.h>
#include <geos/util/GEOSException.h>
#include <geos/util/IllegalArgumentException.h>

#include <iostream>
#include <float.h>
#include <limits>
#include <math.h>
#include <fstream>
#include <chrono>
#include <thread>

#include "terminal_node.hpp"
#include "segment.hpp"
#include "CountHandler.hpp"
#include "IntNode.hpp"

#define ce 5

using namespace std;
using namespace std::this_thread; // sleep_for, sleep_until
using namespace std::chrono; // nanoseconds, system_clock, seconds

using namespace geos;
using namespace geos::geom;
using namespace osmium::geom;

using SKey = long long int;
using SVal = boost::shared_ptr<Segment_creator>;

SKey seg_id = 1000;
extern vector<map<SKey, SVal>> segments;
map<SKey, SKey> relaxed_joins;
vector<SVal> incomplete_roundabouts;
const char key[] = "highway";
const char junc[] = "junction";
int32_t nnn = std::numeric_limits<int>::min();
//double scale = 100000;

extern PrecisionModel* pm = new PrecisionModel();
extern auto global_factory = GeometryFactory::create(pm, 3857);

SVal createSegmentFromWay(const osmium::Way& ways){
    double length = 0;

    boost::shared_ptr<Terminal_node> second{new Terminal_node(ways.nodes().back().ref(), ways.nodes().back().location(), -1)};
    boost::shared_ptr<Terminal_node> first{new Terminal_node(ways.nodes().front().ref(), ways.nodes().front().location(), -1)};

    double minx = DBL_MAX;
    double maxx = -DBL_MAX;
    double miny = DBL_MAX;
    double maxy = -DBL_MAX;
    int size = ways.nodes().size();
    Coordinates c1{ways.nodes().front().location()};

    CoordinateArraySequence* cl = new CoordinateArraySequence();
    vector<IntNode> nodeseq;
    vector<boost::shared_ptr<Terminal_node>> ids;

    for(const auto& n : ways.nodes()){
        Coordinates c_obj{n.location()};
        Coordinates c_mcp = lonlat_to_mercator(c_obj);
        IntNode int_node{n.ref(), n.location()};
        cl->add(Coordinate(pm->makePrecise(c_mcp.x), pm->makePrecise(c_mcp.y)));
        nodeseq.push_back(int_node);

        if(terminal_nodes[dataset].find(n.ref()) != terminal_nodes[dataset].end())
            ids.push_back(terminal_nodes[dataset].find(n.ref())->second);

        length += osmium::geom::haversine::distance(c1, c_obj);
        c1 = c_obj;

        double x = c_mcp.x;
        double y = c_mcp.y;
        assert(x >= -20036376 && x <= 20050000);
        assert(y >= -20050000 && y <= 20050000);

        minx = std::min(minx, x);
        miny = std::min(miny, y);
        maxx = std::max(maxx, x);
        maxy = std::max(maxy, y);
    }

    uint64_t minX = (minx-nnn);
    uint64_t minY = (miny-nnn);
    uint64_t maxX = (maxx-nnn);
    uint64_t maxY = (maxy-nnn);

    std::vector<boost::uint64_t> lower = {minX, minY};
    std::vector<boost::uint64_t> upper = {maxX, maxY};
    boost::shared_ptr<Rectangle> rect{new Rectangle(lower, upper)};

    // if(ways.nodes().front().ref() == ways.nodes().back().ref()){
    //     const Geometry* ls = global_factory->createLinearRing(cl);
    //     Geometry* buff = (ls->buffer(ce));
    //     SVal seg(new Segment_creator((const osmium::object_id_type)seg_id, first, second, length, ways.id(), rect, buff, nodeseq, ce, true, ways.tags()));
    //     seg->insert_tns(ids);
    //     nodeseq.clear();
    //     ++seg_id;
    //     return seg;
    // }
    const Geometry* ls = global_factory->createLineString(cl);
    Geometry* buff = (ls->buffer(ce));
    SVal seg(new Segment_creator((const osmium::object_id_type)seg_id, first, second, max(1.0, length), ways.id(), rect, buff, nodeseq, ce, false, ways.tags()));
    nodeseq.clear();
    ++seg_id;
    return seg;
}

boost::shared_ptr<Rectangle> concatMBRR(boost::shared_ptr<Rectangle> box1, boost::shared_ptr<Rectangle> box2){
    uint64_t minX = std::min(box1->getLower()[0], box2->getLower()[0]);
    uint64_t minY = std::min(box1->getLower()[1], box2->getLower()[1]);
    uint64_t maxX = std::min(box1->getUpper()[0], box2->getUpper()[0]);
    uint64_t maxY = std::min(box1->getUpper()[1], box2->getUpper()[1]);

    std::vector<boost::uint64_t> lower = {minX, minY};
    std::vector<boost::uint64_t> upper = {maxX, maxY};
    boost::shared_ptr<Rectangle> box{new Rectangle(lower, upper)};

    return box;
}

vector<IntNode> concatIntermediateNodesR(vector<IntNode> cl1, vector<IntNode> cl2){
    vector<IntNode> nodeseq;
    nodeseq.insert(nodeseq.end(), cl1.begin(), cl1.end());
    nodeseq.insert(nodeseq.end(), cl2.begin()+1, cl2.end());
    return nodeseq;
}

Geometry* concatBuffR(SVal seg1, SVal seg2, bool flag){
    vector<IntNode> nodeseq = concatIntermediateNodesR(seg1->get_intermediate_nodes(), seg2->get_intermediate_nodes());
    CoordinateArraySequence* cl = new CoordinateArraySequence();
    for(auto& int_node : nodeseq){
        Coordinates it = int_node.get_location();
        Coordinates mer = lonlat_to_mercator(Coordinates(it.x, it.y));
        cl->add(Coordinate(pm->makePrecise(mer.x), pm->makePrecise(mer.y)));
        //delete mer;
    }
    if(flag){
        Geometry* ls = global_factory->createLinearRing(cl);
        Geometry* buff = ls->buffer(ce);
        return buff;
    }
    Geometry* ls = global_factory->createLineString(cl);
    Geometry* buff = ls->buffer(ce);
    return buff;
}


#endif
