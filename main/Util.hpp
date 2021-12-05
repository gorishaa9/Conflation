#ifndef UTIL_H
#define UTIL_H

#include <math.h>
#include <stdlib.h>
#include <cmath>
#include <geos/geom/Geometry.h>

#include "terminal_node.hpp"
#include "segment.hpp"
#include "Tuple.hpp"
#include "Relation.hpp"
#include "IntNode.hpp"
#include "SegmentHandler.hpp"

#define ce 5
#define RADIUS 6378137
#define PI 3.14159265
#define RadianCoeff PI/180
#define DegreesCoeff 180/PI


using namespace osmium::geom;
using namespace geos::geom;
using Coordinate = geos::geom::Coordinate;

using TVal = boost::shared_ptr<Tuple>;

ofstream check_coord("check_coord.txt");

map<SKey, vector<SKey>> hashes; //the key is a terminal node: the value is the list od the id of segments connected to this terminal node.

struct JoinDegree{
    SKey key;
    short degree;
};


ofstream ftest("ftest.txt");
ofstream adjlist1("adjlist1.txt");
ofstream adjlist2("adjlist2.txt");


bool isWithinCE(boost::shared_ptr<Terminal_node> n1, boost::shared_ptr<Terminal_node> n2, double CE){
    Coordinates c1 = n1->get_location();
    Coordinates c2 = n2->get_location();
    double dist = haversine::distance(c1, c2);
    if(dist <= CE)
        return true;
    else
        return false;
}


double similarityScore(boost::shared_ptr<Segment_creator> seg1, boost::shared_ptr<Segment_creator> seg2){
    Geometry* g1 = seg1->get_buffer();
    Geometry* g2 = seg2->get_buffer();
    Geometry* g = g1->intersection(g2);
    Geometry* u = g1->Union(g2);
    double score = g->getArea()/u->getArea();
    return score;
}


boost::shared_ptr<Rectangle> concatMBR(boost::shared_ptr<Rectangle> box1, boost::shared_ptr<Rectangle> box2){
    uint64_t minX = std::min(box1->getLower()[0], box2->getLower()[0]);
    uint64_t minY = std::min(box1->getLower()[1], box2->getLower()[1]);
    uint64_t maxX = std::min(box1->getUpper()[0], box2->getUpper()[0]);
    uint64_t maxY = std::min(box1->getUpper()[1], box2->getUpper()[1]);

    std::vector<boost::uint64_t> lower = {minX, minY};
    std::vector<boost::uint64_t> upper = {maxX, maxY};
    boost::shared_ptr<Rectangle> box{new Rectangle(lower, upper)};

    return box;
}


vector<IntNode> concatIntermediateNodes(vector<IntNode> cl1, vector<IntNode> cl2){
    vector<IntNode> nodeseq;
    nodeseq.insert(nodeseq.end(), cl1.begin(), cl1.end());
    nodeseq.insert(nodeseq.end(), cl2.begin()+1, cl2.end());
    return nodeseq;
}


Geometry* concatBuff(vector<IntNode> nodeseq){
    // vector<Coordinates> nodeseq = concatIntermediateNodes(seg1->get_intermediate_nodes(), seg2->get_intermediate_nodes());
    CoordinateArraySequence* cl = new CoordinateArraySequence();
    for(auto& int_node : nodeseq){
        Coordinates it = int_node.get_location();
        Coordinates mer = lonlat_to_mercator(Coordinates(it.x, it.y));
        cl->add(Coordinate(pm->makePrecise(mer.x), pm->makePrecise(mer.y)));
        //delete mer;
    }
    Geometry* ls = global_factory->createLineString(cl);
    Geometry* buff = ls->buffer(ce);
    return buff;
}


inline double distance(Coordinates c1, Coordinates c2){
    return osmium::geom::haversine::distance(c1, c2);
}


Coordinates getInterpolatedNodeCoordinates(Coordinates c, Coordinates c1, Coordinates c2){
    Coordinates c_m = lonlat_to_mercator(c);
    Coordinates c1_m = lonlat_to_mercator(c1);
    Coordinates c2_m = lonlat_to_mercator(c2);
    if(c1.x == c2.x)
        return Coordinates{c1.x, c.y};
    // double m_prime = (c2_m.y - c1_m.y)/(c2_m.x - c1_m.x);
    double x, y;
    //t=[(Cx-Ax)(Bx-Ax)+(Cy-Ay)(By-Ay)]/[(Bx-Ax)^2+(By-Ay)^2]
    double t = ( ( c_m.x - c1_m.x )*( c2_m.x - c1_m.x ) + ( c_m.y - c1_m.y )*( c2_m.y - c1_m.y ) )/( (c2_m.x - c1_m.x)*(c2_m.x - c1_m.x) + (c2_m.y - c1_m.y)*(c2_m.y - c1_m.y) );

    // x = (c_m.x + (m_prime * m_prime * c2_m.x) + (m_prime * (c_m.y - c2_m.y)))/(1+(m_prime * m_prime));
    // y = (m_prime != 0)? m_prime*(x - c1_m.x) + c1_m.y : c1_m.y;

    x = c1_m.x + ( t * (c2_m.x - c1_m.x) );
    y = c1_m.y + ( t * (c2_m.y - c1_m.y) );

    Coordinates ci{x, y};


    Coordinates ci_m = mercator_to_lonlat(ci);
    assert(ci.x >= -20050000 && ci.x <= 20050000);
    assert(ci.y >= -20050000 && ci.y <= 20050000);
    return ci_m;
}


inline bool liesOnSegment(Coordinates c, Coordinates c1, Coordinates c2){   // checks if c lies on the edge with the endpoints c1 and c2.
    bool flag1 = false, flag2 = false;
    if(c1.x <= c2.x)
        flag1 = c.x >= c1.x && c.x <= c2.x;
    else
        flag1 = c.x >= c2.x && c.x <= c1.x;
    if(c1.y <= c2.y)
        flag2 = c.y >= c1.y && c.y <= c2.y;
    else
        flag2 = c.y >= c2.y && c.y <= c1.y;
    return flag1 && flag2;
}


inline bool liesToLeft(Coordinates c, Coordinates c1, Coordinates c2){
    //bool flag1= false, flag2 = false;
    return (c.x <= c1.x && c.x <= c2.x);
}


double calculateLength(vector<IntNode> nodeseq){
    Coordinates c1 = nodeseq.front().get_location();
    Coordinates c2;
    int length = 0;
    for(auto& int_node : nodeseq){
        Coordinates c2 = int_node.get_location();
        length += osmium::geom::haversine::distance(c1, c2);
        c1 = c2;
    }
    return length;
}


boost::shared_ptr<Rectangle> calculateBox(vector<IntNode> split_segment){
    int32_t nnn = std::numeric_limits<int>::min();
    double minx = DBL_MAX;
    double maxx = -DBL_MAX;
    double miny = DBL_MAX;
    double maxy = -DBL_MAX;
    for(auto& int_node : split_segment){
        Coordinates i = int_node.get_location();
        Coordinates i_m = lonlat_to_mercator(i);
        double x = i_m.x;
        double y = i_m.y;
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
    return rect;
}


Geometry* calculateBuffer(vector<IntNode> segment){
    CoordinateArraySequence* cl = new CoordinateArraySequence();
    for(auto& int_node : segment){
        Coordinates i = int_node.get_location();
        Coordinates mer = lonlat_to_mercator(i);
        assert(mer.x >= -20036376 && mer.x <= 20050000);
        assert(mer.y >= -20050000 && mer.y <= 20050000);
        cl->add(Coordinate(pm->makePrecise(mer.x), pm->makePrecise(mer.y)));
    }
    Geometry* ls = global_factory->createLineString(cl);
    return ls->buffer(ce);
}

// seg1 is within bounds of seg2
bool isWithinBounds(boost::shared_ptr<Segment_creator> seg1, boost::shared_ptr<Segment_creator> seg2){
    Coordinates c_mcp = lonlat_to_mercator(seg1->get_intermediate_nodes()[0].get_location());
    if( !seg2->get_buffer()->intersects( global_factory->createPoint(Coordinate(pm->makePrecise(c_mcp.x), pm->makePrecise(c_mcp.y))) ) )
        return false;
    c_mcp = lonlat_to_mercator(seg1->get_intermediate_nodes().back().get_location());
    if( !seg2->get_buffer()->intersects( global_factory->createPoint(Coordinate(pm->makePrecise(c_mcp.x), pm->makePrecise(c_mcp.y))) ) )
        return false;
    return true;
}

//seg1 is after seg2 --> only the second point of seg1 is within the buffer of seg2.
bool isAfter(SVal seg1, SVal seg2){
    Coordinates c_mcp_s = lonlat_to_mercator(seg1->get_intermediate_nodes().front().get_location());
    Coordinates c_mcp_b = lonlat_to_mercator(seg1->get_intermediate_nodes().back().get_location());
    if( seg2->get_buffer()->intersects( global_factory->createPoint(Coordinate(pm->makePrecise(c_mcp_s.x), pm->makePrecise(c_mcp_s.y))) )
     && !seg2->get_buffer()->intersects( global_factory->createPoint(Coordinate(pm->makePrecise(c_mcp_b.x), pm->makePrecise(c_mcp_b.y))) ) )
        return true;
    return false;
}

bool pointWithinBuffer( boost::shared_ptr<Terminal_node> point, SVal seg ){
    Coordinates c_mcp = lonlat_to_mercator( point->get_location() );
    if( seg->get_buffer()->intersects( global_factory->createPoint(Coordinate(pm->makePrecise(c_mcp.x), pm->makePrecise(c_mcp.y)))) )
        return true;
    return false;
}


map<SKey, vector<SKey>> findJoins(map<SKey, SVal> segments, set<SKey> single_matches){
    //find hash by their end points.
    vector<SKey> roundabout_ids;
    for(auto& it : segments){
        // if(it.first == 1859)
        //     cout << "1859 : " << it.second->is_roundabout() << "\n";
        if(single_matches.find(it.first) != single_matches.end())
            continue;
        if(it.second->is_roundabout()){
            if(it.first == 1859)
                cout << "1859 is a roundabout.\n";
            roundabout_ids.push_back(it.first);
            for(auto& tnode : it.second->get_tns()){
                if(it.first == 1859)
                    cout << "1859: " << tnode->get_id() << "\n";
                // if(it.first == 1859)
                //     cout << "1859: " << tnode->get_id() << "\n";
                if(hashes.find(tnode->get_id()) != hashes.end()){
                    if( find( hashes[ tnode->get_id() ].begin(), hashes[ tnode->get_id() ].end(), it.first ) == hashes[ tnode->get_id() ].end() )
                        hashes.find(tnode->get_id())->second.push_back(it.first);
                }
                else
                    hashes.insert({tnode->get_id(), {it.first}});
            }
            continue;
        }
        SKey id1 = it.second->get_first()->get_id();
        SKey id2 = it.second->get_second()->get_id();
        if(hashes.find(id1) != hashes.end())
            hashes[ id1 ].push_back(it.first);
        else
            hashes.insert({id1, {it.first}});

        if(hashes.find(id2) != hashes.end())
            hashes[ id2 ].push_back(it.first);
        else
            hashes.insert({id2, {it.first}});
    }

    map<SKey, vector<SKey>> joins; //the key is the id of the segment, the value is the list of ids of the other segments it is connected to.

    for(auto& it : segments){
        SKey id = it.first;

        if(single_matches.find(id) != single_matches.end())
            continue;
        if(segments.find(id)->second->is_roundabout()){
            vector<SKey> list;
            for(auto& tnode : it.second->get_tns()){
                vector<SKey> list_joins = hashes.find(tnode->get_id())->second;
                list.insert(list.end(), list_joins.begin(), list_joins.end());
            }
            for(auto item = list.begin(); item != list.end(); item++){
                if(*item == id)
                    list.erase(item--);
            }

            vector< SKey > list_copy = list; //list apparently has duplicate items. Remove them.
            list.clear();
            for( auto l_item : list_copy ){
                if( find( list.begin(), list.end(), l_item ) == list.end() )
                    list.push_back( l_item );
            }
            joins.insert({id, list});
            continue;
        }
        SKey id1 = it.second->get_first()->get_id();
        SKey id2 = it.second->get_second()->get_id();

        if((hashes.find(id1) == hashes.end() || hashes.find(id2) == hashes.end())){
            cout << "Here!" << "\n";
            exit(EXIT_FAILURE);
        }

        vector<SKey> list1 = hashes.find(id1)->second;
        vector<SKey> list2 = hashes.find(id2)->second;
        // if(id == 1312){
        //     cout << "All segments to id1 : " << id1 << "\n";
        //     for(auto ele : list1)
        //         cout << ele << "\n";
        //
        //     cout << "All segments to id2 : " << id2 << "\n";
        //     for(auto ele : list2)
        //         cout << ele << "\n";
        // }

        bool has_roundabout = false;
        vector<SKey> rid;
        for(auto id_seg : list1){
          if( id == 1860 )
              cout << "1860 : " << id_seg << "\n";
            if(segments.find(id_seg)->second->is_roundabout()){
                has_roundabout = true;
                rid.push_back(id_seg);
            }
        }

        if(has_roundabout)  //if the segment is joined to roundabout on one end, only join with roundabout.
            list1 = rid;
        has_roundabout = false;

        rid.clear();

        for(auto id_seg : list2){
            if( id == 1860 )
                cout << "1860 : " << id_seg << "\n";
            if(segments.find(id_seg)->second->is_roundabout()){
                has_roundabout = true;
                rid.push_back(id_seg);
            }
        }
        if(has_roundabout)  //if the segment is joined to roundabout on one end, only join with roundabout.
            list2 = rid;

        for(auto iter = list1.begin(); iter != list1.end(); iter++){
            if(*iter == id)
                list1.erase(iter--);
        }
        for(auto iter = list2.begin(); iter != list2.end(); iter++){
            if(*iter == id || find( list1.begin(), list1.end(), *iter ) != list1.end() )
                list2.erase(iter--);
        }

        vector<SKey> list = list1;
        list.insert(list.end(), list2.begin(), list2.end());

        // list.erase( unique( list.begin(), list.end() ), list.end() );

        joins.insert({id, list});
    }
    return joins;
}


map<SKey, double> calculateLengths(vector<SKey> join_order, map<SKey, SVal> segments){
    map<SKey, double> lengths;
    for(auto it : join_order){
        //SKey id = it.first;
        SVal segment = segments.find(it)->second;
        double length = segment->get_length();
        lengths.insert({it, length});
    }
    return lengths;
}


bool map_compare(map<SKey, TVal> &map1, map<SKey, TVal> &map2){
    if(map1.size() != map2.size())
        return false;
    for(auto& it1 : map1){
        SKey id1 = it1.first;
        TVal tuple1 = it1.second;
        if(map2.find(id1) == map2.end())
            return false;
        if(map2.find(id1)->second != tuple1)
            return false;
    }
    return true;
}


inline bool sameLocation(Coordinates c1, Coordinates c2){
    return c1.x == c2.x && c1.y == c2.y;
}


bool samePattern(SVal delta, SVal cs, Coordinates c1, Coordinates c2, Coordinates cs1, Coordinates cs2){
    if(abs(c1.x - c2.x) > abs(c1.y - c2.y)){
        if(delta->get_id() == 1000 && cs->get_id() == 2911){
            cout << "\nHere!\t";
            cout << "X\n";
        }
        if(c1.x < c2.x && cs1.x < cs2.x)
            return true;
        else if(c1.x > c2.x && cs1.x > cs2.x)
            return true;
    }
    else{
        if(c1.y < c2.y && cs1.y < cs2.y)
            return true;
        else if(c1.y > c2.y && cs1.y > cs2.y)
            return true;
    }
    return false;
}


CoordinateArraySequence* formCA(vector<IntNode> cl){
    CoordinateArraySequence* cseq = new CoordinateArraySequence();
    for(auto& int_node : cl){
        Coordinates c_obj = int_node.get_location();
        Coordinates c_mcp = lonlat_to_mercator(c_obj);
        cseq->add(Coordinate(pm->makePrecise(c_mcp.x), pm->makePrecise(c_mcp.y)));
    }
    return cseq;
}


inline Coordinates degreeToRadians(Coordinates angle){
    return Coordinates{angle.x*PI/180, angle.y*PI/180};
}


inline Coordinates radiansToDegree(Coordinates angle){
    return Coordinates{angle.x*180/PI, angle.y*180/PI};
}

#endif
