#ifndef RTREE_INSERT_HH
#define RTREE_INSERT_HH

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <boost/shared_ptr.hpp>
#include <boost/timer/timer.hpp>
#include <osmium/osm/location.hpp>
#include <osmium/osm/box.hpp>

#include "Rectangle.hh"
#include "RTree.hh"
#include "../map_read/segment.hpp"

using SKey = long long int;
using SVal = boost::shared_ptr<Segment_creator>;

extern RTree tree;
extern std::vector<boost::shared_ptr<Rectangle>> rectangles;
extern std::vector<long long int> ids;

std::list<long long int> rtree_search(boost::shared_ptr<Segment_creator>);
int rtree_insert(std::map<SKey, SVal>);
//int insert_check();

#endif
