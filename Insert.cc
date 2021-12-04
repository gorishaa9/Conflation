#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <boost/shared_ptr.hpp>
#include <boost/timer/timer.hpp>
#include <osmium/osm/location.hpp>
#include <osmium/osm/box.hpp>
#include <fstream>
#include <limits>

#include "Rectangle.hh"
#include "RTree.hh"
#include "Insert.hh"
#include "../map_read/segment.hpp"


#define CE 1000

using boost::shared_ptr;
using std::vector;
using boost::uint64_t;
using SKey = long long int;
using SVal = boost::shared_ptr<Segment_creator>;

ofstream oFile("MBR.txt");

RTree tree{};
uint32_t nn = std::numeric_limits<int>::min();

std::vector<boost::shared_ptr<Rectangle>> rectangles;
std::vector<long long int> ids;

using namespace std;

std::list<long long int> rtree_search(boost::shared_ptr<Segment_creator> query){

    uint64_t loc1x = query->get_box()->getLower()[0];
    uint64_t loc1y = query->get_box()->getLower()[1];
    uint64_t loc2x = query->get_box()->getUpper()[0];
    uint64_t loc2y = query->get_box()->getUpper()[1];
    std::vector<boost::uint64_t> min;
    std::vector<boost::uint64_t> max;
    min.push_back(loc1x-CE);
    min.push_back(loc1y-CE);
    max.push_back(loc2x+CE);
    max.push_back(loc2y+CE);

    boost::shared_ptr<Rectangle> search_rect(new Rectangle(min, max));

    auto begin = std::chrono::high_resolution_clock::now();
    std::list<boost::shared_ptr<NodeEntry>> result = tree.search(search_rect);
    auto end = std::chrono::high_resolution_clock::now();

    std::list<long long int> query_result;

    if(result.empty())
        cout << "No object found! \n";
    else{
        for(auto& it : result){
            long long int id = boost::dynamic_pointer_cast<LeafEntry>(it)->getId();
            query_result.insert(query_result.end(), id);
        }
    }

    std::cout << "Running  query takes "
              << (double)std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "ms" << std::endl;

    return query_result;
}

int rtree_insert(std::map<SKey, SVal> segments)
{

    for(const auto it : segments){

        uint64_t loc1x = it.second->get_box()->getLower()[0];
        uint64_t loc1y = it.second->get_box()->getLower()[1];
        uint64_t loc2x = it.second->get_box()->getUpper()[0];
        uint64_t loc2y = it.second->get_box()->getUpper()[1];

        std::vector<boost::uint64_t> min;
        std::vector<boost::uint64_t> max;
        min.push_back(loc1x-CE);
        min.push_back(loc1y-CE);
        max.push_back(loc2x+CE);
        max.push_back(loc2y+CE);


        boost::shared_ptr<Rectangle> r(new Rectangle(min, max));
        rectangles.push_back(r);
        ids.push_back(it.first);
    }
    boost::timer::cpu_timer timer;

    auto begin = std::chrono::high_resolution_clock::now();

    for (size_t j = 0; j < rectangles.size(); ++j)
    {
        tree.insert(rectangles[j], ids[j]);
    }
    timer.stop();

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Inserting " << rectangles.size() << " segments: "
              << (double)std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "ms" << std::endl;

    return 0;
}
