#ifndef SEGMENT_H
#define SEGMENT_H

#include <osmium/osm/object.hpp>
#include <osmium/osm/way.hpp>
#include <osmium/osm/tag.hpp>
#include <osmium/osm/node_ref.hpp>
#include <boost/shared_ptr.hpp>
#include <geos/geom/Geometry.h>
#include <geos/geom/CoordinateArraySequence.h>

#include "terminal_node.hpp"
#include "Rectangle.hh"
#include "IntNode.hpp"

using namespace std;
using namespace geos::geom;
using Geometry = geos::geom::Geometry;

class Segment_creator{
private:

  osmium::object_id_type _id;
  boost::shared_ptr<Terminal_node> _first;
  boost::shared_ptr<Terminal_node> _second;
  double _length;
  const osmium::object_id_type _wayid;
  boost::shared_ptr<Rectangle> _box;
  Geometry* _buff;
  vector<IntNode> _cl;
  uint32_t _ce;
  bool _roundabout;
  vector<long long int> _parents = {};
  vector<boost::shared_ptr<Terminal_node>> _ids = {};
  const osmium::TagList& _tags;

public:
  Segment_creator(osmium::object_id_type id,
                  boost::shared_ptr<Terminal_node> first,
                  boost::shared_ptr<Terminal_node> second,
                  double length,
                  const osmium::object_id_type way_id,
                  boost::shared_ptr<Rectangle> box,
                  Geometry* buff,
                  vector<IntNode> cl,
                  uint32_t ce,
                  const osmium::TagList& tags) :
      _id(id),
      _first(first),
      _second(second),
      _length(length),
      _wayid(way_id),
      _box(box),
      _buff(buff),
      _cl(cl),
      _ce(ce),
      _roundabout(false),
      _tags(tags) {
      }

    Segment_creator(osmium::object_id_type id,
                    boost::shared_ptr<Terminal_node> first,
                    boost::shared_ptr<Terminal_node> second,
                    double length,
                    const osmium::object_id_type way_id,
                    boost::shared_ptr<Rectangle> box,
                    Geometry* buff,
                    vector<IntNode> cl,
                    uint32_t ce,
                    bool roundabout,
                    const osmium::TagList& tags) :
          _id(id),
          _first(first),
          _second(second),
          _length(length),
          _wayid(way_id),
          _box(box),
          _buff(buff),
          _cl(cl),
          _ce(ce),
          _roundabout(roundabout),
          _tags(tags) {
          }

    osmium::object_id_type get_id(){
      return _id;
    }

    void set_id(osmium::object_id_type id){
      _id = id;
    }

    boost::shared_ptr<Terminal_node> get_first(){
      return _first;
    }

    boost::shared_ptr<Terminal_node> get_second(){
      return _second;
    }

    double get_length(){
      return _length;
    }

    const osmium::object_id_type get_way_id(){
      return _wayid;
    }

    boost::shared_ptr<Rectangle> get_box(){
      return _box;
    }

    Geometry* get_buffer(){
      return _buff;
    }

    vector<IntNode> get_intermediate_nodes(){
      return _cl;
    }

    uint32_t getCircularError(){
      return _ce;
    }

    bool is_roundabout(){
      return _roundabout;
    }

    vector<long long int> parents(){ //for self join, already joined segments.
      return _parents;
    }

    void insert_parent(long long int seg){
      _parents.push_back(seg);
    }

    void insert_parent(vector<long long int> segs){
      _parents.insert(_parents.end(), segs.begin(), segs.end());
    }

    void insert_tn(boost::shared_ptr<Terminal_node> id){ //to keep the reference id of the terminal nodes in the roundabout to determine the join.
      _ids.push_back(id);
    }

    void insert_tns(vector<boost::shared_ptr<Terminal_node>> ids){
      _ids.insert(_ids.end(), ids.begin(), ids.end());
    }

    vector<boost::shared_ptr<Terminal_node>> get_tns(){
      return _ids;
    }

    const osmium::TagList& tags(){
      return _tags;
    }
};

#endif
