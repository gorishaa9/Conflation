#ifndef TERMINAL_NODE_H
#define TERMINAL_NODE_H

#include <osmium/osm/object.hpp>
#include <osmium/osm/location.hpp>
#include <osmium/geom/coordinates.hpp>

using namespace std;
using Coordinates = osmium::geom::Coordinates;

class Terminal_node  {

  const osmium::object_id_type _id;
  const Coordinates _location;
  int _degree;

public:

  Terminal_node (const osmium::object_id_type id,
                 const Coordinates location,
                 int degree)  :
      _id(id),
      _location(location),
      _degree(degree) {
      }

    const osmium::object_id_type get_id(){
      return _id;
    }

    const Coordinates get_location(){
      return _location;
    }

    int get_degree(){
      return _degree;
    }

    void increment_degree(){
      ++_degree;
    }
};

#endif
