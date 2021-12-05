#ifndef INT_NODE
#define INT_NODE

#include <osmium/osm/object.hpp>
#include <osmium/osm/location.hpp>
#include <osmium/geom/coordinates.hpp>

using namespace std;
using Coordinates = osmium::geom::Coordinates;

class IntNode{

    osmium::object_id_type _id;
    Coordinates _location;

public:
    IntNode (){}

    IntNode (const osmium::object_id_type id,
             const Coordinates location) :

        _id(id),
        _location(location){
        }

    const osmium::object_id_type get_id(){
        return _id;
    }

    const Coordinates get_location(){
        return _location;
    }
};

#endif
