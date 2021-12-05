#ifndef MAP_WRITE_H
#define MAP_WRITE_H

#include <utility>
#include <cmath>
#include <osmium/io/any_output.hpp>
#include <osmium/builder/osm_object_builder.hpp>
#include <osmium/osm/node_ref.hpp>
#include <osmium/builder/attr.hpp>

#include "RJReGroup_wprints.hpp"

using namespace osmium::memory;
using namespace osmium::io;
using namespace osmium::builder;
using namespace osmium::builder::attr; // NOLINT(google-build-using-namespace)

long long int new_id = 100000000;

set<SKey> unmatched_segments;
map<SKey, vector<SKey>> matched_segments;
map<SKey, IntNode> node_mapping; //contains node mapping from target terminal nodes to source terminal nodes.
map<SKey, IntNode> node_to_roundabout;

ofstream unm("unmatched_segments.txt");
ofstream check("check.txt");
ofstream notWB("notWB.txt");


struct MatchedCount{
    SKey id;
    uint32_t count;
    MatchedCount(SKey id, uint32_t count){
        this->id = id;
        this->count = count;
    }
};

class MapWrite{
    MatchingUtils mu;
    vector<MatchedCount> unmatched_count;
    map<SKey, SVal> target = segments[1];
    const int buffer_size = 1024000;
    Buffer buffer{buffer_size, Buffer::auto_grow::yes};
    map<SKey, vector<SKey>> target_joins;
    set<SKey> isolated_segments; //not connected to any other segment, for such cases, we don't shift these segments at all.

    void add_tags(Builder* way_builder, const osmium::TagList& tags){
        TagListBuilder tag_builder{buffer, way_builder};
        tag_builder.add_tag( "highway", "footway" );
    }

    void add_tags(Builder* way_builder, const osmium::TagList& new_tags, const osmium::TagList& old_tags){
        TagListBuilder tag_builder{buffer, way_builder};
        for(auto& tag : new_tags){
            if(!(old_tags.has_key(tag.key())))
                tag_builder.add_tag(tag);
            else
                tag_builder.add_tag(tag.key(), old_tags.get_value_by_key(tag.key()));
        }
    }

    void build_way(SVal ref_seg, SVal tar_seg){

        vector<IntNode> source, shifted;

        SVal seg = ref_seg;

        if(!ref_seg)
            seg = tar_seg;

        vector<IntNode> new_int_nodes = seg->get_intermediate_nodes();

        if(ref_seg && tar_seg)
            mapMissingNodes(ref_seg, tar_seg, new_int_nodes);

        if(new_int_nodes.size() <= 1){
            // cout << "Size 1 : " << seg->get_id() << "\n";
            exit(EXIT_FAILURE);
        }

        if(ref_seg || (!ref_seg && tar_seg && isolated_segments.find(tar_seg->get_id()) != isolated_segments.end())){
            source = new_int_nodes;
            shifted = new_int_nodes;

            if( !ref_seg )
                isolated_cnt++;

            for(auto& int_node : new_int_nodes){
                Coordinates c = int_node.get_location();
                osmium::builder::add_node(buffer,
                    _id(int_node.get_id()),
                    _version(1),
                    _timestamp(std::time(nullptr)),
                    _location(osmium::Location{c.x, c.y})
                );
            }
            WayBuilder way_builder{buffer};
            osmium::Way& obj = way_builder.object();

            obj.set_id(new_id++);
            obj.set_version(1);
            way_builder.add_user("fb");
            {
                WayNodeListBuilder wnl_builder{buffer, &way_builder};
                for(auto& int_node : new_int_nodes){
                    Coordinates c = int_node.get_location();
                    wnl_builder.add_node_ref(osmium::NodeRef{int_node.get_id(), osmium::Location{c.x, c.y}});
                }
            }
            if(ref_seg){
                
                add_tags(&way_builder, ref_seg->tags());
                
            }
        }
        else{
            IntNode src1 = new_int_nodes.front(), src2 = new_int_nodes.back();
            IntNode dst1, dst2;
            if(node_to_roundabout.find(src1.get_id()) != node_to_roundabout.end()){
                dst1 = src1; //if roundabout, we don't want it to influence the position of intermediate nodes.
                node_mapping.insert({src1.get_id(), node_to_roundabout.find(src1.get_id())->second});
            }
            else
                dst1 = node_mapping.find(src1.get_id())->second;

            if(node_to_roundabout.find(src2.get_id()) != node_to_roundabout.end()){
                dst2 = src2;
                node_mapping.insert({src2.get_id(), node_to_roundabout.find(src2.get_id())->second});
            }
            else
                dst2 = node_mapping.find(src2.get_id())->second;

            source.push_back( src1 );
            shifted.push_back( dst1 );

            //all this is for the intermediate nodes of the segment.
            Coordinates vec1{dst1.get_location().x - src1.get_location().x, dst1.get_location().y - src1.get_location().y};
            Coordinates vec2{dst2.get_location().x - src2.get_location().x, dst2.get_location().y - src2.get_location().y};

            size_t i=0, size = new_int_nodes.size();
            IntNode prev = src1;
            double d=0, l=seg->get_length();

            for(auto& int_node : new_int_nodes){

                d += distance(prev.get_location(), int_node.get_location());

                if(i > 0 && i < size-1){
                    double r1 = (l-d)/l;
                    double r2 = d/l;
                    double x = (r1*vec1.x + r2*vec2.x) + int_node.get_location().x;
                    double y = (r1*vec1.y + r2*vec2.y) + int_node.get_location().y;
                    Coordinates c{x, y};

                    IntNode inode{new_id++, c};
                    node_mapping.insert({int_node.get_id(), inode});

                }
                // if( int_node.get_id() == 1000283)
                    // cout << "id corresponding to 1000283 ... " << node_mapping.find(int_node.get_id())->second.get_id() << "\n";

                Coordinates c = node_mapping.find(int_node.get_id())->second.get_location();

                source.push_back( int_node );
                shifted.push_back( node_mapping[ int_node.get_id() ] );

                osmium::builder::add_node(buffer,
                    _id(node_mapping.find(int_node.get_id())->second.get_id()),
                    _version(1),
                    _timestamp(std::time(nullptr)),
                    _location(osmium::Location{c.x, c.y})
                );

                prev = int_node;
                ++i;
            }

            WayBuilder way_builder{buffer};
            osmium::Way& obj = way_builder.object();

            // cout << "Reached here in mapwrite ... \n";

            obj.set_id(new_id++);
            obj.set_version(1);
            way_builder.add_user("fb");
            {
                WayNodeListBuilder wnl_builder{buffer, &way_builder};
                for(auto& int_node : new_int_nodes){
                    IntNode inode = node_mapping.find(int_node.get_id())->second;
                    Coordinates c = inode.get_location();
                    wnl_builder.add_node_ref(osmium::NodeRef{inode.get_id(), osmium::Location{c.x, c.y}});
                }
            }
            add_tags(&way_builder, tar_seg->tags());
        }
        calculateMergingAccuracy( source, shifted );
    }

    void calculateMergingAccuracy( vector< IntNode > source, vector< IntNode > shifted ){
        CoordinateArraySequence* cl_source = formCA( source );
        CoordinateArraySequence* cl_shifted = formCA( shifted );

        const Geometry* ls_source = global_factory->createLineString(cl_source);
        Geometry* buff_source = (ls_source->buffer(ce));
        const Geometry* ls_shifted = global_factory->createLineString(cl_shifted);
        // Geometry* buff_shifted = (ls_shifted->buffer(ce));
        if( buff_source->covers( ls_shifted ) )
            covers_cnt++;

        if( !buff_source->intersects( ls_shifted ) )
            not_intersects_cnt++;

        Coordinates c1 = lonlat_to_mercator( shifted.front().get_location() );
        Coordinates c2 = lonlat_to_mercator( shifted.back().get_location() );

        Point* p1 = global_factory->createPoint(Coordinate(pm->makePrecise(c1.x), pm->makePrecise(c1.y)));
        Point* p2 = global_factory->createPoint(Coordinate(pm->makePrecise(c2.x), pm->makePrecise(c2.y)));

        if( buff_source->intersects( p1 ) && buff_source->intersects( p2 ) )
            endpoints_cnt++;

        if( ls_source->buffer(2*ce)->covers( ls_shifted ) )
            covers_cnt2++;

        if( ls_source->buffer(2*ce)->intersects( p1 ) && ls_source->buffer(2*ce)->intersects( p2 ) )
            endpoints_cnt2++;

        if( !ls_source->buffer(2*ce)->intersects( ls_shifted ) )
            not_intersects_cnt2++;

        for(size_t i=1; i<=5; i++){
            buff_source = (ls_source->buffer((double)i*2));
            if( buff_source->covers( ls_shifted ) )
                counts[i-1]++;
        }

        total_cnt++;
    }


    void build_node(Coordinates location){
        NodeBuilder node_builder{buffer};
        node_builder.add_user("fb");
        osmium::Node& obj = node_builder.object();
        // if(new_id == 10001)
        // cout << "build node() ... " << location.x << "\t" << location.y << "\n";
        obj.set_id(new_id++);
        obj.set_location(osmium::Location{location.x, location.y});
    }

    Coordinates getCoordinatesAtFraction(Coordinates c1, double distance, Coordinates c3){

        double bearing = findInitialBearing( c1, c3 );
        return destinationFromStartBearing( c1, bearing, distance );

    }

    Coordinates destinationFromStartBearing( Coordinates from, double bearing, double distance ){
      double angular = distance / RADIUS;
      double bearingRadian = bearing * RadianCoeff;

      double latitudeRadian = from.y * RadianCoeff;
      double longitudeRadian = from.x * RadianCoeff;

      double destLatSine = sin(latitudeRadian) * cos(angular)
          + cos(latitudeRadian) * sin(angular) * cos(bearingRadian);
      double destLatitudeRadian = asin(destLatSine);

      double y = sin(bearingRadian) * sin(angular) * cos(latitudeRadian);
      double x = cos(angular) - sin(latitudeRadian) * destLatSine;
      double destLongitudeRadian = longitudeRadian + atan2(y, x);

      double destLatitude = destLatitudeRadian * DegreesCoeff;
      double destLongitude = destLongitudeRadian * DegreesCoeff;

      return Coordinates{ wrap180(destLongitude), wrap90(destLatitude) };
    }

    double findInitialBearing( Coordinates from, Coordinates to ){
        double fromLatRadian = from.y * RadianCoeff;
        double toLatRadian = to.y * RadianCoeff;
        double longitudeDeltaRadian = (to.x - from.x) * RadianCoeff;

        double x = ( cos(fromLatRadian) * sin(toLatRadian) )
                - ( sin(fromLatRadian) * cos(toLatRadian) * cos(longitudeDeltaRadian) );
        double y = sin(longitudeDeltaRadian) * cos(toLatRadian);
        double t = atan2(y, x);

        double bearing = t * DegreesCoeff;

        return wrap360(bearing);
    }

    double wrap90( double degrees ){
      if (-90 <= degrees && degrees <= 90)
      {
          return degrees;
      }

          // Triangle wave p:360 a:±90 TODO: fix e.g. -315°
      return abs(fmod((fmod(degrees, 360) + 270), 360) - 180) - 90;
    }

    double wrap180( double degrees ){
      if (-180 < degrees && degrees <= 180)
      {
          return degrees;
      }

      // Sawtooth wave p:180, a:±180
      return fmod((degrees + 540), 360) - 180;
    }

    double wrap360( double degrees ){
      if (0 <= degrees && degrees < 360)
      {
          return degrees;
      }

      // Sawtooth wave p:360, a:360
      return fmod((fmod(degrees, 360) + 360), 360);
    }

    Coordinates getCoordinatesAtDistance(Coordinates c1, double distance, Coordinates c3){
        double bearing = findInitialBearing( c1, c3 );
        return destinationFromStartBearing( c1, bearing, distance );
    }

    void mapMissingNodes(SVal ref_seg, SVal tar_seg, vector<IntNode> &new_int_nodes){
        tn qs, qe, ms, me;
        bool reverse = false;

        check_coord << "delta ... " << ref_seg->get_id() << "\n";

        qs = ref_seg->get_first();
        qe = ref_seg->get_second();
        if(distance(qs->get_location(), tar_seg->get_first()->get_location()) < distance(qs->get_location(), tar_seg->get_second()->get_location())){
            ms = tar_seg->get_first();
            me = tar_seg->get_second();
        }
        else{
            ms = tar_seg->get_second();
            me = tar_seg->get_first();
            reverse = true;
        }

        if(!(isWithinCE(qs, ms, ce) && isWithinCE(qe, me, ce)) && !(isWithinCE(qs, me, ce) && isWithinCE(qe, ms, ce)) && !ref_seg->is_roundabout())
            check << ref_seg->get_id() << "\n" << "distance 1 : " << distance(qs->get_location(), ms->get_location()) << "\t" << distance(qe->get_location(), me->get_location())
            << "\t" << "distance 2 : " << distance(qs->get_location(), me->get_location()) << "\t" << distance(qe->get_location(), ms->get_location()) << "\t" << tar_seg->get_id() << " --- \n";

        vector<IntNode> nodeset1 = ref_seg->get_intermediate_nodes(), nodeset2 = tar_seg->get_intermediate_nodes();
        size_t j = 0, cnt = 0;
        size_t size1 = nodeset1.size(), size2 = nodeset2.size();
        if( nodeset1.front().get_id() != nodeset1.back().get_id() ){
          node_mapping.insert({ms->get_id(), ref_seg->get_intermediate_nodes().front()});
          node_mapping.insert({me->get_id(), ref_seg->get_intermediate_nodes().back()});
          j = 1;
        }

        double d = 0, d_prime = 0;
        double l = ref_seg->get_length(), l_prime = tar_seg->get_length();

        new_int_nodes = ref_seg->get_intermediate_nodes();

        if( reverse ){
            std::reverse( nodeset1.begin(), nodeset1.end() );
            std::reverse( new_int_nodes.begin(), new_int_nodes.end() );
        }

        while( j < size2 - 1 ){
            if(terminal_nodes[1].find(nodeset2[j].get_id()) == terminal_nodes[1].end()){
                ++j;
                continue;
            }
            tn source = terminal_nodes[ 1 ][ nodeset2[j].get_id() ];
            size_t pos = mu.getNearestNodesBuffer( source, nodeset1 );
            if( pos < nodeset1.size()-1 ){
                vector< Coordinates > coords;
                coords.push_back( nodeset1[ pos ].get_location() );
                coords.push_back( nodeset1[ pos+1 ].get_location() );
                Coordinates c = getInterpolatedNodeCoordinates(source->get_location(), coords[0], coords[1]);
                mu.insert(source, c, nodeset1, pos, 1);

                IntNode inode{new_id++, c};
                node_mapping.insert({nodeset2[j].get_id(), inode});
            }
            ++j;
        }

    }

public:

    size_t covers_cnt = 0;
    size_t endpoints_cnt = 0;
    size_t not_intersects_cnt = 0;
    size_t total_cnt = 0;
    size_t isolated_cnt = 0;

    int counts[5] = {0, 0, 0, 0, 0};

    size_t covers_cnt2 = 0;
    size_t endpoints_cnt2 = 0;
    size_t not_intersects_cnt2 = 0;

    static bool sortMatchedCount(MatchedCount seg1, MatchedCount seg2){
        return seg1.count > seg2.count;
    }

    void merge(vector<JoinResult> join_result, char* outfile_name){
        osmium::io::File output_file{outfile_name};
        Writer writer{output_file, osmium::io::overwrite::allow};

        findUnmatchedSegments(join_result);

        drawMatched(join_result);

        drawEdgeCases();

        set<SKey> visited;
        for(auto& seg : unmatched_segments){
            if(visited.find(seg) != visited.end())
                continue;

            vector<map<SKey, IntNode>> nodes = findMissingNodes(seg, visited); //these are the terminal nodes that are connected to only all missing segments.

            map<SKey, IntNode> mapped_terminal_nodes = nodes[0];
            map<SKey, IntNode> missing_terminal_nodes = nodes[1];

            if(missing_terminal_nodes.size() == 0)
                continue;
            if(mapped_terminal_nodes.size() == 0){
                for(auto& node : missing_terminal_nodes)
                    node_mapping.insert({node.first, IntNode{new_id++, node.second.get_location()}});
                continue;
            }
            for(auto& node : missing_terminal_nodes){
                double d = 0;
                Coordinates c = node.second.get_location();
                double x = 0, y = 0;
                for(auto& it : mapped_terminal_nodes){
                    Coordinates src = it.second.get_location();
                    Coordinates dst = node_mapping.find(it.first)->second.get_location();
                    double dist = distance(c, src);
                    d += (1/dist);
                    Coordinates vec{dst.x-src.x, dst.y-src.y};
                    x += ((1/dist)*vec.x);
                    y += ((1/dist)*vec.y);
                }
                Coordinates c_final{c.x + (x/d), c.y + (y/d)};
                node_mapping.insert({node.first, IntNode(new_id++, c_final)});
            }
        }

        for(auto& id : unmatched_segments){
            SVal seg = target.find(id)->second;

            build_way(NULL, seg);
            assert(buffer.is_aligned());
        }

        writer(std::move(buffer));
        writer.close();
    }

    void drawEdgeCases(){
        for(auto& seg : edge_cases){
            build_way(segments[0].find(seg)->second, NULL);
            buffer.commit();
        }
    }

    vector<map<SKey, IntNode>> findMissingNodes(SKey seg, set<SKey> &visited){
        map<SKey, IntNode> missing_terminal_nodes;
        map<SKey, IntNode> mapped_terminal_nodes;
        vector<map<SKey, IntNode>> nodes;
        list<SKey> queue;

        if(target_joins.find(seg)->second.size() == 0){
            // cout << "Isolated: " << seg << "\n";
            isolated_segments.insert(seg);
            nodes.push_back(mapped_terminal_nodes);
            nodes.push_back(missing_terminal_nodes);
            return nodes;
        }

        queue.push_back(seg);

        while(!queue.empty()){
            SKey s = queue.front();
            visited.insert(s);
            queue.pop_front();

            SVal segment = target.find(s)->second;
            assert(segment);

            size_t connected_to_roundabout = 0;
            SVal roundabout;

            vector<SKey> joins = target_joins.find(s)->second;

            for(auto& it : joins){
                if(target.find(it)->second->is_roundabout() && unmatched_segments.find(it) == unmatched_segments.end()){ //roundabout is matched but a segment connected to this roundabout is missing.
                    connected_to_roundabout++;
                    SKey roundabout_id = matched_segments.find(it)->second[0];
                    roundabout = segments[0].find(roundabout_id)->second;
                }
                if(unmatched_segments.find(it) == unmatched_segments.end() || visited.find(it) != visited.end())
                    continue;
                queue.push_back(it);
            }

            IntNode first = segment->get_intermediate_nodes().front();
            IntNode second = segment->get_intermediate_nodes().back();

            //TODO: cases when both ends of the missing segment are connected to roundabouts.
            if(connected_to_roundabout == 1){
                assert(roundabout);
                // cout << "Size : " << roundabout->get_id() << "\n";
                double min_d1 = DBL_MAX, min_d2 = DBL_MAX;
                Coordinates c1 = first.get_location(), c2 = second.get_location();
                IntNode n1, n2;
                for(auto& inode : roundabout->get_intermediate_nodes()){
                    double d1 = distance(c1, inode.get_location());
                    double d2 = distance(c2, inode.get_location());
                    if(d1 < min_d1){
                        min_d1 = d1;
                        n1 = inode;
                    }
                    if(d2 < min_d2){
                        min_d2 = d2;
                        n2 = inode;
                    }
                }
                if(min_d1 < min_d2){
                    node_to_roundabout.insert({first.get_id(), n1});
                    if(node_mapping.find(second.get_id()) == node_mapping.end())
                        missing_terminal_nodes.insert({second.get_id(), second});
                    else
                        mapped_terminal_nodes.insert({second.get_id(), second});
                }
                else if(min_d2 < min_d1){
                    node_to_roundabout.insert({second.get_id(), n2});
                    if(node_mapping.find(first.get_id()) == node_mapping.end())
                        missing_terminal_nodes.insert({first.get_id(), first});
                    else
                        mapped_terminal_nodes.insert({first.get_id(), first});
                }
                continue;
            }

            if(node_mapping.find(first.get_id()) == node_mapping.end())
                missing_terminal_nodes.insert({first.get_id(), first});
            else
                mapped_terminal_nodes.insert({first.get_id(), first});

            if(node_mapping.find(second.get_id()) == node_mapping.end())
                missing_terminal_nodes.insert({second.get_id(), second});
            else
                mapped_terminal_nodes.insert({second.get_id(), second});

        }
        nodes.push_back(mapped_terminal_nodes);
        nodes.push_back(missing_terminal_nodes);
        return nodes;
    }

    // in draw matched, add the corresponding terminal nodes from target to the reference database - never mind, we are already doing this in build_way()
    // however, the position of the terminal nodes is apparently not decided by the actual position of corresponding terminal nodes?

    void drawMatched(vector<JoinResult> join_result){

        check << "All segments that have matches not within the bounds : \n\n";

        map<SKey, TVal> output = join_result[0].result;
        for(auto& it : output){
            SKey id = it.first;
            TVal tuple = it.second;
            SVal seg = tuple->segment();

            if(seg->get_id() != 0) //to add the tags accordingly.
                build_way(segments[0].find(id)->second, seg);
            else
                build_way(segments[0].find(id)->second, NULL);
            assert(buffer.is_aligned());
            buffer.commit();
        }
    }

    bool endsWithinBounds(SVal delta, SVal cs){
        if(delta->is_roundabout())
            return true;
        if(cs->get_id() == 0)
            return true;
        tn s1 = delta->get_first(), e1 = delta->get_second(), s2 = cs->get_first(), e2 = cs->get_second();
        if(isWithinCE(s1, s2, ce) && isWithinCE(e1, e2, ce))
            return true;
        if(isWithinCE(s1, e2, ce) && isWithinCE(e1, s2, ce))
            return true;
        notWB << delta->get_id() << "\n";
        return false;
    }

    
    void findUnmatchedSegments(vector<JoinResult> join_result){
        if(join_result.size() == 0)
            return;

        map<SKey, TVal> output = join_result[0].result;
        // map<SKey, vector<SVal>> to_be_added; //key is the segment to be removed from the target database and the corresponding vector are the segments to be added in the target database.
        map<SKey, vector<SVal>> to_be_added; // key is the parent cs id that is interpolated, value is the vector of interpolated node ids on the corresponding cs.
        vector<SKey> to_be_erased;
        // map< SKey, map< SKey, vector<SVal>>> end_points_matches;
        for(auto& it : output){
            SKey delta = it.first;
            SKey matched = it.second->segment()->get_id();

            // if( delta == 1005 ){
                // cout << "int nodes for match of 1005 ... \n";
                // vector< IntNode > all_int_nodes = it.second->segment()->get_intermediate_nodes();
                // for( auto it : all_int_nodes )
                    // cout << it.get_location().x << "\t" << it.get_location().y << "\n";
            // }

            endsWithinBounds(segments[0].find(delta)->second, it.second->segment());

            if(matched_segments.find(matched) != matched_segments.end())
                matched_segments.find(matched)->second.push_back(delta);
            else
                matched_segments.insert({matched, {delta}});

            if(matched >= 0)
                continue;

            if(matched < 0){    //the segment can be splitted or sef-joined. In either case, get the non-negative (original) parents.
                SKey first = it.second->segment()->get_first()->get_id();
                SKey second = it.second->segment()->get_second()->get_id();
                if(first > 0 && second > 0){         //self-joined segment.
                    target.insert({matched, it.second->segment()});
                    vector<SKey> parents = it.second->segment()->parents();
                    for(auto& item : parents){
                        if(matched_segments.find(item) != matched_segments.end())
                            matched_segments.find(item)->second.push_back(delta);
                        else
                            matched_segments.insert({item, {delta}});
                    }
                }
                else{
                    vector<SKey> parents = it.second->segment()->parents();
                    size_t size = parents.size(), i = 0;
                    if(delta == 1960)
                        checking << "Parents for " << delta << " : " << size << "\n";
                    while(i < size){
                        SKey parent_id = parents[i];
                        if(delta == 1960)
                            checking << parent_id << "\n";
                        if(parent_id > 0){
                            // if(parent_id == 3498)
                                // cout << "Ohh no! 3498 inserted in matched, in parent!\n\n";
                            if(matched_segments.find(parent_id) != matched_segments.end())
                                matched_segments.find(parent_id)->second.push_back(delta);
                            else
                                matched_segments.insert({parent_id, {delta}});
                        }
                        else
                            to_be_erased.push_back(parent_id);
                        i++;
                    }
                    checking << "\n";


                    if( first < 0 ){
                        SKey parent_cs = parent_interpolated_segments[ first ];
                        if( to_be_added.find( parent_cs ) == to_be_added.end() )
                            to_be_added[ parent_cs ] = {};

                        vector< SVal > segs = interpolated_segments[ first ];
                        to_be_added[ parent_cs ].insert( to_be_added[ parent_cs ].end(), segs.begin(), segs.end() );
                        to_be_added[ parent_cs ].push_back( it.second->segment() );
                    }

                    if( second < 0 ){
                        SKey parent_cs = parent_interpolated_segments[ second ];
                        if( to_be_added.find( parent_cs ) == to_be_added.end() )
                            to_be_added[ parent_cs ] = {};

                        vector< SVal > segs = interpolated_segments[ second ];
                        to_be_added[ parent_cs ].insert( to_be_added[ parent_cs ].end(), segs.begin(), segs.end() );
                        to_be_added[ parent_cs ].push_back( it.second->segment() );
                    }

                }
            }
        }

        checking << "number of interpolated nodes on each segment ... \n";
        for( auto it : to_be_added ){
            SKey cs_id = it.first;
            vector< SVal > segs = it.second;
            SVal cs = segments[ 1 ][ cs_id ];
            if( segs.size() == 2 ){
                target.insert( { segs[ 0 ]->get_id(), segs[ 0 ] } );
                target.insert( { segs[ 1 ]->get_id(), segs[ 1 ] } );
                break;
            }
            tn first = cs->get_first();
            while( first != cs->get_second() ){
                double min_length = DBL_MAX;
                SVal seg;
                for( auto it : segs ){
                    if( it->get_first() == first && it->get_length() < min_length ){
                        min_length = it->get_length();
                        seg = it;
                    }
                }
                if( seg ){
                    target.insert( { seg->get_id(), seg } );
                    first = seg->get_second();
                }
                else
                    break;
            }
        }

        for(auto& it : to_be_erased){
            target.erase(it);
        }
        for(auto& it : target){
            SKey id = it.first;
            SVal seg = it.second;

            SKey first = seg->get_first()->get_id(), second = seg->get_second()->get_id();
            if(relaxed_joins.find(first) != relaxed_joins.end() && relaxed_joins.find(first)->second == second && matched_segments.find(id) != matched_segments.end())
                continue;
            if(matched_segments.find(id) == matched_segments.end()){
                unmatched_segments.insert(id);
                unm << id << "\t" << seg->get_first()->get_id() << "\t" << seg->get_second()->get_id() << "\n";
            }
        }
        connectedToUnmatched(target);
    }

    void connectedToUnmatched(map<SKey, SVal> &target){
        set<SKey> empty;
        target_joins = findJoins(target, empty);
        for(auto& seg : unmatched_segments){
            uint32_t ctr = 0;
            vector<SKey> all_joins = target_joins.find(seg)->second;
            for(auto& it : all_joins){
                if(matched_segments.find(it) != matched_segments.end())
                    ctr++;
            }
            unmatched_count.push_back(MatchedCount(seg, ctr));
        }
        // cout << "Size of joins : " << target_joins.size() << "\n\n";
        for(auto it1 : target_joins){
            adjlist2 << it1.first << " \n";
            for(auto& it2 : it1.second)
                adjlist2 << "\t" << it2 << "\n";
        }
        sort(unmatched_count.begin(), unmatched_count.end(), sortMatchedCount);
    }
};


#endif
