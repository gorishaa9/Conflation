#include <sstream>
#include <thread>

#include <osmium/io/any_input.hpp>
#include <osmium/visitor.hpp>
#include <osmium/handler.hpp>
#include <osmium/osm/node_ref.hpp>
#include <osmium/util/iterator.hpp>
#include <osmium/object_pointer_collection.hpp>
#include <osmium/osm/object.hpp>
#include <osmium/osm/segment.hpp>
#include <osmium/osm/location.hpp>

#include <osmium/index/map/sparse_mem_array.hpp>
#include <osmium/handler/node_locations_for_ways.hpp>

#include "Initialize.hpp"
#include "terminal_node.hpp"
#include "TerminalNodesHandler.hpp"
#include "CountHandler.hpp"
#include "segment.hpp"
#include "SegmentHandler.hpp"
#include "../index/Insert.hh"
// #include "FindMatches_angle.hpp"
#include "FindMatches_wprints.hpp"

// #include "RJReGroup_wprints.hpp"
#include "RJReGroup_ablation.hpp"
// #include "RJReGroup_noDG.hpp"

#include "MapWrite.hpp"

#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <fstream>
#include <boost/shared_ptr.hpp>

using namespace std;
using namespace std::chrono;

using key_type = osmium::object_id_type;
using value_type = int;

using SKey = long long int;
using SVal = boost::shared_ptr<Segment_creator>;
using RVal = boost::shared_ptr<Relation>;

vector<std::map<SKey, SVal>> segments; //segments[0] is the map of segments from r4efrence database. segment[1] is the map of segments from target database.
int dataset = 0; //reference database.
vector<std::map<key_type, value_type>> nodeset;
vector<std::map<key_type, boost::shared_ptr<Terminal_node>>> terminal_nodes;
map<SKey, RVal> relations;

ofstream result;
// ofstream mergingAccuracy;

int main(int argc, char* argv[]){ //First argument is the reference database, second is the target database.
  CountHandler handler;
  SegmentHandler segHandler;
  TerminalNodesHandler tn_handler;
  RoadLengthHandler road_length_handler;

  Initialize ini;
  ini.initialize();

  string str1(argv[1]);
  string str2(argv[2]);

  total_length.push_back( 0 );
  total_length.push_back( 0 );

  size_t l1 = str1.length(), l2 = str2.length();
  string sub1 = str1.substr(l1-7,3), sub2 = str2.substr(l2-7,3);

  char *source1, *source2;
  source1 = &sub1[ 0 ];
  source2 = &sub2[ 0 ];

  string str( "result" );
  str += argv[5];
  str += ".txt";


  result.open( str );
  // mergingAccuracy.open( str_ma );

  cout << "source1 ... " << source1 << "\tsource2 ... " << source2 << "\n";

  result << "source 1 ... " << str1 << "\n";
  result << "source 2 ... " << str2 << "\n\n";

  namespace map = osmium::index::map;
  using index_type = map::SparseMemArray<osmium::unsigned_object_id_type, osmium::Location>;
  using location_handler_type = osmium::handler::NodeLocationsForWays<index_type>;

  index_type index;
  location_handler_type location_handler{index};

  data_source = source1;

  auto start1 = high_resolution_clock::now();

  osmium::io::Reader readerp{argv[1], osmium::osm_entity_bits::way | osmium::osm_entity_bits::node};
  osmium::apply( readerp, location_handler, road_length_handler );
  readerp.close();

  osmium::io::Reader reader{argv[1], osmium::osm_entity_bits::way | osmium::osm_entity_bits::node};
  osmium::apply(reader, location_handler, handler);
  reader.close();

  osmium::io::Reader reader2{argv[1], osmium::osm_entity_bits::way | osmium::osm_entity_bits::node};
  osmium::apply(reader2, location_handler, tn_handler);
  reader2.close();

  osmium::io::Reader readerW{argv[1], osmium::osm_entity_bits::way | osmium::osm_entity_bits::node};
  osmium::apply(readerW, location_handler, segHandler);
  readerW.close();

  if( strcmp( source1, "OSM") == 0 )
      segHandler.empty_roundabouts();

  auto stop1 = high_resolution_clock::now();
  auto duration1 = duration_cast<milliseconds>(stop1 - start1);

  result << "Time taken to read source1 ... " << duration1.count() << " ms.\n\n";

  cout << "Finished reading source1 ... \n";
  cout << "Total length of roads in source1 ... " << total_length[ 0 ] << "\n\n";
  //
  dataset = 1;  //target database.
  data_source = source2;
  cout << "length of small ways ... " << small_ways.size() << "\n";

  auto start2 = high_resolution_clock::now();

  osmium::io::Reader readerp_d{argv[2], osmium::osm_entity_bits::way | osmium::osm_entity_bits::node};
  osmium::apply( readerp_d, location_handler, road_length_handler );
  readerp_d.close();

  osmium::io::Reader reader_d{argv[2], osmium::osm_entity_bits::way | osmium::osm_entity_bits::node};
  osmium::apply(reader_d, location_handler, handler);
  reader_d.close();

  osmium::io::Reader reader2_d{argv[2], osmium::osm_entity_bits::way | osmium::osm_entity_bits::node};
  osmium::apply(reader2_d, location_handler, tn_handler);
  reader2_d.close();

  osmium::io::Reader readerW_d{argv[2], osmium::osm_entity_bits::way | osmium::osm_entity_bits::node};
  osmium::apply(readerW_d, location_handler, segHandler);
  readerW_d.close();


  if( strcmp( source2, "OSM") == 0 )
      segHandler.empty_roundabouts();

  auto stop2 = high_resolution_clock::now();
  auto duration2 = duration_cast<milliseconds>(stop2 - start2);
  result << "Time taken to read source2 ... " << duration2.count() << " ms.\n\n";

  cout << "length of small ways ... " << small_ways.size() << "\n";

  cout << "Total length of roads in source2 ... " << total_length[ 1 ] << "\n\n";

  std::cout << "Done Reading!" << "\n";

  std::cout << "Size of reference dataset: " << segments[0].size() << "\n";
  std::cout << "Size of target dataset: " << segments[1].size() << "\n\n";

  result << "Size of reference dataset: " << segments[0].size() << "\n";
  result << "Size of target dataset: " << segments[1].size() << "\n\n";

  result << "Total length of roads in reference dataset ... " << total_length[ 0 ] << "\n";
  result << "Total length of roads in target dataset ... " << total_length[ 1 ] << "\n\n";

  result << "Number of terminal nodes in reference dataset ... " << terminal_nodes[ 0 ].size() << "\n";
  result << "Number of terminal nodes in the target dataset ... " << terminal_nodes[ 1 ].size() << "\n\n";

  // return 0;

  // for( auto it : segments[ 0 ])
  //     cout << it.second->get_length() << "\t" << log10( it.second->get_length() ) << "\n";
  // return 0;

  std::cout << "\n" << "Bulk inserting in R tree of segments..." << "\n"; //we have to insert the segments from target database to the RTree.

  auto start3 = high_resolution_clock::now();

  rtree_insert(segments[1]);

  input_file = argv[5];

  auto stop3 = high_resolution_clock::now();
  auto duration3 = duration_cast<milliseconds>(stop3 - start3);
  result << "Time taken to insert segments in r-tree ... " << duration3.count() << " ms.\n\n";

  std::cout << "\n" << "Inserted the segments in rtree." << "\n";
  MatchingUtils mu;
  std::cout << "\n" << "Populating the relation sets of the segments ..." << "\n";

  auto start4 = high_resolution_clock::now();

  mu.populate(segments[0], false);   //the relation sets will be from the query database, the reference database.


  mu.splitQuery(segments[0], segments[ 0 ], relations);
  for( auto it : splitted_segments )
      cout << " splitted segment ... " << it.first << "\n";

  result << "splitted_segments size ... " << splitted_segments.size() << "\n";
  // input_file = static_cast<size_t>(stoul(argv[5]));


  mu.populate(splitted_segments, false);
  size_t ctr = 0;
  // mu.recomputeScores(splitted_querys, relations);
  while( splitted_segments.size() > 0 && ctr <= 2 ){
      splitted_querys.clear();
      splitted_querys.insert( splitted_segments.begin(), splitted_segments.end() );
      splitted_segments.clear();
      mu.splitQuery( segments[ 0 ], splitted_querys, relations );
      mu.populate( splitted_segments, false );
      cout << "splitted segments size ... " << splitted_segments.size() << "\n";
      result << "splitted_segments size ... " << splitted_segments.size() << "\n";
      // return 0;
      ctr++;
  }
  mu.pruneCandidates( relations );

  auto stop4 = high_resolution_clock::now();
  auto duration4 = duration_cast<milliseconds>(stop4 - start4);
  result << "Time taken to find matches ... " << duration4.count() << " ms.\n\n";
  // mu.removeDuplicateNodes( relations );
  // mu.removeDummySegments( relations, 15 );
  mu.finalWrite(segments[0], relations);
  //
  //
  for( auto it : splitted_segments )
      cout << " splitted segment ... " << it.first << "\n";

  cout << "Size of splitted querys ... " << splitted_querys.size() << "\n";
  cout << "Size of splitted segments ... " << splitted_segments.size() << "\n";

  std::cout << "Size of relations ... " << relations.size() << "\n";
  std::cout << "\n" << "Populating done, yay!" << "\n";

  size_t cnt = 0;
  for( auto it : relations )
      if( it.second->get_entries().size() == 1 )
          cnt++;

  cout << "Number of relations with only dummy segment ... " << cnt << "\n";
  result << "Number of relations with only dummy segment ... " << cnt << "\n";

  size_t total_matches_count = 0;
  for( auto it : relations )
      total_matches_count += ( it.second->get_entries().size() - 1 );

  result << "Total number of matches in " << ( relations.size() - cnt ) << " relations : " << total_matches_count << "\n\n";

  // return 0;

  std::cout << "\n" << "Starting Rank Join --- " << "\n";

  auto start5 = high_resolution_clock::now();


  RankJoinCA rj;
  vector<JoinResult> join_result = {};
  std::stringstream ss(argv[6]);
  bool b;
  if(!(ss >> std::boolalpha >> b)) {
    cout << "argv[6] ... " << argv[6] << "\n";
    cout << "Parsing error ... \n";
    return 0;
  }
  auto term_time = static_cast<int>(stoul(argv[7]));

  std::thread rjt( &RankJoinCA::processRankJoin, &rj, std::ref( relations ), std::ref( argv[4] ), std::ref( argv[5] ), std::ref( b ), std::ref( join_result ) );
  // rj.processRankJoin(relations, argv[4], argv[5], b, join_result);
  while( join_result.size() == 0 ){
      if( duration_cast<milliseconds>(high_resolution_clock::now() - start5).count() > term_time ){
          result << "Time limit exceeded!!\n";
          quit_flag = true;
          // std::this_thread::sleep_for(std::chrono::seconds(10));
          rjt.join();
          return 0;
      }
  }
  rjt.join();

  auto stop5 = high_resolution_clock::now();
  auto duration5 = duration_cast<milliseconds>(stop5 - start5);

  result << "Rank Join completed in ... " << duration5.count() << " ms.\n\n";

  cout << "\n\nRank Join completed in " << duration5.count() << " ms.\n";

  result << "Number of groups ... " << join_order_group.size() << "\n\n";

  result << "Group size distrbution ... \n\n";

  for( auto it : relation_groups ){
      result << it.first << "\t" << it.second.size() << "\n";
  }
  result << "\n";

  cout << "Map matching accuracy ... \n\n";

  rj.calculateMatchingAccuracy( join_result );

  result << "Totally covers count within 5 metres ... " << rj.covers_cnt << "\n\n";
  result << "End points within buffer count within 5 metres ... " << rj.endpoints_cnt << "\n\n";
  result << "Not intersects at all within 5 metres ... " << rj.not_intersects_cnt << "\n\n";
  result << "Totally covers count within 10 metres ... " << rj.covers_cnt2 << "\n\n";
  result << "End points within buffer count within 10 metres ... " << rj.endpoints_cnt2 << "\n\n";
  result << "Not intersects at all within 10 metres ... " << rj.not_intersects_cnt2 << "\n\n";
  result << "Total segments count ... " << rj.total_cnt << "\n\n";
  result << "Matched to dummy segments count ... " << rj.match_dummy_cnt << "\n\n";
  result << "Single matches count ... " << single_matches.size() << "\n\n";
  result << "All segments matched to dummy segments and not in single matches ..." << ( not_matched.size() - single_matches.size() ) << "\n\n";
  for( auto it : not_matched ){
      if( single_matches.find( it ) == single_matches.end() )
        result << it << "\t" << segments[ 0 ][ it ]->get_first()->get_id() << "\t" << segments[ 0 ][ it ]->get_second()->get_id() << "\n";
  }

  result << "segments to verify ... " << to_verify_matches.size() << "\n\n";
  for( auto it : to_verify_matches ){
    result << it << "\t" << segments[ 0 ][ it ]->get_first()->get_id() << "\t" << segments[ 0 ][ it ]->get_second()->get_id() << " : " <<
    join_result[ 0 ].result[ it ]->segment()->get_id() << "\t" << join_result[ 0 ].result[ it ]->segment()->get_first()->get_id() << "\t" <<
    join_result[ 0 ].result[ it ]->segment()->get_second()->get_id() << "\n";
  }

  cout << "\n";

  size_t total_tuples = 0, total_seen_tuples = 0;

  for( auto it : relations ){
      total_tuples += it.second->get_entries().size();
      total_seen_tuples += seen_tuples[ it.first ];
  }

  result << "All tuples present for rank join ... " << total_tuples << "\n";

  result << "All seen tuples ... " << total_seen_tuples << "\n";

  result << "\n\n";
  // return 0;
  cout << "\n\nStarting map merging ... \n";
  MapWrite mw;
  auto start6 = high_resolution_clock::now();
  mw.merge(join_result, argv[3]);

  auto stop6 = high_resolution_clock::now();
  auto duration6 = duration_cast<milliseconds>(stop6 - start6);

  result << "Time taken by map merging ... " << duration6.count() << " ms.\n\n";

  result << "Merging accuracy ... \n\n";

  result << "Totally covers count within 5 metres ... " << mw.covers_cnt << "\n\n";
  result << "End points within buffer count within 5 metres ... " << mw.endpoints_cnt << "\n\n";
  result << "Not intersects at all within 5 metres ... " << mw.not_intersects_cnt << "\n\n";
  result << "Totally covers count within 10 metres ... " << mw.covers_cnt2 << "\n\n";
  result << "End points within buffer count within 10 metres ... " << mw.endpoints_cnt2 << "\n\n";
  result << "Not intersects at all within 10 metres ... " << mw.not_intersects_cnt2 << "\n\n";
  result << "Total segments count ... " << mw.total_cnt << "\n\n";
  result << "Isolated segments count ... " << mw.isolated_cnt << "\n\n";
  result << "Segments within 2 meter ... " << mw.counts[ 0 ] << "\n";
  result << "Segments within 4 meter ... " << mw.counts[ 1 ] << "\n";
  result << "Segments within 6 meter ... " << mw.counts[ 2 ] << "\n";
  result << "Segments within 8 meter ... " << mw.counts[ 3 ] << "\n";
  result << "Segments within 10 meter ... " << mw.counts[ 4 ] << "\n";
}
