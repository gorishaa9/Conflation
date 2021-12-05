#ifndef FINDMATCHES_H
#define FINDMATCHES_H

#include "terminal_node.hpp"
//#include "CountHandler.hpp"
#include "segment.hpp"
#include "SegmentHandler.hpp"
#include "../index/Insert.hh"
#include "Util.hpp"
#include "Relation.hpp"
#include "IntNode.hpp"
#include "Tuple.hpp"
#include <geos/geom/Geometry.h>

#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>

#define PI 3.14159265

using namespace std;
using namespace geos::geom;

using TVal = boost::shared_ptr<Tuple>;
using RVal = boost::shared_ptr<Relation>;
using tn = boost::shared_ptr<Terminal_node>;
extern map<SKey, RVal> relations;
char* input_file;

map<SKey, SKey> node_segs;
map<SKey, SVal> splitted_segments; //these are splitted query segments.
map<SKey, SVal> splitted_querys;
// map<SKey, map<SKey, vector<tn>>> modified_nodes; //first SKey is source_nodeId, second SKey is the candidate_segId.
long long int modified_nodeId = -1000;
long long int modified_segId = -1000;
// map<SKey, vector<SKey>> interpolation; //key is the id of original segment in target database, value is the id of interpolated segments.
map<SKey, vector<SVal>> interpolated_segments; //key is the interpolated node id, value are the 2 splitted segments from the interpolated node.
map< SKey, SKey > parent_interpolated_segments; // key is the interpolated node id of the cs and value is the id of the cs that is interpolated.
//In both the interpolated and self-joined segments, the key and value both will be segments from target database.
map< SKey, map< SKey, tn > > modified_nodes; //first SKey - source_nodeID, second SKey - parent csID, tn - modified terminal node.
map< SKey, map< SKey, SVal > > end_points_segments; // two SKeys are two end points of the corresponding segment SVal.

// RTree is formed from segments in "target" database.
ofstream testFile("testFile20.txt");
ofstream itest("itest.txt");
ofstream mnodes("mnodes.txt");
ofstream checking("checking.txt");
// ofstream test_cs("test_cs.txt");
// ofstream test_delta("test_delta.txt");

boost::shared_ptr<Terminal_node> identity_node(new Terminal_node(-1, Coordinates(-1, -1), -1));

class MatchingUtils{

public:
    static bool sortRelationTuples(boost::shared_ptr<Tuple> tuple1, boost::shared_ptr<Tuple> tuple2){
        return tuple1->score() > tuple2->score();
    }

    static bool score_sort( double score1, double score2 ){
        return score1 > score2;
    }

    void populate(map<SKey, SVal> &querys, bool full_pass){     //input is the query database -- OSM.

        long int total_count=0;
        for(auto it : querys){

            cout << "\n\n\nStarting with ------------------------- " << it.first << "\t" << it.second->get_first()->get_id() << "\t" << it.second->get_second()->get_id() << "\n\n\n";

            boost::shared_ptr<Segment_creator> delta = it.second;
            std::list<long long int> candidateMatches = rtree_search(delta);
            boost::shared_ptr<Relation> rel(new Relation(delta->get_id()));
            boost::shared_ptr<Rectangle> box;
            Geometry* buff;
            vector<IntNode> cl;
			      const osmium::TagList tags;
            boost::shared_ptr<Segment_creator> id_seg(new Segment_creator(0, identity_node, identity_node, delta->get_length(), 0, box, buff, cl, 0, tags));
            boost::shared_ptr<Tuple> dummy(new Tuple(id_seg, 0, true));
            rel->insert(dummy);

            map< SKey, map< SKey, vector< SVal > > > map_entries;

      			// if(delta->get_id() == 1005){
        		// 		cout << "All matches : \n";
        		// 		for(auto& i : candidateMatches){
          	// 				SVal cs = segments[1].find(i)->second;
          	// 				cout << cs->get_id() << "\t" << cs->get_first()->get_id() << "\t" << cs->get_second()->get_id() << "\n";
        		// 		}
      			// }
            // cout << "candidates size before ... " << candidateMatches.size() << "\n";
  			    filterCandidates(delta, candidateMatches);
            // cout << "candidates size after ... " << candidateMatches.size() << "\n";
      			// if(delta->get_id() == 1005){
        		// 		cout << "Filtered matches : \n";
        		// 		for(auto& i : candidateMatches){
          	// 				SVal cs = segments[1].find(i)->second;
          	// 				cout << cs->get_id() << "\t" << cs->get_first()->get_id() << "\t" << cs->get_second()->get_id() << "\n";
        		// 		}
      			// }

            if( delta->get_first()->get_id() == delta->get_second()->get_id() ){
                cout << "searching polygon match for ... " << delta->get_id() << "\n";
                SVal seg = polygonMatch( delta, candidateMatches );
                if( seg ){
                  createMatches( delta, seg, rel, candidateMatches, map_entries );
                  // draw_segment( seg );
                  cout << "got a joined roundabout ... " << seg->get_id() << "\t" << seg->get_first()->get_id() << "\t" << seg->get_second()->get_id() << "\n";
                }
                relations.insert({it.first, rel});
                continue;
            }

            // if(delta->get_id() == 1335){
            //     cout << "is within angle 1335 and 2735 ... " << withinAngle( delta, segments[1][2735]) << "\n";
            // }

            for(auto& i : candidateMatches){
                SVal cs = segments[1].find(i)->second;

          				if(delta->is_roundabout()){
          					double score = similarityScore(delta, cs);
          					if(score == 0.0)
          						continue;
          					boost::shared_ptr<Tuple> tuple(new Tuple(cs, score, false));
          					rel->insert(tuple);
          					continue;
          				}

          				if(cs->get_length() >= ce && !withinAngle(delta, cs)){
                      if( delta->get_id() == 1001 && i == 2526 )
                          cout << "not within angle 2526 ... \n";
                      continue;
                  }

                  if( delta->get_id() == 1001 && i == 2526 )
                      cout << "within angle 2526 ... \n";

                  tn qs = delta->get_first(), qe = delta->get_second(), ms = cs->get_first(), me = cs->get_second();

                  if((isWithinCE(qs, ms, ce) && isWithinCE(qe, me, ce))
          				|| (isWithinCE(qs, me, ce) && isWithinCE(qe, ms, ce))){

                      cout << "Wrong place in populate() ... 1" << "\n";
                      createMatches(delta, cs, rel, candidateMatches, map_entries);
                  }

                  else if(!isWithinCE(qs, ms, ce) && !isWithinCE(qe, me, ce) && !isWithinCE(qs, me, ce) && !isWithinCE(qe, ms, ce)){

                      cout << "Correct place in populate() ... 3" << "\n";

                      vector<SVal> cs_primes = determineCase(delta, cs, candidateMatches, map_entries);
                      if( find( cs_primes.begin(), cs_primes.end(), cs ) == cs_primes.end())
                          cs_primes.push_back( cs );
                      for(auto cs_prime : cs_primes){
                          cout << "in ... " << cs_prime->get_id() << "\n";
                          if(cs_prime && isPotential(delta, cs_prime) && withinAngle(delta, cs_prime)){
                              cout << "out ... " << cs_prime->get_id() << "\n";
                              createMatches(delta, cs_prime, rel, candidateMatches, map_entries);
                          }
                      }
                  }
                  else if(isWithinCE(qs, ms, ce) && isWithinCE(qs, me, ce) && !isWithinCE(qe, ms, ce) && !isWithinCE(qe, me, ce)
          					||	isWithinCE(qe, ms, ce) && isWithinCE(qe, me, ce) && !isWithinCE(qs, ms, ce) && !isWithinCE(qs, me, ce)){

          						cout << "not As expected for 2039. \n";
            					relaxed_joins.insert({cs->get_first()->get_id(), cs->get_second()->get_id()});
            					relaxed_joins.insert({cs->get_second()->get_id(), cs->get_first()->get_id()});
                      if( !withinAngle( delta, cs ) )
                          continue;
                      vector< SVal > cs_primes = determineCase(delta, cs, candidateMatches, map_entries);
                      if( find( cs_primes.begin(), cs_primes.end(), cs ) == cs_primes.end())
                          cs_primes.push_back( cs );
                      for(auto cs_prime : cs_primes){
                          cout << "in ... " << cs_prime->get_id() << "\n";
                          if(cs_prime && isPotential(delta, cs_prime) && withinAngle(delta, cs_prime)){
                              cout << "out ... " << cs_prime->get_id() << "\n";
            		              createMatches(delta, cs_prime, rel, candidateMatches, map_entries);
                          }
                      }
                  }

          				else if((isWithinCE(qs, ms, ce) && !isWithinCE(qs, me, ce) && !isWithinCE(qe, me, ce) && !isWithinCE(qe, ms, ce))
          					 || (!isWithinCE(qs, ms, ce) && !isWithinCE(qs, me, ce) && isWithinCE(qe, me, ce) && !isWithinCE(qe, ms, ce))
          					 || (!isWithinCE(qs, ms, ce) && !isWithinCE(qs, me, ce) && !isWithinCE(qe, me, ce) && isWithinCE(qe, ms, ce))
          					 || (!isWithinCE(qs, ms, ce) && isWithinCE(qs, me, ce) && !isWithinCE(qe, me, ce) && !isWithinCE(qe, ms, ce))){

                      vector< SVal > cs_primes = determineCase(delta, cs, candidateMatches, map_entries);
                      if( find( cs_primes.begin(), cs_primes.end(), cs ) == cs_primes.end())
                          cs_primes.push_back( cs );
                      for(auto cs_prime : cs_primes){
                          cout << "in ... " << cs_prime->get_id() << "\n";
                          if(cs_prime && isPotential(delta, cs_prime) && withinAngle(delta, cs_prime)){
                              cout << "out ... " << cs_prime->get_id() << "\n";
            		              createMatches(delta, cs_prime, rel, candidateMatches, map_entries);
                          }
                      }
          				}
  				        //query stub.
          				else if((isWithinCE(qs, ms, ce) && isWithinCE(qe, ms, ce) && !isWithinCE(qs, me, ce) && !isWithinCE(qe, me, ce))
          					 || (isWithinCE(qs, me, ce) && isWithinCE(qe, me, ce) && !isWithinCE(qs, ms, ce) && !isWithinCE(qe, ms, ce))){
          						 // if(delta->get_id() == 1207 && cs->get_id() == 3369)
          						 // 	 cout << "not As expected ! \n";
          				}
  				        //stub in the target database.


                  else {
            					// if(delta->get_id() == 1207 && cs->get_id() == 3369)
            						// cout << "Case 7. \n";

            					if(!isPotential(delta, cs))
            						continue;

            					// if(delta->get_id() == 1207 && cs->get_id() == 3369)
                          // cout << "Wrong place in populate() ... 5" << "\n";

                      double score = similarityScore(delta, cs);
                      if(score != 0.0){
              						Geometry* ls = global_factory->createLineString(formCA(cs->get_intermediate_nodes()));
              						Geometry* buff = delta->get_buffer();
              						if(!buff->intersects(ls))
              							continue;
              						bool flag = buff->covers(ls);
                          boost::shared_ptr<Tuple> tuple(new Tuple(cs, score, flag));
                          rel->insert(tuple);
                      }

                      vector<boost::shared_ptr<Segment_creator>> cs_primes = determineCase(delta, cs, candidateMatches, map_entries);
                      for(auto cs_prime : cs_primes){
                          // cout << "in ... " << cs_prime->get_id() << "\n";
                          if(cs_prime && isPotential(delta, cs_prime)){
                            createMatches(delta, cs_prime, rel, candidateMatches, map_entries);
                            // cout << "out ... " << cs_prime->get_id() << "\n";
                          }
                      }
                  }
              }
              relations.insert({it.first, rel});
          }

          cout << "\n\nstarting with one more pass ... \n\n";
          if(full_pass)
              oneMorePass(querys, relations);

          cout << "\n\nFull pass done ... \n\n";

      		for(auto& it : querys){
      			  RVal rel = relations.find(it.first)->second;
              cout << "duplicate entries size ... " << it.first << "\t" << rel->get_entries().size() << "\n";
              cout << "all entries ... \n\n";
              for( auto en : rel->get_entries() )
                  cout << en->segment()->get_first()->get_id() << "\t" << en->segment()->get_second()->get_id() << "\n";
              removeDuplicates(rel->get_entries());
              sort(rel->get_entries().begin(), rel->get_entries().end(), sortRelationTuples);
      		}

          cout << "removed duplicate queries ... \n\n";
	  }

    void pruneCandidates( map< SKey, RVal > &relations ){
        for( auto it : relations ){
            auto& entries = it.second->get_entries();
            if( entries.size() == 1 )
                continue;
            size_t i = 0;
            for(i = 0; i < entries.size()-1; i++){
                for(size_t j = i+1; j<entries.size()-1; j++){
                    if(subsetOf( entries[ j ], entries[ i ] ) && entries[ j ]->score() <= 0.8 * entries[ i ]->score() ){
                      entries.erase(entries.begin()+j--);
                    }
                }
            }
        }
    }

    bool subsetOf( TVal en1, TVal en2 ){
        if( en1->segment()->get_buffer()->intersection( en2->segment()->get_buffer() )->getArea()/min( en1->segment()->get_buffer()->getArea(), en2->segment()->get_buffer()->getArea() ) >= 0.9999 )
            return true;
        return false;
    }

    SVal polygonMatch( SVal delta, list< SKey > candidateMatches ){
      map< SKey, map< SKey, SVal >> candidates;

      for( auto it : candidateMatches ){
          SVal seg = segments[ 1 ][ it ];
          candidates[ seg->get_first()->get_id() ][ seg->get_second()->get_id() ] = seg;
          candidates[ seg->get_second()->get_id() ][ seg->get_first()->get_id() ] = seg;
      }

      for( auto it1 = candidateMatches.begin(); it1 != candidateMatches.end(); it1++ ){
          SVal seg = segments[ 1 ][ *it1 ];
          cout << "starting with it1 ... " << seg->get_first()->get_id() << "\t" << seg->get_second()->get_id() << "\n";
          if( notPartOfPolygon( seg, candidates ) )
              continue;
          vector< SVal > joined = polygonMatchUtil( seg, candidates );
          size_t i = 0;
          while( i < joined.size() ){
              auto it = joined[ i ];
              if( (it)->get_first()->get_id() == (it)->get_second()->get_id() )
                  return it;
              vector< SVal > joined_seg = polygonMatchUtil( it, candidates );
              if( joined_seg.size() > 0 )
                joined.insert( joined.end(), joined_seg.begin(), joined_seg.end() );
              i++;
          }
      }
      return NULL;
    }


    bool notPartOfPolygon( SVal seg, map< SKey, map< SKey, SVal>> candidates ){
        if( candidates[ seg->get_first()->get_id() ].size() <= 1 || candidates[ seg->get_second()->get_id() ].size() <= 1 )
            return true;
        return false;
    }


    vector< SVal > polygonMatchUtil( SVal seg1, map< SKey, map< SKey, SVal >> candidates ){
        cout << "joining with seg1 ... " << seg1->get_first()->get_id() << "\t" << seg1->get_second()->get_id() << "\n";
        vector< SVal > joined_result = {};

        map< SKey, SVal > firsts = candidates[ seg1->get_first()->get_id() ];
        for( auto it : firsts ){
            SVal seg2 = it.second;
            if( same_entries( seg1, seg2 ) || in_parent( seg2, seg1 ) || notPartOfPolygon( seg2, candidates ) )
                continue;
            SVal joined = concatUtil( seg1, seg2 );
            if( joined )
                joined_result.push_back( joined );
        }

        map< SKey, SVal > seconds = candidates[ seg1->get_second()->get_id() ];
        size_t i = 0, size = joined_result.size(), cnt=0;
        while( i < size ){
            auto seg = joined_result[ i-cnt ];
            for( auto it : seconds ){
                SVal seg2 = it.second;
                if( same_entries( seg, seg2 ) || in_parent( seg2, seg ) || notPartOfPolygon( seg2, candidates ) )
                    continue;
                SVal joined = concatUtil( seg, seg2 );
                if( joined ){
                    joined_result.push_back( joined );
                    joined_result.erase( joined_result.begin() + i - cnt );
                    cnt++;
                }
            }
            i++;
        }

        return joined_result;


        // for( auto seg2 : candidates ){
        //     if( same_entries( seg1, seg2 ) || in_parent( seg2, seg1 ) )
        //         continue;
        //     SVal joined;
        //     if( seg1->get_second()->get_id() == ( seg2 )->get_first()->get_id() ){
        //         joined = concat( seg1, seg2, true );
        //         cout << "joined ... " << joined->get_first()->get_id() << "\t" << joined->get_second()->get_id() << "\n";
        //         // break;
        //     }
        //     else if( ( seg2 )->get_second()->get_id() == seg1->get_first()->get_id() ){
        //         joined = concat( seg2, seg1, true );
        //         cout << "joined ... " << joined->get_first()->get_id() << "\t" << joined->get_second()->get_id() << "\n";
        //
        //         // break;
        //     }
        //     else if( seg1->get_first()->get_id() == ( seg2 )->get_first()->get_id() ||
        //              seg1->get_second()->get_id() == ( seg2 )->get_second()->get_id() ){
        //         joined = concat( seg1, seg2, false );
        //         cout << "joined ... " << joined->get_first()->get_id() << "\t" << joined->get_second()->get_id() << "\n";
        //
        //         // break;
        //     }
        //     if( joined )
        //         joined_result.push_back( joined );
        // }
        // return joined_result;
    }


    SVal concatUtil( SVal seg1, SVal seg2 ){
        SVal joined;
        if( seg1->get_second()->get_id() == ( seg2 )->get_first()->get_id() ){
            joined = concat( seg1, seg2, true );
            cout << "joined ... " << joined->get_first()->get_id() << "\t" << joined->get_second()->get_id() << "\n";
        }
        else if( ( seg2 )->get_second()->get_id() == seg1->get_first()->get_id() ){
            joined = concat( seg2, seg1, true );
            cout << "joined ... " << joined->get_first()->get_id() << "\t" << joined->get_second()->get_id() << "\n";
        }
        else if( seg1->get_first()->get_id() == ( seg2 )->get_first()->get_id() ||
            seg1->get_second()->get_id() == ( seg2 )->get_second()->get_id() ){
            joined = concat( seg1, seg2, false );
            cout << "joined ... " << joined->get_first()->get_id() << "\t" << joined->get_second()->get_id() << "\n";
        }
        return joined;
    }


    void removeDummySegments( map<SKey, RVal> &relations, size_t percentile ){
        vector< double > scores;
        for( auto it : relations ){
            vector< TVal > entries = it.second->get_entries();
            // cout << "entries size ... " << entries.size() << "\n";
            for( auto entry : entries ){
                if( entry->score() > 0 ){
                    scores.push_back( entry->score() );
                }
            }
        }
        sort(scores.begin(), scores.end(), score_sort);
        size_t size = scores.size();
        cout << "size of scores ... " << size << "\n";
        double index = (double)(100 - percentile)/100;
        cout << "index ... " << (int)(( index * size ) - 1) << "\n";
        double percentile_score = scores[ (int)((index * size ) - 1) ];
        cout << "percentile score ... " << percentile_score << "\n";
        size_t dummy_removed = 0, single_matches = 0;
        for( auto it : relations ){
            auto& all_entries = it.second->get_entries();
            if ( all_entries[ 0 ]->score() > percentile_score ){
              all_entries.erase( all_entries.begin() + all_entries.size() - 1);
              dummy_removed++;
            }
            else if(all_entries.size() == 1)
              single_matches++;
        }
        cout << "number of segments where dummy removed ... " << dummy_removed+single_matches << "\n";
    }


    void finalWrite(map<SKey, SVal> &querys, map<SKey, RVal> &relations){
		int total_count = 0;
		for(auto it : querys){
			// if(splitted_segments.find(it.first) != splitted_segments.end())
			// 	continue;
			if(it.first == 1097)
				cout << "Sanity check : " << it.second->get_intermediate_nodes().size() << "\n";
			std::cout.precision(5);
			ofstream test_delta("./ref_data/"+to_string(it.first)+".txt");
			test_delta << "LINESTRING (" << std::fixed;
			CoordinateArraySequence* cl = new CoordinateArraySequence();
			cl = formCA(segments[0].find(it.first)->second->get_intermediate_nodes());
			size_t count=0;
			size_t size = cl->getSize();
			while(count < size){
				Coordinate c = cl->getAt(count);
				if(count == 0)
					test_delta << c.x << " " << c.y;
				else
					test_delta << ", " << c.x << " " << c.y;
				count++;
			}
			test_delta << ")";
			if(relations.find(it.first) == relations.end())
				continue;
            RVal rel = relations.find(it.first)->second;
            removeDuplicates(rel->get_entries());
            sort(rel->get_entries().begin(), rel->get_entries().end(), sortRelationTuples);
            total_count += rel->get_entries().size();
            testFile << it.first << "\n";

            for(auto& item : rel->get_entries()){
                testFile << item->segment()->get_id() << "\t" << item->segment()->get_first()->get_id() << "\t" << item->segment()->get_second()->get_id() << "\t" << item->is_contained() << "\t" << item->score() << "\n";
				// if(it.first == 1189 && item->segment()->get_id() == -1649){
					ofstream test_cs("./tar_data/"+to_string(item->segment()->get_id())+".txt");
					test_cs << "LINESTRING (" << std::fixed;
					CoordinateArraySequence* cl = new CoordinateArraySequence();
					cl = formCA(item->segment()->get_intermediate_nodes());
					size_t count=0;
					size_t size = cl->getSize();
					while(count < size){
						Coordinate c = cl->getAt(count);
						if(count == 0)
							test_cs << c.x << " " << c.y;
						else
							test_cs << ", " << c.x << " " << c.y;
						count++;
					}
					test_cs << ")";
				// }
			}
            testFile << "---------\n";

			// if(it.first == 1189){
			//
			// }
        }

        cout << "Interpolated Nodes: " << (modified_nodeId +1000) << "\n\n";
		    cout << "Size of reference dataset : " << querys.size() << "\n";
        cout << "Average number of candidate matches per relation with CE 20: " << (total_count)/(segments[0].size()) << "\n\n";
        // return relations;
	}


  	void createMatches(SVal delta, SVal cs, RVal &rel, list<SKey> candidateMatches, map< SKey, map< SKey, vector< SVal >>> &map_entries){

        if( delta->get_first()->get_id() == delta->get_second()->get_id()){
            double score = similarityScore(delta, cs);
            if(score <= 0.05)
              return;
            Geometry* ls = global_factory->createLineString(formCA(cs->get_intermediate_nodes()));
            Geometry* buff = delta->get_buffer();
            if(!buff->intersects(ls))
              return;
            bool flag = buff->covers(ls);
            boost::shared_ptr<Tuple> tuple(new Tuple(cs, score, flag));
            rel->insert(tuple);
            return;
        }
    		vector<SVal> results = handleStubs(delta, cs, candidateMatches, map_entries);

    		for(auto& r : results){
            // if( map_entries.find( r->get_first()->get_id() ) != map_entries.end() ){
            //     if( map_entries[ r->get_first()->get_id() ].find( r->get_second()->get_id() ) != map_entries[ r->get_first()->get_id() ].end() )
            //         continue;
            // }
            cout << "in createMatches ... " << r->get_id() << "\t" << r->get_first()->get_id() << "\t" << r->get_second()->get_id() << "\n";

            // if( cs->get_id() == -3527 ){
            //     for( auto it : cs->parents() )
            //         cout << "parent of -3527 ... " << it << "\n";
            // }

      			double score = similarityScore(delta, r);
      			if(score <= 0.05)
      				continue;
      			Geometry* ls = global_factory->createLineString(formCA(r->get_intermediate_nodes()));
      			Geometry* buff = delta->get_buffer();
      			if(!buff->intersects(ls))
      				continue;
      			bool flag = buff->covers(ls);
      			boost::shared_ptr<Tuple> tuple(new Tuple(r, score, flag));
      			rel->insert(tuple);
            if(r->get_length() <= 4*ce){
                relaxed_joins.insert({r->get_first()->get_id(), r->get_second()->get_id()});
                relaxed_joins.insert({r->get_second()->get_id(), r->get_first()->get_id()});
            }
    		}
	   }


    void splitQuery(map<SKey, SVal> &querys, map<SKey, SVal> splitted_querys, map<SKey, RVal> &relations){
    		for(auto& query : splitted_querys ){
      			SVal delta = query.second;

      			if(relations.find(query.first) == relations.end())
      				continue;
      			RVal rel = relations.find(query.first)->second;
      			if(rel->get_entries().size() == 1)
      				continue;

      			SVal cs = rel->get_entries()[0]->segment();

            if( cs->get_id() == 0 )
                continue;

      			tn qs = delta->get_first(), qe = delta->get_second(), ms = cs->get_first(), me = cs->get_second();
      			if((isWithinCE(qs, ms, ce) && isWithinCE(qe, me, ce))
      			|| (isWithinCE(qs, me, ce) && isWithinCE(qe, ms, ce))
      			|| (!isWithinCE(qs, ms, ce) && !isWithinCE(qs, me, ce) && !isWithinCE(qe, me, ce) && !isWithinCE(qe, ms, ce)))
      				continue;

      			tn source_node; bool flag = 0;
      			if(isWithinCE(qs, ms, ce) && !isWithinCE(qe, me, ce)){
      				source_node = me;
      				flag = 0;
      			}
      			else if(isWithinCE(qs, me, ce) && !isWithinCE(qe, ms, ce)){
      				source_node = ms;
      				flag = 0;
      			}
      			else if(!isWithinCE(qs, ms, ce) && isWithinCE(qe, me, ce)){
      				source_node = ms;
      				flag = 1;
      			}
      			else if(!isWithinCE(qs, me, ce) && isWithinCE(qe, ms, ce)){
      				source_node = me;
      				flag = 1;
      			}

      			vector<SVal> segs;
      			bool found = false;
      			// map<size_t, vector<Coordinates>> nearest = nearestNodesForSplit(source_node, delta->get_intermediate_nodes());
            map<size_t, vector<Coordinates>> pair;

            vector<IntNode> cl = delta->get_intermediate_nodes();
            size_t pos=0;

            pos = getNearestNodesBuffer( source_node, cl );
            cout << "\n\npos in interpolate ... " << pos << "\n\n";
            if( pos < cl.size()-1 ){
                vector< Coordinates > coords;
                coords.push_back( cl[ pos ].get_location() );
                coords.push_back( cl[ pos+1 ].get_location() );
                pair.insert( { pos, coords } );
            }
            else
                continue;
            for(auto& item : pair){
        				size_t pos = item.first;
        				Coordinates c1 = item.second[0];
        				Coordinates c2 = item.second[1];
        				Coordinates c = getInterpolatedNodeCoordinates(source_node->get_location(), c1, c2);
        				vector<IntNode> ci = delta->get_intermediate_nodes();
        				vector<vector<IntNode>> split_segments = insert(source_node, c, ci, pos, 1);
        				if(split_segments[flag].size() < 2)
        	          continue;
  	            boost::shared_ptr<Terminal_node> tnode(new Terminal_node(modified_nodeId, c, 1));
  	            --modified_nodeId;
  	            segs = getInterpolatedSegment(source_node, delta, split_segments, tnode);
        				if(!segs[flag])
        					break;
        				else
        					found = true;
      			}
      			if(!found)
      				  continue;

            if( segs[ 0 ]->get_length() == 0 || segs[ 1 ]->get_length() == 0 )
                continue;

      			querys.find(query.first)->second = segs[flag];

      			segs[flag]->set_id(query.first);

      			// splitted_segments.insert({query.first, segs[flag]});

      			segments[0][query.first] = segs[flag];

            recomputeScores( query.first, relations );

            segs[1-flag]->set_id(seg_id++);
            splitted_segments.insert({segs[1-flag]->get_id(), segs[1-flag]});
            querys.insert({segs[1-flag]->get_id(), segs[1-flag]});
            // splitted_querys.insert({segs[1-flag]->get_id(), segs[1-flag]});
            cout << "Inserted : " << delta->get_id() << " " << cs->get_id() << " " << segs[1-flag]->get_id() << "\n";
            cout << "length of split segments ... " << segs[ 0 ]->get_length() << "\t" << segs[ 1 ]->get_length() << "\n";
    		}
  	}


    void recomputeScores(SKey id, map<SKey, RVal> &relations){

    			TVal tuple = relations[ id ]->get_entries()[0];
    			SVal delta = segments[ 0 ][ id ];
    			SVal cs = tuple->segment();

    			double score = similarityScore(delta, cs);
    			if(score == 0.0)
    				return;
    			Geometry* ls = global_factory->createLineString(formCA(cs->get_intermediate_nodes()));
    			Geometry* buff = delta->get_buffer();
    			if(!buff->intersects(ls))
    				return;
    			bool flag = buff->covers(ls);
    			tuple->set_contained(flag);
    			tuple->set_score(score);

  	}


  	void recomputeScores(map<SKey, SVal> &querys, map<SKey, RVal> &relations){
    		for(auto& it : querys){
    			TVal tuple = relations.find(it.first)->second->get_entries()[0];
    			SVal delta = it.second;
    			SVal cs = tuple->segment();

    			double score = similarityScore(delta, cs);
    			if(score == 0.0)
    				continue;
    			Geometry* ls = global_factory->createLineString(formCA(cs->get_intermediate_nodes()));
    			Geometry* buff = delta->get_buffer();
    			if(!buff->intersects(ls))
    				continue;
    			bool flag = buff->covers(ls);
    			tuple->set_contained(flag);
    			tuple->set_score(score);
    		}
  	}


  	bool withinAngle(SVal delta, SVal cs){
    		Coordinates c1 = delta->get_first()->get_location(), c2 = delta->get_second()->get_location(),
    					      c3 = cs->get_first()->get_location(), c4 = cs->get_second()->get_location();
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

        bool stub = delta->get_length() <= ce || cs->get_length() <= ce;

    		// if(delta->get_id() == 1018)
    		// cout << "Angle Difference : " << cs->get_id() << "\t" << cs->get_first()->get_id() << "\t" << cs->get_second()->get_id() << "\t" << abs(atan1 - atan2) << "\t" << PI-abs(atan1 - atan2) << "\t" << PI/6 << "\n";
    		if(stub && min(abs(atan1 - atan2), PI-abs(atan1 - atan2)) > PI/6)
    			// return true;
    		return false;
    		// if(delta->get_length() <= ce && cs->get_length() <= ce)
    		// 	return false;
    		Geometry* g1 = delta->get_buffer();
    		Geometry* g2 = cs->get_buffer();
    		Geometry* g = g1->intersection(g2);
    		double value = g->getArea()/min(g1->getArea(), g2->getArea());
        // cout << "overlap value ... " << value << "\n";
        if( stub && value > 0.8)
            return true;
        if( !stub && value > 0.5)
            return true;
    		// if(value > 0.8)
    		// 	return true;
    		return false;
  	}

    bool endPointsWithinCE( SVal delta, SVal cs ){
        tn qs = delta->get_first(), qe = delta->get_second(), ms = cs->get_first(), me = cs->get_second();
        if( isWithinCE(qs, ms, ce) && isWithinCE(qe, me, ce) )
            return true;
        // if( isWithinCE(qs, ms, ce) && ! isWithinCE(qe, me, ce) )
        //     return true;
        // if( isWithinCE(qe, me, ce) && ! isWithinCE(qs, ms, ce) )
        //     return true;

        if( isWithinCE(qs, me, ce) && isWithinCE(qe, ms, ce) )
            return true;
        // if( isWithinCE(qs, me, ce) && ! isWithinCE(qe, ms, ce) )
        //     return true;
        // if( ! isWithinCE(qs, me, ce) && isWithinCE(qe, ms, ce) )
        //     return true;
        return false;
    }


  	void filterCandidates(SVal delta, list<SKey> &candidateMatches){
    		vector<size_t> toErase;
    		size_t i = 0;
    		for(auto& it : candidateMatches){
      			SVal cs = segments[1].find(it)->second;
            Geometry* ls = global_factory->createLineString(formCA(cs->get_intermediate_nodes()));
      			Geometry* buff = delta->get_buffer();
      			if(!buff->intersects(ls)){
                toErase.push_back(i);
            }
            // else if( cs->get_length() <= ce && !withinAngle( delta, cs ) )
            //     toErase.push_back(i);
      			// if(!isPotentialCandidate(delta, cs)){
      			// 	toErase.push_back(i);
      			// }
      			++i;
    		}
    		size_t cnt = 0;
    		auto it = candidateMatches.begin();
    		for(i = 0; i < toErase.size(); i++){
      			it = candidateMatches.begin();
      			advance(it, toErase[i] - cnt);
      			candidateMatches.erase(it);
      			cnt++;
    		}
  	}


    bool isPotentialCandidate(SVal delta, SVal cs){
    		if(delta->is_roundabout() && cs->is_roundabout())
    			return true;
    		if(delta->is_roundabout() || cs->is_roundabout())
    			return false;
    		if(!withinAngle(delta, cs) && cs->get_length() <= ce)
    			return false;
        if(withinAngle(delta, cs))
            return true;

    		if(delta->get_id() == 1571 && cs->get_id() == 3271)
    			cout << "Right place in potential candidate for 1571. \n";


    		Coordinates c1 = delta->get_first()->get_location(), c2 = delta->get_second()->get_location(),
    					c3 = cs->get_first()->get_location(), c4 = cs->get_second()->get_location();

    		// if((isWithinCE(delta->get_first(), cs->get_first(), ce) && isWithinCE(delta->get_second(), cs->get_second(), ce))
    		// || (isWithinCE(delta->get_first(), cs->get_second(), ce) && isWithinCE(delta->get_second(), cs->get_first(), ce)))
    		// 	return true;

    		if(delta->get_id() == 1571 && cs->get_id() == 3271)
    			cout << "Right place in potential candidate for 1571. \n";

    		double d13 = distance(c1, c3), d23 = distance(c2, c3),
    			   d14 = distance(c1, c4), d24 = distance(c2, c4);

    	   if(delta->get_id() == 1571 && cs->get_id() == 3271)
    			cout << "Distances : " << d13 << "\t" << d23 << "\t" << d14 << "\t" << d24 << "\n";

    		if(d13 < d23){
    			if(d14 < d24){
    				if(d13 < d14){
    					return samePattern(delta, cs, c1, c2, c3, c4);
    				}
    				else{
    					return samePattern(delta, cs, c1, c2, c4, c3);
    				}
    			}
    			else{
    				return samePattern(delta, cs, c1, c2, c3, c4);
    			}
    		}
    		else{
    			if(d14 < d24){
    				return samePattern(delta, cs, c2, c1, c3, c4);
    			}
    			else{
    				if(d23 < d24){
    					return samePattern(delta, cs, c2, c1, c3, c4);
    				}
    				else{
    					return samePattern(delta, cs, c2, c1, c4, c3);
    				}
    			}
    		}
    		return true;
	  }


  	bool isPotential(SVal delta, SVal cs){
        return true;
    		tn ds = delta->get_first(), de = delta->get_second(), css = cs->get_first(), cse = cs->get_second();
    		Coordinates c1 = ds->get_location(), c2 = de->get_location(), cs1 = css->get_location(), cs2 = cse->get_location();
    		if(isWithinCE(ds, css, ce)
    		&& !isWithinCE(ds, cse, ce)
    		&& !isWithinCE(de, css, ce)
    		&& !isWithinCE(de, cse, ce)){

    			if(samePattern(delta, cs, c1, c2, cs1, cs2))
    				return true;
    			return false;
    		}
    		if(isWithinCE(de, cse, ce)
    		&& !isWithinCE(de, css, ce)
    		&& !isWithinCE(ds, css, ce)
    		&& !isWithinCE(ds, cse, ce)){
    			if(samePattern(delta, cs, c1, c2, cs1, cs2))
    				return true;
    			return false;
    		}
    		if(isWithinCE(de, css, ce)
    		&& !isWithinCE(de, cse, ce)
    		&& !isWithinCE(ds, css, ce)
    		&& !isWithinCE(ds, cse, ce)){
    			if(samePattern(delta, cs, c2, c1, cs1, cs2))
    				return true;
    			return false;
    		}
    		if(isWithinCE(ds, cse, ce)
    		&& !isWithinCE(ds, css, ce)
    		&& !isWithinCE(de, css, ce)
    		&& !isWithinCE(de, cse, ce)){
    			if(samePattern(delta, cs, c2, c1, cs1, cs2))
    				return true;
    			return false;
    		}
    		return true;
  	}


    void oneMorePass(map<SKey, SVal> querys, map<SKey, RVal> &relations){
        for(auto it : querys){
            SVal delta = it.second;
            RVal rel = relations[ it.first ];
            if( it.first == 1107 )
                cout << "in one more pass for " << it.first << " ... \n";
            std::list<long long int> candidateMatches = rtree_search( delta );
			      filterCandidates( delta, candidateMatches );
            if( it.first == 1107 )
                cout << "size of matches for 1107 in one more pass ... \n";
            for( auto& i : candidateMatches ){
                SVal cs = segments[ 1 ][ i ];
        				if(delta->get_id() == 1107 && cs->get_id() == 1469){
          					bool check = modified_nodes.find(delta->get_second()->get_id()) != modified_nodes.end();
          					cout << "1107 : " << check << "\n";
          					check = modified_nodes.find(delta->get_first()->get_id()) != modified_nodes.end();
          					cout << "1107 : " << check << "\n";
        				}
                if(modified_nodes.find( delta->get_first()->get_id() ) != modified_nodes.end() || modified_nodes.find(delta->get_second()->get_id()) != modified_nodes.end()){
					           if(delta->get_id() == 1107 && cs->get_id() == 1469 )
						               checking << "Inside 1 got inside ! \n";
          					if(modified_nodes[ delta->get_first()->get_id() ].find(i) != modified_nodes[ delta->get_first()->get_id() ].end()){
          						// if(delta->get_id() == 1647 && cs->get_id() == 2630)
          							// checking << "Inside 1 passed within() ... " << delta->get_id() << "\n";
            						tn node1 = modified_nodes[ delta->get_first()->get_id() ][ i ];
                        tn node2 = cs->get_second();
                        SVal seg;
                        if( end_points_segments.find( node1->get_id() ) != end_points_segments.end() ){
                            if( end_points_segments[ node1->get_id() ].find( node2->get_id() ) != end_points_segments[ node1->get_id() ].end() ){
                                // cout << "checkpoint1 ... \n";
                                seg = end_points_segments[ node1->get_id() ][ node2->get_id() ];
                            }
                        }
                        else if( end_points_segments.find( node2->get_id() ) != end_points_segments.end() ){
                            if( end_points_segments[ node2->get_id() ].find( node1->get_id() ) != end_points_segments[ node2->get_id() ].end() ){
                                // cout << "checkpoint2 ... \n";
                                seg = end_points_segments[ node2->get_id() ][ node1->get_id() ];
                            }
                        }
                        if( seg ){
                            if(isPotentialCandidate( delta, seg ) && withinAngle( delta, seg ) && endPointsWithinCE( delta, seg ) ){
                                double score = similarityScore( delta, seg );
                                if(score != 0.0){
                                    Geometry* ls = global_factory->createLineString( formCA( seg->get_intermediate_nodes() ) );
                                    Geometry* buff = delta->get_buffer();
                                    if(!buff->intersects(ls))
                                        continue;
                                    // bool print = (isWithinCE(delta->get_first(), segs[1]->get_first(), ce) && isWithinCE(delta->get_second(), segs[1]->get_second(), ce)) ||
                                    // 			 (isWithinCE(delta->get_first(), segs[1]->get_second(), ce) && isWithinCE(delta->get_second(), segs[1]->get_first(), ce));
                                    // checking << delta->get_id() << " " << print << "\n";
                                    bool flag = buff->covers(ls);
                                    boost::shared_ptr<Tuple> tuple(new Tuple(seg, score, flag));
                                    rel->insert(tuple);
                                 }
                            }
                        }
                    }

                    if( modified_nodes[ delta->get_second()->get_id() ].find(i) != modified_nodes[ delta->get_second()->get_id() ].end() ){
						        //if(delta->get_id() == 1282)
					              checking << "Inside 2 ... " << delta->get_id() << " " << i << " " << delta->get_first()->get_id() << " " << delta->get_second()->get_id() << "\n";
                        tn node1 = modified_nodes[ delta->get_second()->get_id() ][ i ];
                        tn node2 = cs->get_first();
                        SVal seg;

                        if( end_points_segments.find( node1->get_id() ) != end_points_segments.end() ){
                            if( end_points_segments[ node1->get_id() ].find( node2->get_id() ) != end_points_segments[ node1->get_id() ].end() ){
                                // cout << "checkpoint3 ... \n";
                                seg = end_points_segments[ node1->get_id() ][ node2->get_id() ];
                            }
                        }
                        else if( end_points_segments.find( node2->get_id() ) != end_points_segments.end() ){
                            if( end_points_segments[ node2->get_id() ].find( node1->get_id() ) != end_points_segments[ node2->get_id() ].end() ){
                                // cout << "checkpoint4 ... \n";
                                seg = end_points_segments[ node2->get_id() ][ node1->get_id() ];
                            }
                        }
                        if( seg ){
                            // if(isPotentialCandidate( delta, seg ) && withinAngle( delta, seg ) && endPointsWithinCE( delta, seg ) ){
                              if(isPotentialCandidate( delta, seg ) && withinAngle( delta, seg ) && endPointsWithinCE( delta, seg ) ){
                                double score = similarityScore( delta, seg );
                                if(score != 0.0){
                                    Geometry* ls = global_factory->createLineString( formCA( seg->get_intermediate_nodes() ) );
                                    Geometry* buff = delta->get_buffer();
                                    if(!buff->intersects(ls))
                                        continue;
                                    // bool print = (isWithinCE(delta->get_first(), segs[1]->get_first(), ce) && isWithinCE(delta->get_second(), segs[1]->get_second(), ce)) ||
                                    // 			 (isWithinCE(delta->get_first(), segs[1]->get_second(), ce) && isWithinCE(delta->get_second(), segs[1]->get_first(), ce));
                                    // checking << delta->get_id() << " " << print << "\n";
                                    bool flag = buff->covers(ls);
                                    boost::shared_ptr<Tuple> tuple(new Tuple(seg, score, flag));
                                    rel->insert(tuple);
                                 }
                            }
                        }
                    }
                }
            }
        }
    }


    vector<SVal> handleStubs(SVal delta, SVal seg1, std::list<SKey> candidateMatches, map< SKey, map< SKey, vector< SVal >>> &map_entries){ //when length of delta is <=ce.
        vector<SVal> results = { seg1 };
        // vector<SVal> results = {};
        // return results;
        queue<SVal> result_copy;
        result_copy.push( seg1 );
        bool flag = true;
        while( result_copy.size() > 0 ){
          // cout << delta->get_id() << "  " << seg1->get_id() << " " << result_copy.size() << "\n";
          SVal seg1 = result_copy.front();
          result_copy.pop();
          // if( delta->get_id() == 1207 )
          //     cout << "handle stubs() seg1 ... " << seg1->get_id() << " " << seg1->get_first()->get_id() << " " << seg1->get_second()->get_id() << "\n";
          // cout << result_copy.size() << "\n";
          // if( delta->get_id() == 1207 ){
          //     cout << "parents ... \n";
          //     for( auto it : seg1->parents())
          //         cout << it << "\n";
          // }
          for( auto j : candidateMatches ){
              SVal seg2 = segments[ 1 ].find( j )->second;

              if(seg2->get_length() > 1.5*ce || contains_parent( seg1->parents(), j )){
                  // bool check = contains_parent( seg1->parents(), j );
                  // if(delta->get_id() == 1207)
                  //     cout << "for this seg2 ... " << seg2->get_length() << " " << 2*ce << " " << check << " " << j << " ^\n";
                  continue;
              }
              // if( in_parent( seg2, seg1 ) )
              //     continue;
              if(seg1->get_second()->get_id() == seg2->get_first()->get_id() && isSegWithinCE(seg2->get_second(), delta)){
                  SVal seg_joined = concat(seg1, seg2, true);
                  results.push_back(seg_joined);
                  result_copy.push( seg_joined );
              }
              else if(seg1->get_first()->get_id() == seg2->get_second()->get_id() && isSegWithinCE(seg1->get_first(), delta)){
                  SVal seg_joined = concat(seg2, seg1, true);
                  results.push_back(seg_joined);
                  result_copy.push( seg_joined );
              }
              else if(seg1->get_first()->get_id() == seg2->get_first()->get_id() && isSegWithinCE(seg2->get_second(), delta)){
                  SVal seg_joined = concat(seg2, seg1, false);
                  results.push_back(seg_joined);
                  result_copy.push( seg_joined );
              }
              else if(seg2->get_second()->get_id() == seg1->get_second()->get_id() && isSegWithinCE(seg2->get_first(), delta)){
                  SVal seg_joined = concat(seg1, seg2, false);
                  results.push_back(seg_joined);
                  result_copy.push( seg_joined );
              }
          }
        }
        return results;
    }

    bool contains_parent( vector< SKey > parents, SKey id ){
        for( auto it : parents ){
            if( it == id )
              return true;
        }
        return false;
    }


    inline bool isSegWithinCE(boost::shared_ptr<Terminal_node> tn, SVal seg2){
        if (isWithinCE( tn, seg2->get_first(), ce) && isWithinCE( tn, seg2->get_second(), ce ) )
            return false;
        return (isWithinCE(tn, seg2->get_first(), ce) || isWithinCE(tn, seg2->get_second(), ce));
    }

    void removeDuplicates(vector<TVal> &entries){
        size_t i = 0;
        for(i = 0; i < entries.size(); i++){
            for(size_t j = i+1; j<entries.size(); j++){
                if(same_entries(entries[i], entries[j])){
                   // if( entries[i]->segment()->get_id() == -5982 && entries[j]->segment()->get_id() == 2745)
                   //    cout << "returned true for -5982 and 2745 .... \n";

                  entries.erase(entries.begin()+j--);
                }
                // else if(entries[i]->segment()->get_id() == -5982 && entries[j]->segment()->get_id() == 2745)
                // cout << "returned false for -5982 and 2745 .... \n";
            }
        }
    }

    void removeDuplicates(vector<SVal> &cs_primes){
        size_t i = 0;
        for(i = 0; i < cs_primes.size(); i++){
            for(size_t j = i+1; j<cs_primes.size(); j++){
                if(same_entries(cs_primes[i], cs_primes[j])){
                   // if( entries[i]->segment()->get_id() == -5982 && entries[j]->segment()->get_id() == 2745)
                   //    cout << "returned true for -5982 and 2745 .... \n";

                  cs_primes.erase(cs_primes.begin()+j--);
                }
                // else if(entries[i]->segment()->get_id() == -5982 && entries[j]->segment()->get_id() == 2745)
                // cout << "returned false for -5982 and 2745 .... \n";
            }
        }
    }


    bool same_entries(TVal en1, TVal en2){
        // double score = similarityScore( en1->segment(), en2->segment() );
        // if( score == 1 )
        //     return true;
        // return false;
        // if( en1->segment()->get_id() == -5982 && en2->segment()->get_id() == 2745){
        //     double sim_score = similarityScore(en1->segment(), en2->segment());
        //     cout << "sim score between -5982 and 2745 is ... " << sim_score << "\t" << (sim_score > 0.99999) << "\n";
        // }
        if( en1->segment()->get_id() == 0 || en2->segment()->get_id() == 0 )
            return false;
        // if(en1->get_first()->get_id() == en2->get_first()->get_id() && en1->get_second()->get_id() == en2->get_second()->get_id())
        //     return true;
        //     cout << "entered in similarityScore ... \n";


        // }
        //
    		// if(en1->get_first()->get_id() == en2->get_second()->get_id() && en1->get_second()->get_id() == en2->get_first()->get_id())
        //     return true;
        //     cout << "entered in similarityScore ... \n";
        //     return similarityScore(en1->segment(), en2->segment()) == 1;
        // }
        return ( similarityScore(en1->segment(), en2->segment()) >= 0.99999 );
    			// return similarityScore(en1->segment(), en2->segment()) == 1;
            // return false;
    }

    bool same_entries( SVal seg1, SVal seg2 ){
        return ( similarityScore(seg1, seg2) >= 0.99999 );
    }


    void refine(vector<TVal> &entries){
        if(entries.size() > 0){
            double max_score = entries.front()->score()/2;
            size_t i = 0;
            for(auto en : entries){
                if(en->score() < max_score && en->segment()->get_id() != 0)
                    break;
                ++i;
            }
            if(i == entries.size())
                return;
            for(auto en = entries.begin() + i; en != entries.end()-1; en++)
                entries.erase(en--);
        }
    }


    vector<SVal> determineCase(SVal delta, SVal cs, list<SKey> candidateMatches, map<SKey, map<SKey, vector<SVal>>> &map_entries){

        SVal cs_prime;
        vector<SVal> cs_primes;

        uint32_t ce1 = delta->getCircularError();
        uint32_t ce2 = cs->getCircularError();
        uint32_t CE = max(ce1, ce2);

        if(isWithinCE(delta->get_first(), cs->get_first(), CE) && isWithinCE(delta->get_second(), cs->get_second(), CE)){

            // cout << "Condition 1 ------" <<  cs->get_id() << "\n";
            cs_primes.push_back(cs);
            return cs_primes;
        }
        if(isWithinCE(delta->get_first(), cs->get_second(), CE) && isWithinCE(delta->get_second(), cs->get_first(), CE)){

            // cout << "Condition 2 ------" <<  cs->get_id() << "\n";
            cs_primes.push_back(cs);
            return cs_primes;
        }
        if((isWithinCE(delta->get_first(), cs->get_first(), CE) && !isWithinCE(delta->get_second(), cs->get_second(), CE))
		    || (isWithinCE(delta->get_first(), cs->get_second(), CE) && !isWithinCE(delta->get_second(), cs->get_first(), CE))
        || (!isWithinCE(delta->get_first(), cs->get_first(), CE) && isWithinCE(delta->get_second(), cs->get_second(), CE))
		    || (!isWithinCE(delta->get_first(), cs->get_second(), CE) && isWithinCE(delta->get_second(), cs->get_first(), CE))){


            // cout << "Condition 3 ------ " <<  cs->get_id() << "\n\n";
            int d = delta->get_length();
            int d_prime = cs->get_length();

            if(d > d_prime){

                // cout << "Right place in determineCase () ... " << cs->get_id() << "\t" << cs->get_first()->get_id() << "\t" << cs->get_second()->get_id() << "\n";

                if(candidateMatches.size() == 0)
                    cs_primes.push_back(cs);
                else{
                    vector<SVal> results = selfJoin(delta, cs, candidateMatches, map_entries);
                    // cout << "Results size: " << results.size() << "\n";
                    for(auto& cs_prime : results){
                        if(cs_prime){
              							// cout << "Result : " << cs_prime->get_id() << " " << cs_prime->get_first()->get_id() << " " << cs_prime->get_second()->get_id() << "\n";
              							cs_primes.push_back(cs_prime);
                            vector<SVal> new_results = determineCase(delta, cs_prime, candidateMatches, map_entries);
                            if( new_results.size() > 0 )
                                cs_primes.insert(cs_primes.end(), new_results.begin(), new_results.end());
                        }
                        else{
                            cs_primes.push_back(cs);
                        }
                    }
                }
            }
            else{

                // cout << "Interpolate in determineCase () for ----------------- " << cs->get_id() << "\t" << cs->get_first()->get_id() << "\t" << cs->get_second()->get_id() << "\n\n";

                // draw_segment( cs );

                tn source;

                if( ( isWithinCE(delta->get_first(), cs->get_first(), ce) && !isWithinCE(delta->get_second(), cs->get_second(), ce) )
                ||  ( isWithinCE(delta->get_first(), cs->get_second(), ce) && !isWithinCE(delta->get_second(), cs->get_first(), ce) ) ){
                    source = delta->get_second();
                }
                else if( ( isWithinCE(delta->get_second(), cs->get_second(), ce) && !isWithinCE(delta->get_first(), cs->get_first(), ce) )
                      || ( isWithinCE(delta->get_second(), cs->get_first(), ce) && !isWithinCE(delta->get_first(), cs->get_second(), ce) ) ){

                    // cout << "result from interpolate ... " << cs_primes.size() << "\n\n";
                    source = delta->get_first();
                }

                cs_primes = interpolate(source, delta, cs);

                  // cout << "in interpolate of determine case ... " << cs_primes.size() << "\n\n";
                  for(auto it : cs_primes){
                    cout << it->get_id() << "\t" << it->get_first()->get_id() << "\t" << it->get_second()->get_id() << "\t" << withinAngle(delta, it) << "\n";
                    // draw_segment( it );
                }
            }
            if( cs_primes.size() == 0 )
                cs_primes.push_back( cs );

            // removeDuplicates( cs_primes );
            return cs_primes;
        }


        if(!isWithinCE(delta->get_first(), cs->get_first(), ce) && !isWithinCE(delta->get_second(), cs->get_second(), ce)
        && !isWithinCE(delta->get_first(), cs->get_second(), ce) && !isWithinCE(delta->get_second(), cs->get_first(), ce)){

            // cout << "Condition 4 --- " << cs->get_id() << "\n\n";
            // if delta is within bounds of cs, interpolate on the two end points of cs.
            if( isWithinBounds( delta, cs ) ){
                // cout << "Case 1 --- " << cs->get_id() << "\n\n";
                vector< SVal > int_result = interpolate( delta->get_first(), delta, cs );
                for( auto it : int_result ){
                    cs_primes = interpolate( delta->get_second(), delta, it );
                }
            }

            // if cs is within bounds of delta, self join on two end points of cs.
            else if( isWithinBounds( cs, delta ) ){
              // cout << "Case 2 --- " << cs->get_id() << "\n\n";
                vector< SVal > int_result = selfJoin( delta, cs, candidateMatches, map_entries );
                // cout << "int result size ... " << int_result.size() << "\n";
                for( auto it : int_result ){
                    // cout << "int result ... " << it->get_id() << "\t" << it->get_first()->get_id() << "\t" << it->get_second()->get_id() << "\n";
                    vector< SVal > new_results = determineCase( delta, it, candidateMatches, map_entries );
                    if( new_results.size() > 0 )
                        cs_primes.insert( cs_primes.end(), new_results.begin(), new_results.end() );
                }
            }

            // if one end point is within buffer and the other is not, interpolate from the one that is within the buffer and self join on another.
            else{
                // cout << "Case 3 --- " << cs->get_id() << "\n\n";
                tn qs = delta->get_first(), qe = delta->get_second();
                vector< SVal > int_result = {};
                if( pointWithinBuffer( qs, cs ) ){
                    int_result = interpolate( qs, delta, cs );
                }
                else if( pointWithinBuffer( qe, cs ) ){
                    int_result = interpolate( qe, delta, cs );
                }
                // cout << "int result size ... " << int_result.size() << "\n";
                for( auto it : int_result ){
                    vector< SVal > new_results = determineCase( delta, it, candidateMatches, map_entries );
                    if( new_results.size() > 0 )
                        cs_primes.insert( cs_primes.end(), new_results.begin(), new_results.end() );
                }
            }
            if( cs_primes.size() == 0 )
                cs_primes.push_back( cs );

            // removeDuplicates( cs_primes );
            return cs_primes;
        }

        if((isWithinCE(delta->get_first(), cs->get_first(), ce) && isWithinCE(delta->get_second(), cs->get_first(), ce) && !isWithinCE(delta->get_first(), cs->get_second(), ce) && !isWithinCE(delta->get_second(), cs->get_second(), ce))
			  || (isWithinCE(delta->get_first(), cs->get_second(), ce) && isWithinCE(delta->get_second(), cs->get_second(), ce) && !isWithinCE(delta->get_first(), cs->get_first(), ce) && !isWithinCE(delta->get_second(), cs->get_first(), ce))){

          // cout << "Condition 5 --- " << cs->get_id() << "\n\n";
  				 Coordinates c1 = delta->get_first()->get_location(), c2 = delta->get_second()->get_location(),
  				 			       c3 = cs->get_first()->get_location(), c4 = cs->get_second()->get_location();

  				 boost::shared_ptr<Terminal_node> source;
           int flag = 0;
           if(isWithinCE(delta->get_first(), cs->get_first(), ce) && isWithinCE(delta->get_second(), cs->get_first(), ce)){
  					 flag = 0;
  					 if(abs(c1.x - c2.x) > abs(c1.y - c2.y)){
    						 if(c1.x < c2.x && c3.x < c4.x){
    							 source = delta->get_second();
    							 flag = 0;
    						 }
    						 else if(c1.x < c2.x && c3.x > c4.x){
    							 source = delta->get_first();
    							 flag = 1;
    						 }
    						 else if(c1.x > c2.x && c3.x > c4.x){
    							 source = delta->get_first();
    							 flag = 1;
    						 }
    						 else if(c1.x > c2.x && c3.x < c4.x){
    							 source = delta->get_second();
    							 flag = 0;
    						 }
  					 }
  					 else{
    						 if(c1.y < c2.y && c3.y < c4.y){
    							 source = delta->get_second();
    							 flag = 0;
    						 }
    						 else if(c1.y < c2.y && c3.y > c4.y){
    							 source = delta->get_first();
    							 flag = 1;
    						 }
    						 else if(c1.y > c2.y && c3.y > c4.y){
    							 source = delta->get_first();
    							 flag = 1;
    						 }
    						 else if(c1.y > c2.y && c3.y < c4.y){
    							 source = delta->get_second();
    							 flag = 0;
    						 }
    					}
           }
  				 else if(isWithinCE(delta->get_first(), cs->get_first(), ce) && isWithinCE(delta->get_second(), cs->get_first(), ce)){
    					 flag = 0;
    					 if(abs(c1.x - c2.x) > abs(c1.y - c2.y)){
    						 if(c1.x < c2.x && c3.x < c4.x){
    							 source = delta->get_second();
    							 flag = 0;
    						 }
    						 else if(c1.x < c2.x && c3.x > c4.x){
    							 source = delta->get_first();
    							 flag = 1;
    						 }
    						 else if(c1.x > c2.x && c3.x > c4.x){
    							 source = delta->get_first();
    							 flag = 1;
    						 }
    						 else if(c1.x > c2.x && c3.x < c4.x){
    							 source = delta->get_second();
    							 flag = 0;
    						 }
    						 else
    						 	return cs_primes;
    					 }
    					 else{
    						 if(c1.y < c2.y && c3.y < c4.y){
    							 source = delta->get_second();
    							 flag = 0;
    						 }
    						 else if(c1.y < c2.y && c3.y > c4.y){
    							 source = delta->get_first();
    							 flag = 1;
    						 }
    						 else if(c1.y > c2.y && c3.y > c4.y){
    							 source = delta->get_first();
    							 flag = 1;
    						 }
    						 else if(c1.y > c2.y && c3.y < c4.y){
    							 source = delta->get_second();
    							 flag = 0;
    						 }
    						 else
    						 	return cs_primes;
    					 }
            }
            cs_primes = interpolate(source, delta, cs );
  			  }
        // removeDuplicates( cs_primes );
        return cs_primes;
    }


    SKey parentSegmentToBeInterpolated( tn source_node, SVal cs ){
        SKey parent_cs = cs->get_id();
        if( cs->get_id() >= 0 ){
          // cout << "\nParent segment to be interpolated for source node id ... " << source_node->get_id() << "\t" << parent_cs << "\n\n";
          return parent_cs;
        }


        Coordinates c_mcp = lonlat_to_mercator( source_node->get_location() );
        Geometry* point_buffer = global_factory->createPoint(Coordinate(pm->makePrecise(c_mcp.x), pm->makePrecise(c_mcp.y)))->buffer( ce );
        double overlap = 0;

        for( auto it : cs->parents() ){
            if( it < 0 )
                continue;
            SVal csp = segments[ 1 ][ it ];
            double o = csp->get_buffer()->intersection( point_buffer )->getArea();
            if( o > overlap ){
                overlap = o;
                parent_cs = it;
            }
        }
        // cout << "\nParent segment to be interpolated for source node id ... " << source_node->get_id() << "\t" << parent_cs << "\n\n";
        return parent_cs;
    }

    SVal joinConcatUtil( SVal base, SVal other, map<SKey, map<SKey, vector<SVal>>> map_entries ){
        if( base->get_second()->get_id() == other->get_first()->get_id() ){
          if( map_entries.find( base->get_first()->get_id() ) != map_entries.end() &&
              map_entries[ base->get_first()->get_id() ].find( other->get_second()->get_id() ) != map_entries[ base->get_first()->get_id() ].end() )
                return NULL;
          // cout << "case 1 in joinConcatUtil ... \n";
          return concat( base, other, true );
        }

        if( other->get_second()->get_id() == base->get_first()->get_id() ){
          if( map_entries.find( other->get_first()->get_id() ) != map_entries.end() &&
              map_entries[ other->get_first()->get_id() ].find( base->get_second()->get_id() ) != map_entries[ other->get_first()->get_id() ].end() )
                return NULL;
          // cout << "case 2 in joinConcatUtil ... \n";
          return concat( other, base, true );
        }

        if( base->get_second()->get_id() == other->get_second()->get_id() ){
          if( map_entries.find( base->get_first()->get_id() ) != map_entries.end() &&
              map_entries[ base->get_first()->get_id() ].find( other->get_first()->get_id() ) != map_entries[ base->get_first()->get_id() ].end() )
                return NULL;
          // cout << "case 3 in joinConcatUtil ... \n";
          // cout << "base ... " << base->get_first()->get_id() << "\t" << base->get_second()->get_id() << "\n";
          // cout << "other ... " << other->get_first()->get_id() << "\t" << other->get_second()->get_id() << "\n";
          return concat( base, other, false );
        }

        if( base->get_first()->get_id() == other->get_first()->get_id() ){
          if( map_entries.find( base->get_second()->get_id() ) != map_entries.end() &&
              map_entries[ base->get_second()->get_id() ].find( other->get_second()->get_id() ) != map_entries[ base->get_second()->get_id() ].end() )
                return NULL;
          // cout << "case 4 in joinConcatUtil ... \n";
          return concat( base, other, false );
        }
        // cout << "case 5 in joinConcatUtil ... \n";
        return NULL;
    }

    vector<SVal> selfJoin(SVal delta, SVal cs, list<SKey> candidateMatches, map<SKey, map<SKey, vector<SVal>>> &map_entries){

        vector<SVal> joined;

        if((isWithinCE(delta->get_first(), cs->get_first(), ce) && !isWithinCE(delta->get_second(), cs->get_second(), ce))
		      || (isWithinCE(delta->get_second(), cs->get_first(), ce) && !isWithinCE(delta->get_first(), cs->get_second(), ce))){

			        // if(delta->get_id() == 1026)
                  // cout << "Right place in self join 1 ... 1004 " << "\n";
              for(auto& i : candidateMatches){
                  SVal delta_prime = segments[1].find(i)->second;
          				if( same_entries(delta_prime, cs) || in_parent(delta_prime, cs)){
                    // if( cs->get_id() == -1086 && i == 1023 )
                        // cout << "oops wrong place ... \n";
                    continue;
                  }


                  if(cs->get_second() == delta_prime->get_first()){
                    // if( cs->get_id() == -1130 && i == 1022 )
                        // cout << "right place ... \n";
                    // continue;

                      if( map_entries.find( cs->get_first()->get_id() ) != map_entries.end() &&
                          map_entries[ cs->get_first()->get_id() ].find( delta_prime->get_second()->get_id() ) != map_entries[ cs->get_first()->get_id() ].end() )
                            continue;

                      // if(delta->get_id() == 1026)
                      // cout << "Possible concat id : " << i << "\n";
                      SVal joined_seg = concat(cs, delta_prime, true);
                      // if( map_entries.find( joined_seg->get_first()->get_id() ) != map_entries.end() &&
                      //     map_entries[ joined_seg->get_first()->get_id() ].find( joined_seg->get_second()->get_id() ) != map_entries[ joined_seg->get_first()->get_id() ].end() )
                      //       continue;

                      // if( withinAngle( delta, joined_seg ))
            						joined.push_back(joined_seg);
                        insertInMapEntries( joined_seg, map_entries );
    				          // if(delta->get_id() == 1026)
    		                cout << "Joined ... " << "\t" << cs->get_id() << "\t" << delta_prime->get_id() << "\t" << joined_seg->get_id() << "\t" << joined_seg->get_first()->get_id() << "\t" << joined_seg->get_second()->get_id() << "\n";
                      // if(delta->get_id() == 1869 && cs->get_id() == 2944)
                      //     cout << "Joined : " << joined_seg->get_first()->get_id() << " : " << joined_seg->get_second()->get_id() << "\n";
                      // return joined;
                  }
          				else if(cs->get_second() == delta_prime->get_second()){
                    // if(delta->get_id() == 1026)
                        // cout << "Right place in subcondition : " << i << "\n";
                      if( map_entries.find( cs->get_first()->get_id() ) != map_entries.end() &&
                          map_entries[ cs->get_first()->get_id() ].find( delta_prime->get_first()->get_id() ) != map_entries[ cs->get_first()->get_id() ].end() ){
                            // if(delta->get_id() == 1026)
                                // cout << "Wrong place in map_entries : " << i << "\n";
                            continue;
                          }

                        // if(delta->get_id() == 1026)
                            // cout << "Possible concat id after : " << i << "\n";
                      SVal joined_seg = concat(cs, delta_prime, false);
                      // if( map_entries.find( joined_seg->get_first()->get_id() ) != map_entries.end() &&
                      //     map_entries[ joined_seg->get_first()->get_id() ].find( joined_seg->get_second()->get_id() ) != map_entries[ joined_seg->get_first()->get_id() ].end() )
                      //       continue;

            					// if(isPotentialCandidate(delta, joined_seg) && withinAngle(delta, joined_seg)){
                      // if( withinAngle( delta, joined_seg )) {
            						joined.push_back(joined_seg);
                        insertInMapEntries( joined_seg, map_entries );

            						// if(delta->get_id() == 1026)
            			        // cout << "Joined 2 ... " << "\t" << cs->get_id() << "\t" << delta_prime->get_id() << "\t" << joined_seg->get_id() << "\t" << joined_seg->get_first()->get_id() << "\t" << joined_seg->get_second()->get_id() << "\n";
                  }
              }
          }

        else if((isWithinCE(delta->get_second(), cs->get_second(), ce) && !isWithinCE(delta->get_first(), cs->get_first(), ce))
			     ||	(isWithinCE(delta->get_first(), cs->get_second(), ce) && !isWithinCE(delta->get_second(), cs->get_first(), ce))){

             // if(delta->get_id() == 1026)
                 // cout << "Wrong place in self join 1 ... " << cs->get_id() << "\t" << cs->get_first()->get_id() << "\t" << cs->get_second()->get_id() << "\n";
            for(auto& i : candidateMatches){
                boost::shared_ptr<Segment_creator> delta_prime = segments[1].find(i)->second;
                if( same_entries(delta_prime, cs) || in_parent(delta_prime, cs))
                  continue;

                if(cs->get_first() == delta_prime->get_second()){
                  if( map_entries.find( cs->get_second()->get_id() ) != map_entries.end() &&
                      map_entries[ cs->get_second()->get_id() ].find( delta_prime->get_first()->get_id() ) != map_entries[ cs->get_second()->get_id() ].end() )
                        continue;

                    SVal joined_seg = concat(delta_prime, cs, true);
                    // if( map_entries.find( joined_seg->get_first()->get_id() ) != map_entries.end() &&
                    //     map_entries[ joined_seg->get_first()->get_id() ].find( joined_seg->get_second()->get_id() ) != map_entries[ joined_seg->get_first()->get_id() ].end() )
                    //       continue;

                    joined.push_back(joined_seg);
                    insertInMapEntries( joined_seg, map_entries );
					         // if(delta->get_id() == 1026)
			                cout << "Joined ... " << "\t" << joined_seg->get_first()->get_id() << "\t" << joined_seg->get_second()->get_id() << "\n";
                }
				        else if(cs->get_first() == delta_prime->get_first()){

                    if( map_entries.find( cs->get_second()->get_id() ) != map_entries.end() &&
                        map_entries[ cs->get_second()->get_id() ].find( delta_prime->get_second()->get_id() ) != map_entries[ cs->get_second()->get_id() ].end() )
                          continue;

                    SVal joined_seg = concat(delta_prime, cs, false);
                    // if( map_entries.find( joined_seg->get_first()->get_id() ) != map_entries.end() &&
                    //     map_entries[ joined_seg->get_first()->get_id() ].find( joined_seg->get_second()->get_id() ) != map_entries[ joined_seg->get_first()->get_id() ].end() )
                    //       continue;

          					// if(isPotentialCandidate(delta, joined_seg) && withinAngle(delta, joined_seg)){
          					// 	joined.push_back(joined_seg);
          					// }
                    // if(withinAngle(delta, joined_seg)){
          						joined.push_back(joined_seg);
                      insertInMapEntries( joined_seg, map_entries );
          					// }
          					// if(delta->get_id() == 1026)
          			       // cout << "Joined ... " << "\t" << joined_seg->get_first()->get_id() << "\t" << joined_seg->get_second()->get_id() << "\n";
                }
            }
        }

        else{   //case when none of the terminal node coincides. //Missing some cases.
            // cout << "in self join 1 ...  " << cs->get_id() << "\n";
            for(auto& i : candidateMatches){
                SVal delta_prime = segments[1].find(i)->second;

                if( same_entries( delta_prime, cs ) || in_parent( delta_prime, cs ) )
                    continue;

                // cout << "cm ... " << delta_prime->get_id() << "\t" << delta_prime->get_first()->get_id() << "\t" << delta_prime->get_second()->get_id() << "\n";

                SVal joined_seg = joinConcatUtil( cs, delta_prime, map_entries );
                if( joined_seg ){
                  // if( map_entries.find( joined_seg->get_first()->get_id() ) != map_entries.end() &&
                  //     map_entries[ joined_seg->get_first()->get_id() ].find( joined_seg->get_second()->get_id() ) != map_entries[ joined_seg->get_first()->get_id() ].end() )
                  //       continue;

                    joined.push_back( joined_seg );
                    insertInMapEntries( joined_seg, map_entries );
                    // cout << "joined ... " << joined_seg->get_id() << "\t" << joined_seg->get_first()->get_id() << "\t" << joined_seg->get_second()->get_id() << "\n";
                }
            }
        }

        return joined;
    } // old_version.

    void insertInMapEntries( SVal r, map<SKey, map<SKey, vector<SVal>>> &map_entries ){
        if( map_entries.find( r->get_first()->get_id() ) != map_entries.end() ){
            if( map_entries[ r->get_first()->get_id() ].find( r->get_second()->get_id() ) != map_entries[ r->get_first()->get_id() ].end() )
                map_entries[ r->get_first()->get_id() ][ r->get_second()->get_id() ].push_back( r );
            else
                map_entries[ r->get_first()->get_id() ][ r->get_second()->get_id() ] = { r };
        }
        else {
            map_entries[ r->get_first()->get_id() ][ r->get_second()->get_id() ] = { r };
        }
        if( map_entries.find( r->get_second()->get_id() ) != map_entries.end() ){
            if( map_entries[ r->get_second()->get_id() ].find( r->get_first()->get_id() ) != map_entries[ r->get_second()->get_id() ].end() )
                map_entries[ r->get_second()->get_id() ][ r->get_first()->get_id() ].push_back( r );
            else
                map_entries[ r->get_second()->get_id() ][ r->get_first()->get_id() ] = { r };
        }
        else {
            map_entries[ r->get_second()->get_id() ][ r->get_first()->get_id() ] = { r };
        }
    }

    bool in_parent( SVal delta, SVal cs ){
        for( auto it : cs->parents() ){
            if( it == delta->get_id() )
                return true;
        }
        return false;
    }


    SVal concat(SVal base, SVal other, bool same_order){
        osmium::object_id_type id = modified_segId;
    		double length = concatLength(base, other);
    		osmium::object_id_type way_id = base->get_way_id();
    		const osmium::TagList& tags = base->tags();
    		boost::shared_ptr<Rectangle> box = concatMBR(base->get_box(), other->get_box());
    		--modified_segId;

    		boost::shared_ptr<Terminal_node> first, second;
    		vector<IntNode> cl;
    		Geometry* buff;

        // cout << "in concat ... \n";
        // cout << "base ... " << base->get_first()->get_id() << "\t" << base->get_second()->get_id() << "\n";
        // cout << "other ... " << other->get_first()->get_id() << "\t" << other->get_second()->get_id() << "\n";

    		if(same_order){
      			first = base->get_first();
    		    second = other->get_second();
    		    cl = concatIntermediateNodes(base->get_intermediate_nodes(), other->get_intermediate_nodes());
    		    buff = concatBuff(cl);
            // cout << "inside concat 1 ... \n";
            // cout << "base ... " << base->get_first()->get_id() << "\t" << base->get_second()->get_id() << "\n";
            // cout << "other ... " << other->get_first()->get_id() << "\t" << other->get_second()->get_id() << "\n";
    		}
    		else{
      			if(base->get_first()->get_id() == other->get_first()->get_id()){
        				first = base->get_second();
        				second = other->get_second();
        				vector<IntNode> cl_reverse = base->get_intermediate_nodes();
        				reverse(cl_reverse.begin(), cl_reverse.end());

        				cl = concatIntermediateNodes(cl_reverse, other->get_intermediate_nodes());
        				buff = concatBuff(cl);
                // cout << "inside concat 2 ... \n";
      			}
      			else if(base->get_second()->get_id() == other->get_second()->get_id()){
        				first = base->get_first();
        				second = other->get_first();
        				vector<IntNode> cl_reverse = other->get_intermediate_nodes();
        				reverse(cl_reverse.begin(), cl_reverse.end());
        				cl = concatIntermediateNodes(base->get_intermediate_nodes(), cl_reverse);
        				buff = concatBuff(cl);
                // cout << "inside concat 3 ... \n";
      			}
            // else{
            //     cout << base->get_
            //     cout << "none true in concat ... \n";
            // }

    		}
        assert( first );
        assert( second );
        boost::shared_ptr<Segment_creator> joined(new Segment_creator(id, first, second, length, way_id, box, buff, cl, ce, tags));
    		joined->insert_parent(base->get_id());
    		joined->insert_parent(other->get_id());
    		if(base->get_id() < 0)
    			joined->insert_parent(base->parents());
    		if(other->get_id() < 0)
    			joined->insert_parent(other->parents());
        return joined;
    }

    void draw_segment( SVal seg ){
      string location("./verify_data/");
      location += input_file;
      location += "_";
      // cout << "location ... " << location << "\n\n\n";
      ofstream test_cs(location +to_string(seg->get_id())+".txt");
      test_cs << "LINESTRING (" << std::fixed;
      CoordinateArraySequence* cl = new CoordinateArraySequence();
      cl = formCA(seg->get_intermediate_nodes());
      size_t count=0;
      size_t size = cl->getSize();
      while(count < size){
        Coordinate c = cl->getAt(count);
        if(count == 0)
          test_cs << c.x << " " << c.y;
        else
          test_cs << ", " << c.x << " " << c.y;
        count++;
      }
      test_cs << ")";
      // cout << "\n\n";
    }


    vector<SVal> interpolate(tn source, SVal delta, SVal cs ){
        //find all intermediate nodes within circular error of the node from the query segment to be interpolated.

        vector< SVal > result;
        vector<SVal> segs;
        SVal clipper;
        size_t index;

        SKey parent_cs = parentSegmentToBeInterpolated( source, cs );
        if( parent_cs < 0 )
            return result;

        if( modified_nodes.find(source->get_id()) != modified_nodes.end()){
            // if( delta->get_id() == 1018 && parent_cs == 1164 )
            //     cout << "expected place in interpolate() ... \n";
            if( modified_nodes[ source->get_id() ].find( parent_cs ) != modified_nodes[ source->get_id() ].end() ){
                // if( delta->get_id() == 1018 && parent_cs == 1164 )
                //   cout << "Inside modified nodes check for ... " << delta->get_id() << "\t" << cs->get_id() << "\n";
                vector< SVal > result;
                tn node1 = modified_nodes[ source->get_id() ][ parent_cs ];
                segs = interpolated_segments[ node1->get_id() ];
                // if( delta->get_id() == 1018 && parent_cs == 1164 )
                    // cout << "entry in modified nodes ... " << source->get_id() << "\t" << parent_cs << "\t" << node1->get_id() << "\t" << segs.size() << "\n";
                tn node2;
                index = getClipper( delta, segs );
                clipper = segs[ index ];

                assert( clipper->get_id() != 0 );
            }
        }

        SVal actual_cs = segments[ 1 ][ parent_cs ];

        if( !clipper ){
            if( delta->get_id() == 1018 && parent_cs == 1164)
                cout << "not previously interpolated ... \n";
            map<size_t, Coordinates> nearest; // vector that stores the nearest two nodes of the data segment, to serve for the interpolation.
            map<size_t, vector<Coordinates>> pair;

            vector<IntNode> cl = actual_cs->get_intermediate_nodes();
            size_t pos=0;

            pos = getNearestNodesBuffer( source, cl );
            cout << "\n\npos in interpolate ... " << pos << "\n\n";
            if( pos < cl.size()-1 ){
                vector< Coordinates > coords;
                coords.push_back( cl[ pos ].get_location() );
                coords.push_back( cl[ pos+1 ].get_location() );
                pair.insert( { pos, coords } );
            }
            else
                return segs;

            if(pair.size() > 0){
                segs = interpolatePairs(delta, source, pair, actual_cs, parent_cs ); //two consectuive nodes are within circular error.
                // return segs;
            }
            if( segs.size() > 0 ){
                index = getClipper( delta, segs );
                clipper = segs[ index ];
            }

        }

        if( !clipper )
            return result;

        cout << "delta ... " << delta->get_id() << "\t" << delta->get_first()->get_id() << "\t" << delta->get_second()->get_id() << "\n";
        cout << "cs ... " << cs->get_id() << "\t" << cs->get_first()->get_id() << "\t" << cs->get_second()->get_id() << "\n";
        cout << "clipper ... " << clipper->get_id() << "\t" << clipper->get_first()->get_id() << "\t" << clipper->get_second()->get_id() << "\n";
        tn inode = modified_nodes[ source->get_id() ][ parent_cs ];
        // clipper = interpolated_segments[ inode->get_id() ][ 1 - flag ];
        tn first, second;

        vector< IntNode > cl = cs->get_intermediate_nodes();
        vector< IntNode > cll = clipper->get_intermediate_nodes();

        vector< IntNode > seq = {};

        if( clipper->get_first()->get_id() == cs->get_first()->get_id() ){
            cout << "interpolate 1 ... \n";
            vector< IntNode > clipped( cl.begin() + cll.size()-1, cl.end() );
            clipped.insert( clipped.begin(), cll.back() );
            first = inode;
            second = cs->get_second();
            seq = clipped;
        }
        else if( clipper->get_second()->get_id() == cs->get_second()->get_id() ){
            cout << "interpolate 2 ... \n";
            vector< IntNode > clipped( cl.begin(), cl.end() - cll.size() + 2 );
            clipped.push_back( cll.front() );
            first = cs->get_first();
            second = inode;
            seq = clipped;
        }
        else if( clipper->get_first()->get_id() == cs->get_second()->get_id() ){
            cout << "interpolate 3 ... \n";
            vector< IntNode > clipped( cl.begin(), cl.end() - cll.size() + 2 );
            clipped.push_back( cll.back() );
            first = cs->get_first();
            second = inode;
            seq = clipped;
        }
        else if( clipper->get_second()->get_id() == cs->get_first()->get_id() ){
            cout << "interpolate 4 ... \n";
            if( cl.begin() + cll.size() - 1 >= cl.end() )
                return result;
            vector< IntNode > clipped( cl.begin() + cll.size()-1, cl.end() );
            clipped.insert( clipped.begin(), cll.front() );
            first = inode;
            second = cs->get_second();
            seq = clipped;
        }
        // assert( seq.size() > 1 );
        if( seq.size() <= 1 ){
            cout << "seq size is less than 1!! \n";
            return result;
        }

        cout << "first second ... " << first->get_id() << "\t" << second->get_id() << "\n";
        cout << "start and abck of seq ... " << seq.front().get_id() << "\t" << seq.back().get_id() << "\n";
        // draw_segment( clipper );

        osmium::object_id_type id = modified_segId--;
        double length = calculateLength( seq );
        osmium::object_id_type way_id = segments[ 1 ][ parent_cs ]->get_way_id();
        boost::shared_ptr<Rectangle> box = calculateBox( seq );
        Geometry* buff = calculateBuffer( seq );
        uint32_t CE = segments[ 1 ][ parent_cs ]->getCircularError();
        const osmium::TagList& tags = segments[ 1 ][ parent_cs ]->tags();
        SVal iseg( new Segment_creator(id, first, second, length, way_id, box, buff, seq, CE, tags) );
        iseg->insert_parent( parent_cs );
        iseg->insert_parent( segs[ 1-index ]->get_id() );
        iseg->insert_parent( cs->parents() );

        cout << "id of iseg ... " << iseg->get_id() << "\n";

        // draw_segment( iseg );
        // segs.clear();
        result.push_back( iseg );
        return result;
    }


    size_t getClipper( SVal delta, vector< SVal > segs ){
        double min_overlap = DBL_MAX;
        size_t index, i=0;
        auto buffer = delta->get_buffer();
        for( auto it : segs){
            double o = buffer->intersection( it->get_buffer() )->getArea();
            cout << "in getClipper() ... seg ... " << it->get_first()->get_id() << "\t" << it->get_second()->get_id() << "\n";
            cout << "in getClipper() ... o and min_overlap ... " << o << "\t" << min_overlap << "\n";
            if( o < min_overlap ){
                min_overlap = o;
                index = i;
            }
            i++;
        }
        return index;
    }


    size_t getNearestNodesBuffer( tn source_node, vector< IntNode > cl ){
        Coordinates c_mcp = lonlat_to_mercator( source_node->get_location() );
        Geometry* point_buffer = global_factory->createPoint(Coordinate(pm->makePrecise(c_mcp.x), pm->makePrecise(c_mcp.y)))->buffer( ce );

        size_t pos=cl.size()+1, i=1, size = cl.size();
        double overlap = 0;

        CoordinateArraySequence* clt = new CoordinateArraySequence();
        c_mcp = lonlat_to_mercator( cl[ 0 ].get_location() );
        auto p1 = Coordinate(pm->makePrecise(c_mcp.x), pm->makePrecise(c_mcp.y));

        while( i < size ){
            c_mcp = lonlat_to_mercator( cl[ i ].get_location() );
            auto p2 = Coordinate(pm->makePrecise(c_mcp.x), pm->makePrecise(c_mcp.y));
            clt->add( p1 );
            clt->add( p2 );
            // cout << "clt size ... " << clt->size() << "\n";
            const Geometry* ls = global_factory->createLineString(clt);
            // cout << "length of linestring ... " << ls->getLength() << "\n";
            Geometry* line_buffer = ls->buffer( ce );
            // cout << "area of linestring ... " << line_buffer->getArea() << "\n";
            double o = line_buffer->intersection( point_buffer )->getArea();
            // cout << "area in getNearestNodesBuffer ... " << o << "\n";
            if( o > overlap ){
                overlap = o;
                pos = i;
            }
            p1 = p2;
            clt->clear();
            i++;
        }
        cout << "\n\n";
        return pos-1;
    }


    vector<SVal> interpolatePairs(SVal delta, boost::shared_ptr<Terminal_node> source, map<size_t, vector<Coordinates>> pair, SVal cs, SKey parent_cs ){

        vector<SVal> result;
        assert( pair.size() == 1 );
        auto node_seg_pair = modified_nodes.find(source->get_id());
        if(node_seg_pair != modified_nodes.end()){
            auto item = node_seg_pair->second.find(parent_cs);
            if(item != node_seg_pair->second.end()){
                tn node1 = item->second;
                tn node2;
                if( flag == 1 )
                    node2 = cs->get_second();
                else
                    node2 = cs->get_first();
                if( end_points_segments.find( node1->get_id() ) != end_points_segments.end() ){
                    if( end_points_segments[ node1->get_id() ].find( node2->get_id() ) != end_points_segments[ node1->get_id() ].end() ){
                        result.push_back( end_points_segments[ node1->get_id() ][ node2->get_id() ] );
                        return result;
                    }
                }
                if( end_points_segments.find( node2->get_id() ) != end_points_segments.end() ){
                    if( end_points_segments[ node2->get_id() ].find( node1->get_id() ) != end_points_segments[ node2->get_id() ].end() ){
                        result.push_back( end_points_segments[ node2->get_id() ][ node1->get_id() ] );
                        return result;
                    }
                }
            }
        }
        // if( delta->get_id() == 1018 )
        //     cout << "interpolating for 1018 ... using " << cs->get_id() << " ... \n";
        Coordinates c_query = source->get_location();
        vector<boost::shared_ptr<Terminal_node>> allTNodes;
        for(auto& item : pair){
            vector<Coordinates> coord = item.second;
            Coordinates c = getInterpolatedNodeCoordinates(c_query, coord[0], coord[1]);
            vector<IntNode> ci = cs->get_intermediate_nodes();
            vector<vector<IntNode>> split_segments = insert(source, c, ci, item.first, 1); //to insert c at the right position in ci, it will check if c lies on the segment or on the extension of the segment. It returns the truncated CoordinateArraySequence.
            if( split_segments[ 0 ].size() < 2 || split_segments[ 1 ].size() < 2 )
                continue;
            boost::shared_ptr<Terminal_node> tnode(new Terminal_node(modified_nodeId, c, 1));

            allTNodes.push_back(tnode);
            --modified_nodeId;
            vector<SVal> segs = getInterpolatedSegment(source, cs, split_segments, tnode);
            int flag = 1 - getClipper( delta, segs );
            if( segs[ flag ] ){

                // cout << "add in modified nodes ... " << source->get_id() << "\t" << parent_cs << "\t" << tnode->get_id() << "\n";

                modified_nodes[ source->get_id() ][ parent_cs ] = tnode;
                end_points_segments[ segs[ flag ]->get_first()->get_id() ][ segs[ flag ]->get_second()->get_id() ] = segs[ flag ];
                end_points_segments[ segs[ 1-flag ]->get_first()->get_id() ][ segs[ 1-flag ]->get_second()->get_id() ] = segs[ 1-flag ];
                // mnodes << tnode->get_id() << "\t" << c.x << "\t" << c.y << "\t" << delta->get_id() << "\t" << source->get_id() << "\t" << cs->get_id() <<  "\n";
                interpolated_segments.insert({tnode->get_id(), segs});

                parent_interpolated_segments[ tnode->get_id() ] = parent_cs;
                return segs;
            }
        }
        return result;
    }


    vector<SVal> getInterpolatedSegment(boost::shared_ptr<Terminal_node> source, SVal source_segment, vector<vector<IntNode>> split_segments, boost::shared_ptr<Terminal_node> other){
        vector<SVal> int_segs;
        SVal seg1, seg2;
        boost::shared_ptr<Terminal_node> first, second, third;
        first = source_segment->get_first();
        second = other;
        third = source_segment->get_second();
        if(split_segments[0].size() > 1 && split_segments[1].size() > 1){
            osmium::object_id_type id1 = modified_segId--;
            osmium::object_id_type id2 = modified_segId--;
            double length1 = calculateLength(split_segments[0]);
            double length2 = calculateLength(split_segments[1]);
            osmium::object_id_type way_id = source_segment->get_way_id();
            boost::shared_ptr<Rectangle> box1 = calculateBox(split_segments[0]);
            boost::shared_ptr<Rectangle> box2 = calculateBox(split_segments[1]);
            Geometry* buff1 = calculateBuffer(split_segments[0]);
            Geometry* buff2 = calculateBuffer(split_segments[1]);
            uint32_t CE = source_segment->getCircularError();
			      const osmium::TagList& tags = source_segment->tags();
            SVal ptr1(new Segment_creator(id1, first, second, length1, way_id, box1, buff1, split_segments[0], CE, tags));
            SVal ptr2(new Segment_creator(id2, second, third, length2, way_id, box2, buff2, split_segments[1], CE, tags));
            seg1 = ptr1;
            seg2 = ptr2;
				seg1->insert_parent(source_segment->get_id());
				seg2->insert_parent(source_segment->get_id());
				seg1->insert_parent(source_segment->parents());
				seg2->insert_parent(source_segment->parents());
        }
        int_segs.push_back(seg1);
        int_segs.push_back(seg2);
        return int_segs;
    }


    inline double concatLength(boost::shared_ptr<Segment_creator> seg1, boost::shared_ptr<Segment_creator> seg2){
        return seg1->get_length() + seg2->get_length();
    }

    //insert Coordinate c in CoordinateArraySequence cl at position pos.
    vector<vector<IntNode>> insert(tn source, Coordinates c, vector<IntNode> &cl, size_t pos, bool pair){ //truncate the CoordinateArraySequence after the interpolated node. The interpolated node will be the terminal node.
        size_t final_pos=0;

        if(pair){
            Coordinates c1 = cl[pos].get_location();
            Coordinates c2 = cl[pos+1].get_location();
            if(liesOnSegment(c, c1, c2)){
                cl.insert(cl.begin()+pos+1, IntNode{modified_nodeId, c});
                final_pos = pos+1;
            }
            else{
                if(liesToLeft(c, c1, c2)){
                    cl.insert(cl.begin()+pos, IntNode{modified_nodeId, c});
                    final_pos = pos;
                }
                else{
                    if(pos+1 == cl.size())
                        cl.push_back(IntNode{modified_nodeId, c});
                    else
                        cl.insert(cl.begin()+pos+1, IntNode{modified_nodeId, c});
                    final_pos = pos+1;
                }
            }
        }
        else{ //Right now assuming the interpolated node will fall on the edge. Correct this maybe later.
            cl.insert(cl.begin()+pos, IntNode{modified_nodeId, c});
            final_pos = pos;
        }
        if( source->get_id() == 4062883790 )
            cout << "final pos ... " << final_pos << "\n";
        size_t total = cl.size();
		// if((delta->get_id() == 1189 || delta->get_id() == 1188) && source->get_id() == 7050936898)
		// 	cout << "Total size after interpolation ! " << delta->get_id() << "\t" << cl.size() << "\n";
        // vector<Coordinates> split_segment;
        vector<vector<IntNode>> split_segments;
        vector<IntNode> split_segment1;
        vector<IntNode> split_segment2;

        size_t i=0;
        //flag is to determine whether to keep the first half or the second half of the segment.
        while(i<=final_pos){
            split_segment1.push_back(cl[i]);
            ++i;
        }
        split_segment2.push_back(cl[i-1]); //so that the split segments are continuous and not broken.
        while(i < total){
            split_segment2.push_back(cl[i]);
            ++i;
        }
        split_segments.push_back(split_segment1);
        split_segments.push_back(split_segment2);
        return split_segments;
    }
};


#endif
