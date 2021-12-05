#ifndef IRANKJOINCA_HPP
#define IRANKJOINCA_HPP

#include "SegmentHandler.hpp"
#include "Util.hpp"
// #include "FindMatches_angle.hpp"
#include "FindMatches_wprints.hpp"
#include "RJUtilReGroup_ablation.hpp"

#include <map>
#include <vector>
#include <set>
#include <fstream>

ofstream score_ana("score_analysis.txt");
ofstream avg_gap("avg_gap.txt");

set< SKey > not_matched;
set< SKey > to_verify_matches;
set<SKey> single_matches;

class RankJoinCA{
private:

    ofstream outfile;
    string str_o = "output_ca";
    ofstream filedebug;
    string str_d = "debug_ca";
    ofstream analysis;
    string str_a = "ana_ca";
    ofstream jo;

public:

    size_t covers_cnt = 0;
    size_t endpoints_cnt = 0;
    size_t not_intersects_cnt = 0;
    size_t total_cnt = 0;
    size_t match_dummy_cnt = 0;

    size_t covers_cnt2 = 0;
    size_t endpoints_cnt2 = 0;
    size_t not_intersects_cnt2 = 0;

    static bool sort_output(JoinResult result1, JoinResult result2){
        return result1.score > result2.score;
    }

    void processRankJoin(map<SKey, RVal> relations, char* ck, char* count, bool la, vector<JoinResult> &join_result){
        // exit( EXIT_FAILURE );
        RankJoinUtil rj_utils;
        k = static_cast<size_t>(stoul(ck));

        cout << "k ... " << *ck << "\t" << k << "\n";


        rj_utils.findSingleMatches(single_matches, relations);
        joins                       = findJoins(segments[0], single_matches);
        cout << "Size of joins : " << joins.size() << "\n\n";
        cout << "total relations ... " << relations.size() << "\n\n";
        cout << "segments[0] length ... " << segments[ 0 ].size() << "\n\n";
        for(auto it1 : joins){
            adjlist1 << it1.first << " \n";
            for(auto& it2 : it1.second)
                adjlist1 << "\t" << it2 << "\n";
        }
        // return output;

        if(*count != 0){
            size_t                          join_size               = static_cast<size_t>(stoul(count));
                                                                      // rj_utils.formQueryGroupsUsingBFS( joins, join_size );
                                                                      // rj_utils.formFixedSizeBFSGroups( joins, 5 );
                                                                      rj_utils.unitSizeGroups( joins );
        }
        else{
            // exit( EXIT_FAILURE );
            size_t join_size = relations.size();
            rj_utils.formFixedSizeBFSGroups( joins, 10 );
        }

        // n                       = join_order_group.size();             //number of relations.
        n = relation_groups.size();
        cout << "Join order group size : " << n << "\n\n";

        size_t n_prime = rj_utils.verifySize( );
        cout << "Join size from group order join : " << n_prime << "\n";

        str_o.append(count).append(".txt");
        str_d.append(count).append(".txt");
        str_a.append(count).append(".txt");
        outfile.open(str_o);
        filedebug.open(str_d);
        analysis.open(str_a);


        analysis << "Number of relations groups = " << n << "\n\n";

        intermediate_results          = rj_utils.initializeVector( n );
        intermediate_results_groups   = rj_utils.initializeVectorGroups( n-1 );
        size_of_relations             = rj_utils.getSize( relations );           //size of each relation.
        start_indices                 = rj_utils.initializeStartIndices( n-1 ) ;
                                        rj_utils.calculateLengths( segments[0], la );
                                        rj_utils.initializeInputBuffersGroup( n );
                                        rj_utils.calculateInitialScoreBounds();
                                        rj_utils.accumulate_dp_max_scores( dp_max_scores );
        min_output_score              = -DBL_MAX;

        cout << "Intermediate result size : " << intermediate_results.size() << "\n";
        cout << "Intermediate results group size : " << intermediate_results_groups.size() << "\n";
        cout << "t relations bounds : " << t_relation_bounds.size() << "\n";
        cout << "t bounds size : " << t_bounds.size() << "\n";
        cout << "Size of score bounds ... " << tau_score_bounds.size() << "\n";
        cout << "Tau ... " << tau << "\n";
        for( auto it : join_order_group )
            cout << it << "\t" << tau_score_bounds[ it ] << "\n";
        // return output;
        size_t ctr = 0, i = 0;
        bool reGrouped = false;
        while( output.size() < k || int(tau*1000000) > int(min_output_score*1000000) ){

            cout << "Pull Bound starting ... ctr : " << ctr++ << "\n";
            cout << "Pass : " << ctr << "\n";

            // if( !reGrouped )
            i = rj_utils.chooseInput( join_order_group, tau_score_bounds );

            if( output.size() < 1 && i > 1 && input_buffers_groups[ join_order_group[ i ] ].size() > 0 ){
                ir_tuple ir = input_buffers_groups[ join_order_group[ i ] ][ 0 ];
                auto prev_irs = intermediate_results_groups[ i-2 ];
                cout << "Size of intermediate results in each group ... " << intermediate_results[ i ].size() << "\n";
                cout << "prev irs size ..... " << prev_irs.size() << "\n";
                map< size_t, size_t > group_failing;
                for( auto prev_ir = prev_irs.begin(); prev_ir != prev_irs.end(); prev_ir++ ){
                    double penalty = 0;
                    if( rj_utils.joinsWith( ( *prev_ir ).join_result, ir.join_result, i, penalty, group_failing ) ){
                        cout << "joinsWith() succeeded ... \n";
                        map< SKey, TVal > new_ir = ( *prev_ir ).join_result;
                        new_ir.insert( ir.join_result.begin(), ir.join_result.end() );
                        ir_tuple res{ new_ir, ( *prev_ir ).score + ir.score, ( *prev_ir ).penalty + ir.penalty + penalty };
                        intermediate_results_groups[ i-1 ].push_back( res);
                        cout << "\nscore of int result ... " << ( res.score + dp_max_scores.back() - dp_max_scores[ i ] )/rj_utils.log_total_length << "\n\n";
                        // assert( ( res.score + dp_max_scores.back() - dp_max_scores[ i ] )/rj_utils.log_total_length >= tau );
                    }
                }
                assert( intermediate_results_groups[ i-1 ].size() == 0 );

                size_t max_fail_group = -1, max_fail = 0;
                for( auto it : group_failing ){
                    if( it.second > max_fail ){
                        max_fail = it.second;
                        max_fail_group = it.first;
                    }
                }
                assert( max_fail_group != -1 );
                cout << "max fail group ... " << max_fail_group << "\n";
                rj_utils.reGroup( join_order_group[ i ], max_fail_group, i );
                cout << "i after group ... " << i << "\n";
                cout << "size of int results after regroup ... " << intermediate_results_groups[ i-1 ].size() << "\n";
                reGrouped = true;

            }

            reGrouped = false;
            size_t group_number = join_order_group[ i ];
            cout << "Join order index picked ... " << i << "\n";
            cout << "Relation group picked ... " << group_number << "\t";

            cout << "Relations in the group:\n";
            for(auto r : relation_groups[ group_number ] )
                cout << r << "\n";

            cout << i << "\n";
            if(i == -1)
                break;

            ir_tuple ir = rj_utils.getNextJoinTuple( i, la );

            if( ir.score == -1 ){
              if( quit_flag )
                  return;
                assert( tau_score_bounds[ group_number ] == -DBL_MAX );
                continue;
            }

            seen_tuples_groups[ group_number ] += 1;

            input_buffers_groups[ group_number ].push_back( ir );

            cout << "Score of result obtained ... " << ir.score << "\n";

            // cout << "Previous score ... " << prev_scores[ group_number ] << "\n";

            cout << "tau before ... " << tau_score_bounds[ group_number ] << "\n";

            if( seen_tuples_groups[ group_number ] == 1 )
                rj_utils.updateAllTaus( i );

            double la_score = 0;
            if( la ){
              double la_score = t_bounds[ group_number ];
              if( intermediate_results[ i ].size() > seen_tuples_groups[ group_number ] ){
                  ir_tuple la_tuple = intermediate_results[ i ][ seen_tuples_groups[ group_number ] ];
                  double score = ( la_tuple.score - la_tuple.penalty );  //there no penalty associated with lookahead tuple yet.
                  if( score > la_score )
                     la_score = score;
              }
            }
            else{
                la_score = ir.score - ir.penalty;
            }

            rj_utils.updateGroupBound( i, la_score );

            cout << "Look ahead score ... " << la_score << "\n";

            cout << "New tau bound for group ... " << tau_score_bounds[ group_number ] << "\n";
            cout << "New tau bound ... " << tau << "\n";

            if( i == 0 ){
                size_t start_index = intermediate_results_groups[ 0 ].size();
                start_indices[ 0 ] = start_index;
                cout << "Size of intermediate results in each group ... " << intermediate_results[ i ].size() << "\n";
                if( n > 1 ){
                    if( regrouped_groups.find( join_order_group[ 1 ] ) != regrouped_groups.end() ){
                      intermediate_results_groups[ 0 ].push_back( ir );
                      cout << "\nscore of int result ... " << ( ir.score + dp_max_scores.back() - dp_max_scores[ i ] )/rj_utils.log_total_length << "\n\n";
                      // assert( ( ir.score + dp_max_scores.back() - dp_max_scores[ i ] )/rj_utils.log_total_length >= tau );
                    }
                    else{
                      for( auto prev_ir : input_buffers_groups[ join_order_group[ 1 ] ] ){
                          double penalty = 0;
                          if( rj_utils.joinsWith( prev_ir.join_result, ir.join_result, i, penalty ) ){ //joins() will also check if there are no joins found in the previous ir to simply take the cross product.
                              // cout << "inside i = 0 ... \n";
                              map< SKey, TVal > new_ir = prev_ir.join_result;
                              new_ir.insert( ir.join_result.begin(), ir.join_result.end() );
                              ir_tuple res{ new_ir, prev_ir.score + ir.score, prev_ir.penalty + ir.penalty + penalty };
                              intermediate_results_groups[ 0 ].push_back( res );
                              cout << "\nscore of int result ... " << ( res.score + dp_max_scores.back() - dp_max_scores[ i ] )/rj_utils.log_total_length << "\n\n";
                              // assert( ( res.score + dp_max_scores.back() - dp_max_scores[ i ] )/rj_utils.log_total_length >= tau );
                          }
                      }
                    }
                }

                rj_utils.removeDuplicates( intermediate_results_groups[ i ] );
                rj_utils.completeJoin( i+2 );
                // rj_utils.formOutputs( output );
                // if( output.size() >= k )
                //     min_output_score = output[ k-1 ].score;
            }

            else if( i == 1 ){
                size_t start_index = intermediate_results_groups[ i-1 ].size();
                start_indices[ i-1 ] = start_index;
                cout << "Size of intermediate results in each group ... " << intermediate_results[ i ].size() << "\n";
                for( auto prev_ir : input_buffers_groups[ join_order_group[ 0 ] ] ){
                    double penalty = 0;
                    if( rj_utils.joinsWith( prev_ir.join_result, ir.join_result, i, penalty ) ){ //joins() will also check if there are no joins found in the previous ir to simply take the cross product.
                        map< SKey, TVal > new_ir = prev_ir.join_result;
                        new_ir.insert( ir.join_result.begin(), ir.join_result.end() );
                        ir_tuple res{ new_ir, prev_ir.score + ir.score, prev_ir.penalty + ir.penalty + penalty };
                        intermediate_results_groups[ i-1 ].push_back( res );
                        cout << "\nscore of int result ... " << ( res.score + dp_max_scores.back() - dp_max_scores[ i ] )/rj_utils.log_total_length << "\n\n";
                        // assert( ( res.score + dp_max_scores.back() - dp_max_scores[ i ] )/rj_utils.log_total_length >= tau );
                    }
                }
                if( intermediate_results_groups[ 0 ].size() == 0 ){
                    cout << "re-grouping ... \n";
                    rj_utils.reGroup( group_number, join_order_group[ 0 ], i );
                    reGrouped = true;
                    continue;
                }
                rj_utils.removeDuplicates( intermediate_results_groups[ i-1 ] );
                rj_utils.completeJoin( i+1 );
            }

            else{
                size_t start_index = start_indices[ i-2 ];
                start_indices[ i-1 ] = intermediate_results_groups[ i-1 ].size();

                auto prev_irs = intermediate_results_groups[ i-2 ];
                cout << "Size of intermediate results in each group ... " << intermediate_results[ i ].size() << "\n";
                cout << "prev irs size ..... " << prev_irs.size() << "\n";
                map< size_t, size_t > group_failing;
                for( auto prev_ir = prev_irs.begin(); prev_ir != prev_irs.end(); prev_ir++ ){
                    double penalty = 0;
                    if( rj_utils.joinsWith( ( *prev_ir ).join_result, ir.join_result, i, penalty, group_failing ) ){
                        cout << "joinsWith() succeeded ... \n";
                        map< SKey, TVal > new_ir = ( *prev_ir ).join_result;
                        new_ir.insert( ir.join_result.begin(), ir.join_result.end() );
                        ir_tuple res{ new_ir, ( *prev_ir ).score + ir.score, ( *prev_ir ).penalty + ir.penalty + penalty };
                        intermediate_results_groups[ i-1 ].push_back( res);
                        cout << "\nscore of int result ... " << ( res.score + dp_max_scores.back() - dp_max_scores[ i ] )/rj_utils.log_total_length << "\n\n";
                        // assert( ( res.score + dp_max_scores.back() - dp_max_scores[ i ] )/rj_utils.log_total_length >= tau );
                    }
                }
                if( intermediate_results_groups[ i-1 ].size() == 0 ){
                    size_t max_fail_group = -1, max_fail = 0;
                    for( auto it : group_failing ){
                        if( it.second > max_fail ){
                            max_fail = it.second;
                            max_fail_group = it.first;
                        }
                    }
                    assert( max_fail_group != -1 );
                    cout << "max fail group ... " << max_fail_group << "\n";
                    rj_utils.reGroup( group_number, max_fail_group, i );
                    cout << "i after group ... " << i << "\n";
                    cout << "size of int results after regroup ... " << intermediate_results_groups[ i-1 ].size() << "\n";
                    reGrouped = true;
                }
                else{
                  rj_utils.removeDuplicates( intermediate_results_groups[ i-1 ] );
                  rj_utils.completeJoin( i+1 );
                  cout << "Intermediate results at " << i-1 << " size ... " << intermediate_results_groups[ i-1 ].size() << "\n";
                }
                if( quit_flag )
                    return;


                size_t j;
                if( i > 0 )
                    j = i-1;
                else
                    j = 0;
                while( j < intermediate_results_groups.size() ){
                  cout << "regroup intermediate_results size at ... " << j << " ... " << intermediate_results_groups[ j ].size() << "\n";
                  j++;
                }
            }
            rj_utils.formOutputs( output );
            if( output.size() >= k )
                min_output_score = output[ k-1 ].score;

            cout << "The value of tau : " << tau << "\n";
            cout << "Minimum output score : " << min_output_score << "\n";
            cout << "Result size : " << intermediate_results_groups.back().size() << "\t" << output.size() << "\t" << k << "\n\n";
            // if( output.size() >= k )
            //     break;
        }
        cout << "number of single matches ... " << single_matches.size() << "\n";
        cout << "printing dummy ... \n";
        cout << dummy << "\t" << dummy->segment()->get_id() << "\n";
        JoinResult res = output[0];
        for( auto it : single_matches ){
            output[ 0 ].result.insert( { it, dummy } );
        }
        outfile << "Score : " << res.score << "\n\n";
        for(auto it : res.result){
            SKey id = it.first;
            outfile << it.first << " " << it.second->segment()->get_id() << "\t" << it.second->is_contained() << "\t"
                    << it.second->segment()->get_first()->get_id() << "\t" << it.second->segment()->get_second()->get_id() << "\n";
            if( relations.find( id ) == relations.end() )
                continue;
            if(relations[ id ]->get_entries()[0]->segment()->get_id() != it.second->segment()->get_id())
                cout << id << "\t" << it.second->segment()->get_id() << "\n";
        }
        cout << "size of result ... " << res.result.size() << "\n\n";
        join_result = output;
        // return output;
    }

    void calculateMatchingAccuracy( vector< JoinResult > res ){
        MatchingUtils mu;
        for( auto it : res[ 0 ].result ){
            SKey id = it.first;
            SVal ref = segments[ 0 ][ id ];
            SVal match = it.second->segment();

            if( it.second->segment()->get_id() == 0 ){
                match_dummy_cnt++;
                not_matched.insert( id );
                continue;
            }
            vector< IntNode > source = ref->get_intermediate_nodes(), shifted = match->get_intermediate_nodes();

            CoordinateArraySequence* cl_source = formCA( source );
            CoordinateArraySequence* cl_shifted = formCA( shifted );

            const Geometry* ls_source = global_factory->createLineString(cl_source);
            Geometry* buff_source = (ls_source->buffer(ce));
            const Geometry* ls_shifted = global_factory->createLineString(cl_shifted);
            // Geometry* buff_shifted = (ls_shifted->buffer(ce));
            if( buff_source->covers( ls_shifted ) )
                covers_cnt++;
            else{
                to_verify_matches.insert( id );
                mu.draw_segment( ref );
                mu.draw_segment( match );
            }


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

            total_cnt++;
        }

    }

    void print_scores(map<SKey, double> score_bounds){
        analysis << "Score bounds : " << "\n";
        for(auto it : score_bounds)
            analysis << it.first << "\t:\t" << it.second << "\n";
    }
};

#endif
