#ifndef RANKJOIN_UTIL_H
#define RANKJOIN_UTIL_H

#include "SegmentHandler.hpp"
#include "Util.hpp"
// #include "FindMatches_angle.hpp"
#include "FindMatches_wprints.hpp"

#include <map>
#include <vector>
#include <set>
#include <fstream>
#include <type_traits>
#include <typeinfo>
#include <algorithm>

ofstream checkj("check.txt");
double t = 1;
SKey input_relation = -1;
// vector<SKey> edge_cases = {3721,3720,3719,3736,1343,1282,2039,2022,1330,2066,1324,1114,1333,1636,1323,1583,1843};
//1282 - seg fault.
//check 1843.
//2022 - seg fault.
//3720 - infinite loop, out of removing duplicates.
//3738 - infinite loop.
//3750 - seg fault -- no seg fault now.
//1324 - seg fault.
//1636 - seg fault.
//1327 - seg fault.
// vector<SKey> edge_cases = {3738,2022,1324,1636,1759,3720,3763,1327,1366}; //- formQueryGroupsUsingBFS() - 1.7 seconds


//1330,1843,1583,3721,3750,1323,3720,1324,2039 - works
//1239,1114,3719,1282,2022,1327,1636 problem.
// - formStarGroups() - ~2.3 seconds
// vector<SKey> edge_cases = {1239,1636,1759,3763,1282,1327,1366};  // - formStarGroups()

//1330,1843,1583,3721,3750,1323,3720,1324,2039,1114,3719,2022,3738 - works
//1239,1327,1636,1282 problem.
// - formStarGroups() with singleton removed - ~7.5 seconds (result size 16)

//below - approximate join implemented. It is giving result for all the segments in edge_cases.
vector<SKey> edge_cases = {};
//1114,1330 - works, 3738,2039,2022,1843,1583,3719,3750,3720,3721
//1636,1366,1759,1327,1239,3763,1282

// vector<SKey> edge_cases = {3738,2039,1239,2022,1324,1636,1323,1843,1759,1330,1114,1583,3719,3750,1282,3720,3721,3763,1327,1366};

struct ir_tuple{
    map< SKey, TVal > join_result;
    double score;
    double penalty;
};

struct JoinParameters{
    vector<SKey> join_order;
    map<SKey, SKey> query_graph;
};

struct JoinResult{                                                                              //represents one join result.
    map<SKey, TVal> result;
    double score; //score willbe the unnormalized sum of scores in the join tuple.
};

using TVal = boost::shared_ptr<Tuple>;
using IRs = vector< ir_tuple >;

TVal                            dummy;
map<SKey, vector<SKey>>         joins;
map<SKey, vector<TVal>>         input_buffers;
map< size_t, IRs >              input_buffers_groups;
vector<IRs>                     intermediate_results; //for each group.
vector<JoinResult>              output;                                                 //the value is again a list of fixed size.
map<SKey, size_t>               seen_tuples;                                            //the number of tuples seen for each relation.
set<SKey>                       exhausted_relations;                                    //the relations for which all tuples are seen.
map<SKey, size_t>               size_of_relations;           //size of each relation.
vector< size_t>                 start_indices;
map<SKey, double>               lengths;
double                          min_output_score;
size_t                          n;
double                          tau = 1;
map< size_t, vector<SKey> >     relation_groups;
map< size_t, size_t >           query_graph_group;
vector< size_t >                join_order_group;
map< size_t, double >           tau_score_bounds;   //the score bound of all relations.
map< size_t, double >           t_bounds;           //the bound of score result for the relation group.
map< SKey, double >             t_relation_bounds; //the usual t bound for each relation.
map< size_t, double >           max_scores_groups;
map< size_t, double >           max_scores_lengths;
vector< IRs >                   intermediate_results_groups; //combining the groups in left deep join manner.
map< size_t, size_t >           seen_tuples_groups;
map< size_t, double >           prev_scores;
map< SKey, size_t >             degree;
map< SKey, SKey >               approximate_joins; //map of node ids to calculate the penalty for a pair of nodes only once.
size_t                          k;
vector< double >                dp_max_scores;
// map< size_t, double >           group_score_weights;
map< SKey, size_t >             group_number;
set< SKey >                     regrouped_groups;


struct compareDegree {
    bool operator()(SKey id1, SKey id2 ){
        return degree[ id1 ] > degree[ id2 ];
    }
};

class RankJoinUtil{

    double total_length = 0;
    double total_max_score = 0;
    // double total_group_weights = 0;

public:
    double log_total_length = 0;

    static bool sort_output(JoinResult result1, JoinResult result2){
        return result1.score > result2.score;
    }

    static bool sort_ir( ir_tuple ir1, ir_tuple ir2 ){
        return (ir1.score - ir1.penalty) > (ir2.score - ir2.penalty);
    }

    void reGroup( SKey curr_group, SKey prev_group, size_t &i ){
        regrouped_groups.insert( curr_group );

        double prev_score_sum = ( input_buffers_groups[ prev_group ][ 0 ].score - input_buffers_groups[ prev_group ][ 0 ].penalty )
                              + ( input_buffers_groups[ curr_group ][ 0 ].score - input_buffers_groups[ curr_group ][ 0 ].penalty );


        double prev_group_scores = 0, curr_group_scores = 0;
        for( auto it : relation_groups[ prev_group ] ){
            double length = lengths[ it ];
            prev_group_scores += ( log10( length ) * relations[ it ]->get_entries()[ 0 ]->score() );
        }
        for( auto it : relation_groups[ curr_group ] ){
            double length = lengths[ it ];
            curr_group_scores += ( log10( length ) * relations[ it ]->get_entries()[ 0 ]->score() );
        }

        // recalculate t_relation_bounds
        for( auto it : relation_groups[ prev_group ] )
            t_relation_bounds[ it ] += curr_group_scores;

        for( auto it : relation_groups[ curr_group ] ){
            t_relation_bounds[ it ] += prev_group_scores;
            group_number[ it ] = prev_group;
        }

        vector< SKey > merged_group = relation_groups[ prev_group ];
        merged_group.insert( merged_group.end(), relation_groups[ curr_group ].begin(), relation_groups[ curr_group ].end() );

        relation_groups[ prev_group ] = merged_group;
        relation_groups[ curr_group ].clear();

        input_buffers_groups[ curr_group ].clear();
        input_buffers_groups[ prev_group ].clear();

        tau_score_bounds[ curr_group ] = -DBL_MAX;
        t_bounds[ curr_group ] = -DBL_MAX;



        // group_score_weights[ curr_group ] = 0;
        // group_score_weights[ prev_group ] = groupWeights( prev_group );
        max_scores_groups[ prev_group ] = prev_group_scores + curr_group_scores;
        max_scores_groups[ curr_group ] = 0;

        seen_tuples_groups[ curr_group ] = 0;
        seen_tuples_groups[ prev_group ] = 0;


        size_t j = i;           // j is the join_order_group index of curr_group.
        i = 0;                  // i is the join_order_group index of the prev_group.
        for( auto it : join_order_group ){
            if( it == prev_group )
                break;
            i++;
        }

        intermediate_results[ j ].clear();
        intermediate_results[ i ].clear();

        // recalculate dp_max_scores
        size_t k = i;
        while( k < dp_max_scores.size() ){
            // dp_max_scores[ k ] = dp_max_scores[ k-1 ] + ( max_scores_groups[ join_order_group[ k ] ] * group_score_weights[ join_order_group[ k ] ] );
            dp_max_scores[ k ] = dp_max_scores[ k-1 ] + max_scores_groups[ join_order_group[ k ] ];
            k++;
        }

        // recalculate t_bounds
        double max_score_bound = -DBL_MAX;
        for( auto it : merged_group ){
            max_score_bound = max( max_score_bound, t_relation_bounds[ it ] );
        }
        t_bounds[ prev_group ] = max_score_bound;

        // empty all the intermediate_results_groups ?
        if( i > 0 ){
          assert( j != 1 );
          for( size_t j=i-1; j < intermediate_results_groups.size(); j++ ){
              if( intermediate_results_groups[ j ].size() == 0 )
                  break;
              intermediate_results_groups[ j ].clear();
          }
        }
        else
            intermediate_results_groups[ 0 ].clear();


        // update tau_score_bounds for the prev_group.

        double score = 0; k = 0;
        while( k < join_order_group.size() ){
            if( input_buffers_groups[ join_order_group[ k ] ].size() > 0 )
                score += input_buffers_groups[ join_order_group[ k ] ][ 0 ].score - input_buffers_groups[ join_order_group[ k ] ][ 0 ].penalty;
            else
                score += max_scores_groups[ join_order_group[ k ] ];
            k++;
        }
        score = score/log_total_length;
        tau_score_bounds[ prev_group ] = score;

        cout << "Prev tau ... " << tau << "\n";
        // tau = max( tau, tau_score_bounds[ prev_group ] );
        tau = tau_score_bounds[ prev_group ];

        // update tau_score_bounds for all other groups.

        for( auto it : join_order_group ){
            if( it == prev_group || it == curr_group || tau_score_bounds[ it ] == -DBL_MAX )
                continue;
            score = tau_score_bounds[ it ];
            score = score * log_total_length;
            score -= prev_score_sum;
            score += max_scores_groups[ prev_group ];
            score = score/log_total_length;
            tau_score_bounds[ it ] = score;
            tau = max( tau, score );
        }

        cout << "New tau in regroup() ... " << tau << "\n";

    }

    void updateGroupBound( size_t i, double la_score ){
        if( la_score == -DBL_MAX ){
            tau_score_bounds[ join_order_group[ i ] ] = -DBL_MAX;
        }
        else{
            double score = la_score;
            for( auto it : join_order_group ){
                if( it == join_order_group[ i ] )
                    continue;
                size_t size = input_buffers_groups[ it ].size();
                if( size >= 1 )
                    score += ( input_buffers_groups[ it ][ 0 ].score - input_buffers_groups[ it ][ 0 ].penalty );
                else
                    // score += group_score_weights[ it ] * max_scores_groups[ it ];
                    score += max_scores_groups[ it ];
            }
            // tau_score_bounds[ join_order_group[ i ] ] = score / total_group_weights;
            tau_score_bounds[ join_order_group[ i ] ] = score / log_total_length;
        }

        double max_tau = -DBL_MAX;
        for( auto sb : tau_score_bounds )
            max_tau = max( sb.second, max_tau );
        tau = max_tau;
    }

    void accumulate_dp_max_scores( vector< double > &dp_max_scores ){
        for( size_t i=1; i < dp_max_scores.size(); i++ )
            dp_max_scores[ i ] = dp_max_scores[ i ] + dp_max_scores[ i-1 ];
    }

    void updateAllTaus( size_t i ){
        double score_to_add = input_buffers_groups[ join_order_group[ i ] ][ 0 ].score - input_buffers_groups[ join_order_group[ i ] ][ 0 ].penalty;
        // if( score_to_add == group_score_weights[ join_order_group[ i ] ] * max_scores_groups[ join_order_group[ i ] ] )
        //     return;
        if( score_to_add == max_scores_groups[ join_order_group[ i ] ] )
            return;
        cout << "updating all taus ... \n";
        cout << "\nscore to add ... " << score_to_add << "\n\n";
        cout << "\nscore to subtract ... " << max_scores_groups[ join_order_group[ i ] ] << "\n\n";
        for( auto it : join_order_group ){
            if( it == join_order_group[ i ] || tau_score_bounds[ it ] == -DBL_MAX )
                continue;
            // double score = tau_score_bounds[ it ] * total_group_weights;
            double score = tau_score_bounds[ it ] * log_total_length;
            // score -= ( group_score_weights[ join_order_group[ i ] ] * max_scores_groups[ join_order_group[ i ] ] );
            score -= max_scores_groups[ join_order_group[ i ] ];
            score += ( score_to_add );
            // score = score / total_group_weights;
            score = score / log_total_length;
            tau_score_bounds[ it ] = score;
        }
        cout << "updated all taus ... \n";
        // double max_score = group_score_weights[ join_order_group[ i ] ] * max_scores_groups[ join_order_group[ i ] ];
        double max_score = max_scores_groups[ join_order_group[ i ] ];
        for(; i < dp_max_scores.size(); i++ )
            dp_max_scores[ i ] += ( score_to_add -  max_score );
    }

    void initializeInputBuffersGroup( size_t n ){
        for( auto it : join_order_group )
            input_buffers_groups.insert( { it, {} } );
    }

    vector<IRs> initializeVectorGroups(size_t n){
        vector<IRs> intermediate_results_groups( n );
        return intermediate_results_groups;
    }

    void calculateLengths( map< SKey, SVal > query ){
        total_max_score = 0;

        for( auto it : join_order_group ){
            double max_t_group = 0;
            double max_score_group = 0;
            vector< SKey > relations_group = relation_groups[ it ];

            for( auto rc : relations_group ){
                double length = query[ rc ]->get_length();
                double score = relations[ rc ]->get_entries()[ 0 ]->score();
                lengths.insert( { rc, length } );
                cout << "log10 length of ... " << rc << "\t" << log10( length ) << "\n";
                max_t_group += ( log10( length ) * score );
                max_score_group += log10( length );
            }

            t_bounds.insert( { it, max_t_group } );

            log_total_length += max_score_group;
            total_max_score += max_t_group;
            cout << "log_total_length ... " << log_total_length << "\n";
            cout << "total_max_score ... " << total_max_score << "\n\n";

            // group_score_weights.insert( { it, groupWeights( it ) } );

            // dp_max_scores.push_back( group_score_weights[ it ] * max_t_group );
            dp_max_scores.push_back( max_t_group );

            for( auto rc : relations_group )
                t_relation_bounds.insert( { rc, max_t_group } );

            max_scores_groups.insert( { it, max_t_group } );
            // group_score_weights.insert( { it, max_score_group } );

            // total_group_weights += group_score_weights[ it ];

            // if( it == 0 )
            //     cout << "Prev score of rel group 0 inside calculateLengths() ... " << max_t_group << "\n";
            // prev_scores.insert( { it, max_t_group } );
        }
    }

    void calculateInitialScoreBounds(){
        cout << "total_max_score ... " << total_max_score << "\n";
        cout << "log_total_length ... " << log_total_length << "\n";
        // double score = total_max_score/log_total_length;
        double score = 0;
        for( auto it : join_order_group ){
            // score += ( group_score_weights[ it ] * max_scores_groups[ it ] );
            score += max_scores_groups[ it ];
        }
        for( auto it : join_order_group ){
            // tau_score_bounds.insert( { it, score/total_group_weights } );
            tau_score_bounds.insert( { it, score/log_total_length } );
            seen_tuples_groups.insert( { it, 0 } );
        }
    }

    size_t chooseInput( vector< size_t > join_order_group, map< size_t, double > score_bounds ){
        double relation_group_index = -1, max_score = -DBL_MAX;
        size_t i = 0;

        // while( output.size() < 1 && i < join_order_group.size() ){
        //     if( i == 0 && intermediate_results_groups[ i ].size() == 0 ){
        //       if( score_bounds[ join_order_group[ 0 ] ] >= score_bounds[ join_order_group[ 1 ] ] )
        //           return 0;
        //       return 1;
        //     }
        //     else if( intermediate_results_groups[ i-1 ].size() == 0 ){
        //         return i;
        //     }
        //     i++;
        // }

        if( output.size() < 1 ){
            while( i < join_order_group.size() ){
                if( i == 0 ){
                    if( intermediate_results_groups[ i ].size() == 0 ){
                        if( intermediate_results[ i ].size() == 0 )
                            return i;
                        if( intermediate_results[ i+1 ].size() == 0 )
                            return i+1;
                        if( score_bounds[ join_order_group[ 0 ] ] >= score_bounds[ join_order_group[ 1 ] ] )
                            return 0;
                        return 1;
                    }
                }
                else if( i == 1 ){
                  if( intermediate_results_groups[ i-1 ].size() == 0 ){
                      if( intermediate_results[ i ].size() == 0 )
                          return i;
                      if( score_bounds[ join_order_group[ 0 ] ] >= score_bounds[ join_order_group[ 1 ] ] )
                          return 0;
                      return 1;
                  }
                }
                else if( intermediate_results_groups[ i-1 ].size() == 0 ){
                    return i;
                }
                i++;
            }
        }
        i=0;

        while( i < join_order_group.size() ){
            if( score_bounds[ join_order_group[ i ] ] > max_score ){
                max_score = score_bounds[ join_order_group[ i ] ];
                relation_group_index = i;
            }
            i++;
        }

        return relation_group_index;
    }

    size_t chooseInput( map< size_t, double > score_bounds ){
        double relation_group = -1, max_score = -DBL_MAX;
        for( auto it : score_bounds ){
            if( it.second > max_score ){
                max_score = it.second;
                relation_group = it.first;
            }
        }
        return relation_group;
    }

    ir_tuple getNextJoinTuple( size_t i ){

        size_t group_number = join_order_group[ i ];
        vector< SKey > relations_group = relation_groups[ group_number ];
        cout << "Inside getNextJointuple() for relation group ... " << i << "\n";
        cout << "Size of relations in the group ... " << relations_group.size() << "\n";


        if( tau_score_bounds[ group_number ] == -DBL_MAX )
            return ir_tuple{ { { } }, -1, 0 };

        cout << "Initial tau score bound ... " << tau_score_bounds[ group_number ] << "\n";

        size_t index = seen_tuples_groups[ group_number ];
        size_t size = intermediate_results[ i ].size();

        // seen_tuples_groups[ group_number ] = seen_tuples_groups[ group_number ] + 1;

        cout << "Index : " << index << "\t" << "Size : " << size << "\n";

        //tau_i should be max_j t_j, for all relations j in group i.
        //when updating after each iteration, we need to update the t_relation_bound and the t_bound = max_j t_relation_bound_j.
        //t_bounds_i = sum_j t_j for all relations j in group i.
        // while( size < index + 1 || ( group_score_weights[ group_number ] * t_bounds[ group_number ] ) >= ( intermediate_results[ i ][ index ].score - intermediate_results[ i ][ index ].penalty ) ){
        while( size < index + 1 || t_bounds[ group_number ] >= ( intermediate_results[ i ][ index ].score - intermediate_results[ i ][ index ].penalty ) ){

            SKey new_relation = -1; double max_t = -DBL_MAX;
            for( auto rel : relations_group ){
                if( t_relation_bounds[ rel ] > max_t ){
                    max_t = t_relation_bounds[ rel ];
                    new_relation = rel;
                }
            }

            //if all relations in this subgroup are exhausted, then \tau_i = -\infty.
            if( new_relation == -1 || t_bounds[ group_number ] == -DBL_MAX ){
                tau_score_bounds[ group_number ] = -DBL_MAX;
                return ir_tuple{ { { } }, -1, 0 };
            }

            TVal tuple1 = getNextTuple( group_number, new_relation, t_bounds[ group_number ] );

            if( !tuple1 )
                break;

            IRs all_irs; //all join tuples for the relation group i.

            // size_t start = size;
            cout << "Tuple segment picked ... " << tuple1->segment()->get_id() << " from ... " << new_relation << "\n";
            cout << "Score of tuple picked ... " << tuple1->score() << "\n";

            map< SKey, TVal > jr = { { new_relation, tuple1 } };
            ir_tuple ir{ jr, 0, 0 };
            all_irs.push_back( ir );
            bool is_joining = true;
            for( auto rel_to_join : relations_group ){
                if( rel_to_join == new_relation )
                    continue;

                is_joining = join( all_irs, rel_to_join ); //there can be more than one join.
                if( !is_joining )
                    break;
            }
            cout << "The value of is_joining ... " << is_joining << "\n";

            cout << "The size of all_irs ... " << all_irs.size() << "\n";

            if( !is_joining )
                continue;

            for( auto it : all_irs ){
                if( it.join_result.size() == relations_group.size() ){

                    // double score = group_score_weights[ group_number ] * computeScore( it.join_result );
                    double score = computeScore( it.join_result );
                    ir_tuple res{ it.join_result, score, it.penalty };
                    cout << "total penalty ... " << it.penalty << "\n";
                    double max_score = score - it.penalty;
                    if( i > 0 )
                        max_score += dp_max_scores[ i-1 ];
                    max_score += ( dp_max_scores.back() - dp_max_scores[ i ] );
                    // max_score = max_score/total_group_weights;
                    max_score = max_score/log_total_length;
                    if( output.size() >= k && max_score <= min_output_score )
                        continue;
                    intermediate_results[ i ].push_back( res );
                    // cout << "Inside ... inserted result in intermediate_results[ i ] ... \n";
                }
            }

            // removeDuplicates( intermediate_results[ i ] );

            size = intermediate_results[ i ].size();
            cout << "Total size of results ... " << intermediate_results[ i ].size() << "\n";
            cout << "group number ... " << group_number << "\n";
            // cout << "t_bound[ i ] ... " << ( group_score_weights[ group_number ] * t_bounds[ group_number ] ) << "\n";
            cout << "t_bound[ i ] ... " << t_bounds[ group_number ] << "\n";
            cout << "Score of the result to be returned ... " << intermediate_results[ i ][ index ].score - intermediate_results[ i ][ index ].penalty << "\n";

            sort( intermediate_results[ i ].begin(), intermediate_results[ i ].end(), sort_ir );

        }
        cout << "Score of the result returned ... " << intermediate_results[ i ][ index ].score - intermediate_results[ i ][ index ].penalty << "\n";
        return intermediate_results[ i ][ index ];
    }

    //when all relations in a subgroup i are exhausted, then tau_i = -\infty.
    TVal getNextTuple( size_t i, SKey new_relation, double &t ){
        TVal tuple1;
        size_t index = 0;
        if( seen_tuples.find( new_relation )->second == size_of_relations.find( new_relation )->second )
            return tuple1;
        if ( seen_tuples.find( new_relation ) == seen_tuples.end() ){
            tuple1 = relations.find( new_relation )->second->get_entries().front();
            seen_tuples.insert( { new_relation, 1 } );
            input_buffers.insert( { new_relation, { tuple1 } } );
        }
        else{
            index = seen_tuples.find( new_relation )->second;
            tuple1 = relations.find( new_relation )->second->get_entries().at( index );
            seen_tuples.find( new_relation )->second = index + 1;
            input_buffers.find( new_relation )->second.push_back( tuple1 );
        }
        index = seen_tuples.find( new_relation )->second;
        if(index == size_of_relations.find( new_relation )->second){
            exhausted_relations.insert( new_relation );
            t = updateScoreBounds( i, new_relation, -DBL_MAX, tuple1->score(), lengths );
        }
        else{
            TVal look_ahead_tuple = relations.find( new_relation )->second->get_entries()[ index ];
            t = updateScoreBounds( i, new_relation, look_ahead_tuple->score(), tuple1->score(), lengths );
        }
        cout << "The value of t inside getNextTuple() ... " << t << "\n";
        return tuple1;
    }

    bool join( IRs &all_irs, SKey rel_to_join ){
        bool join_found = false;
        IRs all_irs_copy = all_irs;
        all_irs.clear();
        vector< TVal > all_tuples = input_buffers.find( rel_to_join )->second;
        if( all_tuples.size() == 0 ){
            return false;
        }
        for( auto ir : all_irs_copy ){
            for( auto tuple : all_tuples ){
                double penalty = 0;
                if( joinsAll( tuple, ir.join_result, rel_to_join, joins, penalty ) ){
                    map< SKey, TVal > new_ir = ir.join_result;
                    new_ir.insert( { rel_to_join, tuple } );
                    ir_tuple res{ new_ir, computeScore( new_ir ), ir.penalty + penalty };
                    cout << "inside joinAll() of join() ... " << "\n";
                    all_irs.push_back( res );
                    join_found = true;
                }
            }
        }
        return join_found;
    }


    bool joinsWith( map< SKey, TVal > prev_ir, map< SKey, TVal > ir, size_t i, double &penalty, map< size_t, size_t > &group_failing ){
        // cout << "Inside new joinsWith() ... \n";
        // map< size_t, size_t > group_failing; // count (value) of the relations from the group number (key) failing to join with ir.
        map< SKey, TVal > new_ir = prev_ir;
        for( auto tup : ir ){
            SKey id = tup.first;
            TVal tuple1 = tup.second;
            vector< SKey > all_joins = joins.find( id )->second;
            for( auto rel : all_joins ){
                if( prev_ir.find( rel ) != prev_ir.end() ){
                    // if( !joinsAll( tuple1, prev_ir, id, joins ) )
                    //     return false;
                    if( !isAJoin( tuple1, prev_ir[ rel ], penalty, id, rel ) ){
                        // if( i == 158 )
                            // cout << "Join failed in joinsWith() for relations ... " << rel << "\t" << id << "\n";
                        SKey grp_number = group_number[ rel ];
                        if( group_failing.find( grp_number ) != group_failing.end() )
                            group_failing[ grp_number ] += 1;
                        else
                            group_failing[ grp_number ] = 1;
                        cout << "Join failed at ... " << rel << "\t" << id << "\t" << prev_ir[ rel ]->segment()->get_id() << "\t" << tuple1->segment()->get_id() << "\n";
                        return false;
                    }
                }
            }
            new_ir.insert( { id, tuple1 } );
            // if( i == 49 )
            //     cout << "Inserted a tuple in joinsWith() ... \n";
        }
        if( new_ir.size() == prev_ir.size() + ir.size() )
            return true;
    }

    bool joinsWith( map< SKey, TVal > prev_ir, map< SKey, TVal > ir, size_t i, double &penalty ){
        // cout << "Inside new joinsWith() ... \n";
        // map< size_t, size_t > group_failing; // count (value) of the relations from the group number (key) failing to join with ir.
        map< SKey, TVal > new_ir = prev_ir;
        for( auto tup : ir ){
            SKey id = tup.first;
            TVal tuple1 = tup.second;
            vector< SKey > all_joins = joins.find( id )->second;
            for( auto rel : all_joins ){
                if( prev_ir.find( rel ) != prev_ir.end() ){
                    // if( !joinsAll( tuple1, prev_ir, id, joins ) )
                    //     return false;
                    if( !isAJoin( tuple1, prev_ir[ rel ], penalty, id, rel ) ){
                        // if( i == 158 )
                            // cout << "Join failed in joinsWith() for relations ... " << rel << "\t" << id << "\n";
                        cout << "Join failed at ... " << rel << "\t" << id << "\t" << prev_ir[ rel ]->segment()->get_id() << "\t" << tuple1->segment()->get_id() << "\n";
                        return false;
                    }
                }
            }
            new_ir.insert( { id, tuple1 } );
            // if( i == 49 )
            //     cout << "Inserted a tuple in joinsWith() ... \n";
        }
        if( new_ir.size() == prev_ir.size() + ir.size() )
            return true;
    }


    double computeScore( map< SKey, TVal > ir){
        double score = 0;
        for( auto it : ir ){
            SKey id = it.first;
            double length = lengths[ id ];
            score += ( log10( length ) * it.second->score() );
        }
        return score;
    }

    void findSingleMatches(set<SKey> &single_matches, map<SKey, RVal> &relations){
        for(auto it : relations){
            SKey id = it.first;
            RVal rel = it.second;
            if(rel->get_entries().size() == 1 && rel->get_entries()[ 0 ]->segment()->get_id() == 0){
                single_matches.insert(id);
                dummy = rel->get_entries()[0];
            }
        }
        for(auto it : single_matches){
            relations.erase(it);
        }
    }

    void completeJoin( size_t i ) {
        for(auto item = join_order_group.begin() + i; item != join_order_group.end(); ++i, ++item){
            size_t start_index = start_indices[ i-2 ];
            start_indices[ i-1 ] = intermediate_results_groups[ i-1 ].size();
            IRs prev_irs = intermediate_results_groups[ i-2 ];
            vector< SKey > relations_group = relation_groups.find( *item )->second;

            if( regrouped_groups.find( *item ) != regrouped_groups.end() ){
                intermediate_results_groups[ i-1 ].insert( intermediate_results_groups[ i-1 ].end(), prev_irs.begin() + start_index, prev_irs.end() );
                removeDuplicates( intermediate_results_groups[ i-1 ] );
                continue;
            }
            if( input_buffers_groups[ *item ].size() == 0 )
                return;
            IRs all_tuples = input_buffers_groups[ *item ];

            assert(start_index <= prev_irs.size());

            for( auto prev_ir = prev_irs.begin() + start_index; prev_ir != prev_irs.end(); ++prev_ir ){
                for( auto ir : all_tuples ){
                    double penalty = 0;
                    if( joinsWith( ( *prev_ir ).join_result, ir.join_result, *item, penalty ) ){
                        map< SKey, TVal > new_ir = ( *prev_ir ).join_result;

                        double score = ( *prev_ir ).score + ir.score;
                        double total_score = score + dp_max_scores.back() - dp_max_scores[ i ];
                        double total_penalty = ( *prev_ir ).penalty + ir.penalty + penalty;
                        // if( output.size() >= k && ( total_score - total_penalty )/total_group_weights <= min_output_score ){
                        if( output.size() >= k && ( total_score - total_penalty )/log_total_length <= min_output_score ){
                            cout << "in complete join () ... " << min_output_score << "\t" << ( total_score - total_penalty )/log_total_length << "\n";
                            continue;
                        }

                        // cout << "Max. possible score of intermediate result ... " << ( total_score - total_penalty )/log_total_length << "\n";
                        // cout << "Minimum output score ... " << min_output_score << "\n";
                        new_ir.insert( ir.join_result.begin(), ir.join_result.end() );
                        ir_tuple res{ new_ir, score, total_penalty };
                        intermediate_results_groups[ i-1 ].push_back( res );
                    }
                }
            }
            removeDuplicates( intermediate_results_groups[ i-1 ] );
        }
    }

    void formOutputs( vector< JoinResult > &output){
        cout << "Intermediate results groups back size() in formOutput() ... " << intermediate_results_groups.back().size() << "\n";
        IRs results = intermediate_results_groups.back();
        size_t start_index = start_indices[ n-2 ];
        size_t initial_size = output.size();
        cout << "Start index in formOutputs() ... " << start_index << "\n";
        for( auto it = results.begin() + start_index; it != results.end(); ++it ){
            map< SKey, TVal > res = ( *it ).join_result;
            // assert( res.size() == joins.size() );
            cout << "Join result size in the formOutputs() ... " << res.size() << "\n";
            double score = ( *it ).score;
            double penalty = ( *it ).penalty;
            cout << "Penalty in formOutput() ... " << penalty << "\n";
            JoinResult result;
            result.result = res;
            // result.score = ( score - penalty ) / total_group_weights;
            result.score = ( score - penalty ) / log_total_length;
            output.push_back( result );
        }
        // if( initial_size == 0 && output.size() > 0 )
        //     accumulate_dp_max_scores( dp_max_scores );
        start_indices[ n-2 ] = results.size();
        sort( output.begin(), output.end(), sort_output );
    }

    bool joinsAll(TVal tuple, map<SKey, TVal> ir, SKey relation, map<SKey, vector<SKey>> adjlist, double &penalty ){
        bool flag = true;
        for( auto v : adjlist.find( relation )->second ){
            if( ir.find(v) != ir.end() ){
                if( !isAJoin( ir.find(v)->second, tuple, penalty, v, relation ) ){
                    flag = false;
                    break;
                }
            }
        }

        return flag;
    }

    bool present(TVal tuple, map<SKey, TVal> ir){
        return false;
        if(tuple->segment()->get_id() == 0)
            return false;
        for(auto& i : ir){
            if(i.second->segment()->get_id() == tuple->segment()->get_id())
                return true;
        }
        return false;
    }

    double updateScoreBounds( size_t i, SKey new_relation, double score, double prev_score, map<SKey, double> &lengths){
        double length = lengths[ new_relation ];
        double new_score = -DBL_MAX;
        double current_score = t_relation_bounds.find( new_relation )->second;
        if(score != -DBL_MAX)
            new_score = current_score + ( log10( length ) * ( score - prev_score ) );

        t_relation_bounds[ new_relation ] = new_score;
        assert( new_score <= current_score );

        double new_bound = -DBL_MAX;
        vector< SKey > relations_group = relation_groups[ i ];
        for( auto sb : relations_group )
            new_bound = max( new_bound, t_relation_bounds[ sb ] );

        return new_bound;
    }

    vector<IRs> initializeVector(size_t n){
        vector<IRs> intermediate_results(n);
        return intermediate_results;
    }

    vector< size_t > initializeStartIndices(size_t n){
        vector< size_t > start_indices;
        for(size_t i=0; i<n; i++)
            start_indices.push_back( 0 );
        return start_indices;
    }

    void removeDuplicates(IRs &ir){
        if( ir.size() <= 1 )
            return;
        vector<size_t> toErase;
        size_t i = 0;
        for(auto it1 = ir.begin(); it1 != ir.end(); ++it1, ++i){
            for(auto it2 = ir.begin(); it2 < it1; ++it2){
                if( ( *it1 ).score == ( *it2 ).score && map_compare( ( *it1 ).join_result, ( *it2 ).join_result ) ){
                    toErase.push_back( i );
                    break;
                }
            }
        }
        cout << "duplicates size ... " << toErase.size() << "\n";
        // if( toErase.size() > 1 ){
        //     for( auto i : toErase ){
        //         auto itr = ir[ toErase[ i ] ];
        //         for( auto it : itr.join_result ){
        //             cout << it.first << ": " << it.second->segment()->get_id() << "\n";
        //         }
        //         cout << "----------\n\n";
        //     }
        //     exit( EXIT_SUCCESS );
        // }
        size_t cnt = 0;
        for( auto i=0; i<toErase.size(); ++i, ++cnt )
            ir.erase( ir.begin() + toErase[ i-cnt ] );
    }

    map<SKey, size_t> getSize(map<SKey, RVal> &relations){
        map<SKey, size_t> size_of_relations;
        for(auto& it : relations){
            SKey id = it.first;
            size_t n = it.second->get_entries().size();
            size_of_relations.insert({id, n});
        }
        return size_of_relations;
    }

    bool isAJoin(TVal &tuple1, TVal &tuple2, double &penalty, SKey relation1, SKey relation2 ){
        if(tuple1->segment()->get_id() == 0 && tuple2->segment()->get_id() == 0){
            // cout << "passed ... 1 \n";
            return true;
        }

        // cout << "failed ... 1 \n";
        if(tuple1->segment()->get_id() == tuple2->segment()->get_id() && (tuple1->segment()->get_length() <= ce || tuple2->segment()->get_length() <= ce)){
          // cout << "passed ... 2 \n";
          return true;
        }
        //     return true;
        // cout << "failed ... 2 \n";
        if(tuple1->segment()->is_roundabout() && tuple2->segment()->is_roundabout()){
          // cout << "passed ... 3 \n";
          return true;
        }
        //     return true;
        // cout << "failed ... 3 \n";
        if(tuple1->segment()->get_id() == 0 || tuple2->segment()->get_id() == 0){
          // cout << "passed ... 4 \n";
          return true;
        }
        //     return true;
        // cout << "failed ... 4 \n";
        if(tuple1->segment()->is_roundabout() || tuple2->segment()->is_roundabout()){

            SVal seg1, seg2;
            if(tuple1->segment()->is_roundabout()){
                seg1 = tuple1->segment();
                seg2 = tuple2->segment();
            }
            else{
                seg2 = tuple1->segment();
                seg1 = tuple2->segment();
            }

            for(auto& it : seg1->get_intermediate_nodes()){
                if(sameLocation(it.get_location(), seg2->get_first()->get_location()) || sameLocation(it.get_location(), seg2->get_second()->get_location())){
                  // cout << "passed ... 5 \n";
                  return true;
                }
                    // return true;
            }
            // cout << "failed ... 5 \n";
            return false;
        }
        if(tuple1->get_first()->get_id() == tuple2->get_first()->get_id()){
          // cout << "passed ... 6 \n";
          return true;
        }
        //     return true;
        // cout << "failed ... 6 \n";
        if(tuple1->get_second()->get_id() == tuple2->get_second()->get_id()){
          // cout << "passed ... 7 \n";
          return true;
        }
        //     return true;
        // cout << "failed ... 7 \n";
        if(tuple1->get_first()->get_id() == tuple2->get_second()->get_id()){
          // cout << "passed ... 8 \n";
          return true;
        }
        //     return true;
        // cout << "failed ... 8 \n";
        if(tuple1->get_second()->get_id() == tuple2->get_first()->get_id()){
          // cout << "passed ... 9 \n";
          return true;
        }
        //     return true;
        // cout << "failed ... 9 \n";
        if(possible_approximate( tuple1, relation1, tuple2, relation2 ) && isRelaxedJoin(tuple1, tuple2)){
            return true;
        }

        if( possible_approximate( tuple1, relation1, tuple2, relation2 ) && isApproximateJoin( tuple1, tuple2, penalty ) ){
            return true;
        }
        return false;
    }

    //approximate join is feasible only when there are no higher scored tuples that join "strictly" with each other in the corresponding relations.
    bool possible_approximate( TVal tuple1, SKey relation1, TVal tuple2, SKey relation2 ){
        // return true;
        vector< TVal > all_tuples1 = relations[ relation1 ]->get_entries();
        vector< TVal > all_tuples2 = relations[ relation2 ]->get_entries();
        for( auto tup1 : all_tuples1 ){
            if( tup1->segment()->get_id() == tuple1->segment()->get_id() )
                break;
            if( tup1->segment()->get_first()->get_id() == tuple2->segment()->get_second()->get_id() )
                return false;
            if( tup1->segment()->get_first()->get_id() == tuple2->segment()->get_first()->get_id() )
                return false;
            if( tup1->segment()->get_second()->get_id() == tuple2->segment()->get_second()->get_id() )
                return false;
            if( tup1->segment()->get_second()->get_id() == tuple2->segment()->get_first()->get_id() )
                return false;
        }
        for( auto tup2 : all_tuples2 ){
            if( tup2->segment()->get_id() == tuple2->segment()->get_id() )
                break;
            if( tup2->segment()->get_first()->get_id() == tuple1->segment()->get_second()->get_id() )
                return false;
            if( tup2->segment()->get_first()->get_id() == tuple1->segment()->get_first()->get_id() )
                return false;
            if( tup2->segment()->get_second()->get_id() == tuple1->segment()->get_second()->get_id() )
                return false;
            if( tup2->segment()->get_second()->get_id() == tuple1->segment()->get_first()->get_id() )
                return false;
        }
        // for( auto tup1 : all_tuples1 ){
        //       if( tup1->segment()->get_id() == tuple1->segment()->get_id() )
        //           break;
        //       for( auto tup2 : all_tuples2 ){
        //             if( tup2->segment()->get_id() == tuple2->segment()->get_id() )
        //                 break;
        //             if( tup2->segment()->get_first()->get_id() == tuple1->segment()->get_second()->get_id() )
        //                 return false;
        //             if( tup2->segment()->get_first()->get_id() == tuple1->segment()->get_first()->get_id() )
        //                 return false;
        //             if( tup2->segment()->get_second()->get_id() == tuple1->segment()->get_second()->get_id() )
        //                 return false;
        //             if( tup2->segment()->get_second()->get_id() == tuple1->segment()->get_first()->get_id() )
        //                 return false;
        //       }
        // }
        return true;
    }

    bool isApproximateJoin( TVal tuple1, TVal tuple2, double &penalty ){

        Coordinates c1 = tuple1->segment()->get_first()->get_location();
        Coordinates c2 = tuple1->segment()->get_second()->get_location();
        Coordinates c3 = tuple2->segment()->get_first()->get_location();
        Coordinates c4 = tuple2->segment()->get_second()->get_location();

        double d11 = haversine::distance( c1, c3 );
        double d22 = haversine::distance( c2, c4 );
        double d12 = haversine::distance( c1, c4 );
        double d21 = haversine::distance( c2, c3 );

        size_t cases = 5;

        if( d11 <= ce && d11 < d22 && d11 < d12 && d11 < d21 )
            cases = 1;

        else if( d12 <= ce && d12 < d11 && d12 < d21 && d12 < d22 )
            cases = 2;

        else if( d21 <= ce && d21 < d11 && d21 < d12 && d21 < d22 )
            cases = 3;

        else if( d22 <= ce && d22 < d11 && d22 < d12 && d22 < d21 )
            cases = 4;

        cout << "case in approx. join ... " << cases << " " << ce << " " << d11 << " " << d12 << " " << d21 << " " << d22 << "\n";

        if( cases == 5 )
            return false;

        Coordinates c1_m = lonlat_to_mercator( c1 );
        Coordinates c2_m = lonlat_to_mercator( c2 );
        Coordinates c3_m = lonlat_to_mercator( c3 );
        Coordinates c4_m = lonlat_to_mercator( c4 );

        const Geometry* p1 = global_factory->createPoint( Coordinate( pm->makePrecise( c1_m.x ), pm->makePrecise( c1_m.y ) ) );
        const Geometry* p2 = global_factory->createPoint( Coordinate( pm->makePrecise( c2_m.x ), pm->makePrecise( c2_m.y ) ) );
        const Geometry* p3 = global_factory->createPoint( Coordinate( pm->makePrecise( c3_m.x ), pm->makePrecise( c3_m.y ) ) );
        const Geometry* p4 = global_factory->createPoint( Coordinate( pm->makePrecise( c4_m.x ), pm->makePrecise( c4_m.y ) ) );

        Geometry *g_intersection, *g_union, *buff1, *buff2;
        SKey node_id1, node_id2;

        if( cases == 1 ){
            buff1 = p1->buffer( ce );
            buff2 = p3->buffer( ce );
            node_id1 = tuple1->segment()->get_first()->get_id();
            node_id2 = tuple2->segment()->get_first()->get_id();
        }
        else if( cases == 2 ){
            buff1 = p1->buffer( ce );
            buff2 = p4->buffer( ce );
            node_id1 = tuple1->segment()->get_first()->get_id();
            node_id2 = tuple2->segment()->get_second()->get_id();
        }
        else if ( cases == 3 ){
            buff1 = p2->buffer( ce );
            buff2 = p3->buffer( ce );
            node_id1 = tuple1->segment()->get_second()->get_id();
            node_id2 = tuple2->segment()->get_first()->get_id();
        }
        else{
            buff1 = p2->buffer( ce );
            buff2 = p4->buffer( ce );
            node_id1 = tuple1->segment()->get_second()->get_id();
            node_id2 = tuple2->segment()->get_second()->get_id();
        }

        // if( node_id1 < 0 || node_id2 < 0 )
        //     return false;

        if(end_points_segments.find( node_id1 ) != end_points_segments.end()){
            if(end_points_segments[ node_id1 ].find( node_id2 ) != end_points_segments[ node_id1 ].end() ){
               cout << "approximate join avoided ... \n";
               return false;
            }

        }
        if(end_points_segments.find( node_id2 ) != end_points_segments.end()){
            if(end_points_segments[ node_id2 ].find( node_id1 ) != end_points_segments[ node_id2 ].end() ){
              cout << "approximate join avoided ... \n";
              return false;
            }
        }

        if(relaxed_joins.find( node_id1 ) != relaxed_joins.end()){
            if(relaxed_joins[ node_id1 ] == node_id2 ){
              cout << "approximate join avoided ... \n";
              return false;
            }
        }
        if(relaxed_joins.find( node_id2 ) != relaxed_joins.end()){
            if(relaxed_joins[ node_id2 ] == node_id1 ){
              cout << "approximate join avoided ... \n";
              return false;
            }
        }

        // if( approximate_joins.find( node_id1 ) != approximate_joins.end() ){
        //     if( approximate_joins[ node_id1 ] == node_id2 )
        //       return true;
        // }
        //
        // if( approximate_joins.find( node_id2 ) != approximate_joins.end() ){
        //     if( approximate_joins[ node_id2 ] == node_id1 )
        //       return true;
        // }

        g_intersection = buff1->intersection( buff2 );
        g_union = buff1->Union( buff2 );

        double overlap = g_intersection->getArea()/g_union->getArea();
        // cout << "Previous penalty ... " << penalty << "\n";
        // penalty += ( 1 - overlap ) * ( 1 - overlap );
        // cout << "Returning true in approximate join ... " << penalty << "\n";
        // cout << "joining ... " << tuple1->segment()->get_id() << " " << tuple1->segment()->get_first()->get_id() << " " << tuple1->segment()->get_second()->get_id() << "\n";
        // cout << "joining ... " << tuple2->segment()->get_id() << " " << tuple2->segment()->get_first()->get_id() << " " << tuple2->segment()->get_second()->get_id() << "\n";
        // if( approximate_joins.find( tuple1->segment()->get_id() ) == approximate_joins.end() )
        //     approximate_joins[ tuple1->segment()->get_id() ] = {};
        // approximate_joins[ tuple1->segment()->get_id() ].push_back( tuple2->segment()->get_id() );
        approximate_joins.insert( { node_id1, node_id2 } );
        penalty = 0;
        return true;
    }

    bool isRelaxedJoin(TVal tuple1, TVal tuple2){
        if(relaxed_joins.find(tuple1->get_first()->get_id()) == relaxed_joins.end() && relaxed_joins.find(tuple1->get_second()->get_id()) == relaxed_joins.end())
            return false;
        if(relaxed_joins.find(tuple2->get_first()->get_id()) == relaxed_joins.end() && relaxed_joins.find(tuple2->get_second()->get_id()) == relaxed_joins.end())
            return false;
        if(relaxed_joins.find(tuple1->get_first()->get_id())->second == tuple2->get_second()->get_id())
            return true;
        if(relaxed_joins.find(tuple1->get_first()->get_id())->second == tuple2->get_first()->get_id())
            return true;
        if(relaxed_joins.find(tuple1->get_second()->get_id())->second == tuple2->get_second()->get_id())
            return true;
        if(relaxed_joins.find(tuple1->get_second()->get_id())->second == tuple2->get_first()->get_id())
            return true;
        return false;
    }

    void formQueryGroupsUsingBFS( map< SKey, vector<SKey> > joins, size_t count ){
        set< SKey > visited;
        // map< SKey, size_t > group_number; //group number that each relation belongs to.
        size_t index_of_group = 0;
        list< SKey > queue;
        for( auto ec : edge_cases )
            visited.insert( ec );

        count = count - edge_cases.size();
        for(auto& it : joins){
            if( count <= 0 )
                break;

            SKey start = it.first;
                                                                                             //if the segment is not conencted to any segment, the corresponding list in the join will be empty. Insert "-1" for such cases.
            if( visited.find( start ) == visited.end() ){
                visited.insert( start );
                group_number.insert( { start, index_of_group } );
                count--;                                                  //add the "start" vertex because it will not be added during the BFS traversal in the query_order.

                if( it.second.size() == 0 ){
                    query_graph_group.insert( { index_of_group, -1 } );
                    relation_groups.insert( { index_of_group, { start } } );
                    join_order_group.push_back( index_of_group );
                    index_of_group++;
                    continue;
                }                                                       //isolated segment.

                vector< SKey > group = { start };

                for( auto rc : joins.find( start )->second ){
                    if( visited.find( rc ) != visited.end() )
                        continue;
                    group_number.insert( { rc, index_of_group } );
                    queue.push_back( rc );
                    group.push_back( rc );
                    visited.insert( rc );
                }
                query_graph_group.insert( { index_of_group, index_of_group + 1 } );
                relation_groups.insert( { index_of_group, group } );
                join_order_group.push_back( index_of_group++ );
                group.clear();

                while( !queue.empty() ){

                    if( count == 0 )
                        break;
                    SKey s = queue.front();
                    queue.pop_front();

                    size_t index_of_group_to_join = group_number.find( s )->second;
                    for( auto rc : joins.find( s )->second ){
                        if( visited.find( rc ) != visited.end() )
                            continue;

                        group.push_back( rc );
                        visited.insert( rc );
                        for( auto rc2 : joins[ rc ] ){
                           if( visited.find( rc2 ) != visited.end() )
                              continue;

                            group_number.insert( { rc2, index_of_group } );
                            queue.push_back( rc2 );
                            group.push_back( rc2 );
                            visited.insert( rc2 );
                            count--;
                            if( count == 0 )
                                break;
                        }
                        bool inserted = false;
                        if( group.size() == 1 ){
                            for( auto rel_join : joins[ group[ 0 ] ] ){
                                if( group_number.find( rel_join ) != group_number.end() ){
                                    size_t index_rel_join = group_number[ rel_join ];
                                    relation_groups[ index_rel_join ].push_back( group[ 0 ] );
                                    inserted = true;
                                    group_number.insert( { group[ 0 ], index_rel_join } );
                                    break;
                                }
                            }
                        }
                        if( group.size() > 0 && !inserted ){
                            query_graph_group.insert( { index_of_group, index_of_group_to_join } );
                            relation_groups.insert( { index_of_group, group } );
                            join_order_group.push_back( index_of_group++ );
                        }
                        group.clear();
                    }
                }
            }
        }
    }

    void formFixedSizeBFSGroups( map< SKey, vector< SKey > > joins, size_t size ){
        set< SKey > visited;
        // map< SKey, size_t > group_number;

        for( auto ec : edge_cases )
            visited.insert( ec );

        size_t index_of_group = 0;

        for( auto it : joins ){
            SKey id = it.first;

            if( visited.find( id ) != visited.end() )
                continue;

            visited.insert( id );

            list< SKey > queue;

            vector< SKey > group = { id };
            visited.insert( id );
            size_t ctr = size-1;
            queue.push_back( id );

            while( ctr > 0 && !queue.empty() ){
                SKey start = queue.front();
                queue.pop_front();

                for( auto rel_connected : joins[ start ] ){
                    if( visited.find( rel_connected ) != visited.end() )
                        continue;
                    group.push_back( rel_connected );
                    visited.insert( rel_connected );
                    queue.push_back( rel_connected );
                    ctr--;
                    if( ctr == 0 )
                        break;
                }
            }

            bool found = false;
            if( group.size() == 1 ){
                for( auto rel_in_group : joins[ group[ 0 ] ] ){
                    if( group_number.find( rel_in_group ) != group_number.end() ){
                        found = true;
                        relation_groups[ group_number[ rel_in_group ] ].push_back( group[ 0 ] );
                        group_number.insert( { group[ 0 ], group_number[ rel_in_group ] } );
                        break;
                    }
                }
            }
            if( group.size() > 0 && !found ){
                for( auto group_rel : group )
                    group_number.insert( { group_rel, index_of_group } );
                relation_groups.insert( { index_of_group++, group } );
            }
        }

        list< SKey > queue;
        visited.clear();

        for( auto it1 : relation_groups ){
            if( find( join_order_group.begin(), join_order_group.end(), it1.first ) != join_order_group.end() )
                continue;
            join_order_group.push_back( it1.first );
            for( auto it2 : it1.second ){
                visited.insert( it2 );
                for( auto it3: joins[ it2 ] )
                    queue.push_back( it3 );
            }

            while( !queue.empty() ){
                SKey id = queue.front();
                visited.insert( id );
                queue.pop_front();
                SKey grp_number = group_number[ id ];
                if( find( join_order_group.begin(), join_order_group.end(), grp_number ) == join_order_group.end() )
                    join_order_group.push_back( grp_number );
                for( auto it2 : joins[ id ] ){
                    if( visited.find( it2 ) == visited.end() )
                        queue.push_back( it2 );
                }
            }
        }

    }

    void formStarGroups( map< SKey, vector<SKey> > joins ){

        set< SKey > visited;
        // map< SKey, size_t > group_number;

        for( auto ec : edge_cases )
            visited.insert( ec );

        for( auto it : joins ){
            if( visited.find( it.first ) != visited.end() )
                continue;

            size_t ctr = 0;
            for( auto it2 : joins[ it.first ] ){
                if( visited.find( it2 ) != visited.end() )
                    ctr++;
            }
            degree.insert( { it.first, it.second.size() - ctr } );
            assert( it.second.size() - ctr >= 0 );
        }

        size_t index_of_group = 0;

        while( degree.size() > 0 ){
            SKey id = findLargest( degree );

            if( id == -1 )
                break;
            if( visited.find( id ) != visited.end() )
                continue;
            vector< SKey > group = { id };
            visited.insert( id );

            for( auto rel_connected : joins[ id ] ) {
                if( visited.find( rel_connected ) != visited.end() || degree.find( rel_connected ) == degree.end() || degree[ rel_connected ] == 0 )
                    continue;
                group.push_back( rel_connected );
                visited.insert( rel_connected );

                for( auto degree_dec_rel : joins[ rel_connected ] ) {
                    if( degree.find( degree_dec_rel ) == degree.end() || visited.find( degree_dec_rel ) != visited.end() )
                        continue;
                    degree[ degree_dec_rel ]--;
                    if( degree[ degree_dec_rel ] == 0 && visited.find( degree_dec_rel ) != visited.end() )
                        degree.erase( degree_dec_rel );
                }
                degree.erase( rel_connected );
            }
            //comment below the if-block if you do not want to put relation from singleton group in some other group.
            bool found = false;
            if( group.size() == 1 ){
                for( auto rel_in_group : joins[ group[ 0 ] ] ){
                    if( group_number.find( rel_in_group ) != group_number.end() ){
                        found = true;
                        relation_groups[ group_number[ rel_in_group ] ].push_back( group[ 0 ] );
                        group_number.insert( { group[ 0 ], group_number[ rel_in_group ] } );
                        break;
                    }
                }
            }
            if( group.size() > 0 && !found ){
                for( auto group_rel : group )
                    group_number.insert( { group_rel, index_of_group } );
                relation_groups.insert( { index_of_group++, group } );
            }
            degree.erase( id );
        }

        // put the singleton groups in some other group.

        list< SKey > queue;
        visited.clear();

        for( auto it1 : relation_groups ){
            if( find( join_order_group.begin(), join_order_group.end(), it1.first ) != join_order_group.end() )
                continue;
            join_order_group.push_back( it1.first );
            for( auto it2 : it1.second ){
                visited.insert( it2 );
                for( auto it3: joins[ it2 ] )
                    queue.push_back( it3 );
            }

            while( !queue.empty() ){
                SKey id = queue.front();
                visited.insert( id );
                queue.pop_front();
                SKey grp_number = group_number[ id ];
                if( find( join_order_group.begin(), join_order_group.end(), grp_number ) == join_order_group.end() )
                    join_order_group.push_back( grp_number );
                for( auto it2 : joins[ id ] ){
                    if( visited.find( it2 ) == visited.end() )
                        queue.push_back( it2 );
                }
            }
        }

        cout << "Size of join_order_group in star groups ... " << join_order_group.size() << "\n";
        // size_t gn = group_number[ 1327 ];
        // cout << "Verification ... " << gn << "\t" << relation_groups[ gn ].size() << "\n";
    }

    void unitSizeGroups( map< SKey, vector< SKey > > joins ){
      set< SKey > visited;
      size_t index_of_group = 0;

      for( auto it : joins ){
          SKey id = it.first;
          if( visited.find( id ) != visited.end() )
              continue;

          visited.insert( id );
          vector< SKey > group = { id };
          group_number.insert( { id, index_of_group } );
          relation_groups.insert( { index_of_group++, group } );
      }

      list< SKey > queue;
      visited.clear();

      for( auto it1 : relation_groups ){
          if( find( join_order_group.begin(), join_order_group.end(), it1.first ) != join_order_group.end() )
              continue;
          join_order_group.push_back( it1.first );
          for( auto it2 : it1.second ){
              visited.insert( it2 );
              for( auto it3: joins[ it2 ] )
                  queue.push_back( it3 );
          }

          while( !queue.empty() ){
              SKey id = queue.front();
              visited.insert( id );
              queue.pop_front();
              SKey grp_number = group_number[ id ];
              if( find( join_order_group.begin(), join_order_group.end(), grp_number ) == join_order_group.end() )
                  join_order_group.push_back( grp_number );
              for( auto it2 : joins[ id ] ){
                  if( visited.find( it2 ) == visited.end() )
                      queue.push_back( it2 );
              }
          }
      }
    }

    void formFixedSizeStarGroups( map< SKey, vector<SKey> > joins, size_t size ){

        set< SKey > visited;
        // map< SKey, size_t > group_number;

        for( auto ec : edge_cases )
            visited.insert( ec );

        for( auto it : joins ){
            if( visited.find( it.first ) != visited.end() )
                continue;

            size_t ctr = 0;
            for( auto it2 : joins[ it.first ] ){
                if( visited.find( it2 ) != visited.end() )
                    ctr++;
            }
            degree.insert( { it.first, it.second.size() - ctr } );
            assert( it.second.size() - ctr >= 0 );
        }

        size_t index_of_group = 0;

        while( degree.size() > 0 ){
            SKey id = findLargest( degree );

            if( id == -1 )
                break;
            if( visited.find( id ) != visited.end() )
                continue;
            vector< SKey > group = { id };
            visited.insert( id );
            size_t ctr = size-1;

            list< SKey > queue;
            queue.push_back( id );

            while( ctr > 0 && !queue.empty() ){
                SKey start = queue.front();
                queue.pop_front();
                for( auto rel_connected : joins[ start ] ) {
                    if( visited.find( rel_connected ) != visited.end() || degree.find( rel_connected ) == degree.end() || degree[ rel_connected ] == 0 )
                        continue;
                    group.push_back( rel_connected );
                    visited.insert( rel_connected );
                    queue.push_back( rel_connected );

                    for( auto degree_dec_rel : joins[ rel_connected ] ) {
                        if( degree.find( degree_dec_rel ) == degree.end() || visited.find( degree_dec_rel ) != visited.end() )
                            continue;
                        degree[ degree_dec_rel ]--;
                        if( degree[ degree_dec_rel ] == 0 && visited.find( degree_dec_rel ) != visited.end() )
                            degree.erase( degree_dec_rel );
                    }
                    degree.erase( rel_connected );
                    ctr--;
                    if( ctr == 0 )
                        break;
                }
            }

            //comment below the if-block if you do not want to put relation from singleton group in some other group.
            bool found = false;
            if( group.size() == 1 ){
                for( auto rel_in_group : joins[ group[ 0 ] ] ){
                    if( group_number.find( rel_in_group ) != group_number.end() ){
                        found = true;
                        relation_groups[ group_number[ rel_in_group ] ].push_back( group[ 0 ] );
                        group_number.insert( { group[ 0 ], group_number[ rel_in_group ] } );
                        break;
                    }
                }
            }
            if( group.size() > 0 && !found ){
                for( auto group_rel : group )
                    group_number.insert( { group_rel, index_of_group } );
                relation_groups.insert( { index_of_group++, group } );
            }
            degree.erase( id );
        }

        // put the singleton groups in some other group.

        list< SKey > queue;
        visited.clear();

        for( auto it1 : relation_groups ){
            if( find( join_order_group.begin(), join_order_group.end(), it1.first ) != join_order_group.end() )
                continue;
            join_order_group.push_back( it1.first );
            for( auto it2 : it1.second ){
                visited.insert( it2 );
                for( auto it3: joins[ it2 ] )
                    queue.push_back( it3 );
            }

            while( !queue.empty() ){
                SKey id = queue.front();
                visited.insert( id );
                queue.pop_front();
                SKey grp_number = group_number[ id ];
                if( find( join_order_group.begin(), join_order_group.end(), grp_number ) == join_order_group.end() )
                    join_order_group.push_back( grp_number );
                for( auto it2 : joins[ id ] ){
                    if( visited.find( it2 ) == visited.end() )
                        queue.push_back( it2 );
                }
            }
        }

        cout << "Size of join_order_group in star groups ... " << join_order_group.size() << "\n";
        // size_t gn = group_number[ 1327 ];
        // cout << "Verification ... " << gn << "\t" << relation_groups[ gn ].size() << "\n";
    }

    SKey findLargest( map< SKey, size_t > degree ){
        int max_degree = 0, max_deg_rel = -1;
        for( auto it : degree ){
            if( it.second >= max_degree ){
                max_degree = it.second;
                max_deg_rel = it.first;
            }
        }
        return max_deg_rel;
    }

    size_t verifySize(){
        size_t total = 0;
        for( auto it : relation_groups ){
            total += it.second.size();
            checkj << it.first << ":\n";
            vector< SKey > group = it.second;
            for( auto r : group )
                checkj << "\t" << r << "\n";
            checkj << "\n";
        }

        return total;
    }

};

#endif
