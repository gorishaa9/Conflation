#ifndef TUPLE_H
#define TUPLE_H

class Tuple{

private:
    boost::shared_ptr<Segment_creator> _seg;
    double _score;
    double _is_contained;

public:
    // Tuple(boost::shared_ptr<Segment_creator> seg,
    //       double score) :
    //       _seg(seg),
    //       _score(score),
    //       _is_contained(false) {
    //       }

    Tuple(boost::shared_ptr<Segment_creator> seg,
          double score,
          bool contained) :
          _seg(seg),
          _score(score),
          _is_contained(contained) {
          }

    boost::shared_ptr<Segment_creator> segment(){
        return _seg;
    }

    bool is_interpolated(){
        return _seg->get_id()<0;
    }

    double score(){
        return _score;
    }

    void set_score(double score){
        _score = score;
    }

    inline boost::shared_ptr<Terminal_node> get_first(){
        return _seg->get_first();
    }

    inline boost::shared_ptr<Terminal_node> get_second(){
        return _seg->get_second();
    }

    bool is_contained(){
        return _is_contained;
    }

    void set_contained(bool contained){
        _is_contained = contained;
    }
};

#endif
