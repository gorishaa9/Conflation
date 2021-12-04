#include "LeafEntry.hh"

LeafEntry::LeafEntry(boost::shared_ptr<Rectangle> mbr, boost::shared_ptr<HilbertValue> lhv, long long int id)
{
    this->mbr=mbr;
    this->lhv=lhv;
    this->id=id;
}

LeafEntry::~LeafEntry()
{

}

boost::shared_ptr<HilbertValue> LeafEntry::getLHV()
{
    assert(!!(this->lhv));
    return this->lhv;
}

boost::shared_ptr<Rectangle> LeafEntry::getMBR()
{
    assert(!!(this->mbr));
    return this->mbr;
}

long long int LeafEntry::getId(){
    return this->id;
}

bool LeafEntry::isLeafEntry()
{
    return true;
}
