#include <cstdlib>

#include "Kmer.hpp"
#include "KmerIntPair.hpp"


KmerIntPair::KmerIntPair(const Kmer &km, unsigned int val) {
  SetKey(km);
  SetVal(val);
}

void KmerIntPair::SetVal(unsigned int val) {
  char val8 = (val > 0xFF) ?  0xFF : (char)val;
  //memcpy(&this->v + KmerIntPair::IntOffset, &val8, sizeof(uint8_t));
  this->v[KmerIntPair::IntOffset] = val8;
}


// use:  b = p.ParallelIncrement()
// pre:  0 <= p.GetVal() <= 255
// post: if b is true then p has been atomically incremented 
//       and p.GetValue() < 255 prior to the increment
//       if b is false, the value was 255
bool KmerIntPair::ParallelIncrement() {
  unsigned int val = this->v[KmerIntPair::IntOffset];
  unsigned int newval;
  do {
    val = this->v[KmerIntPair::IntOffset]; // read new value
    if (val == 255) {
      return false; // someone else beat us to it
    }
    newval = val + 1;
  } while (!__sync_bool_compare_and_swap(this->v+KmerIntPair::IntOffset,val, newval));
  return true;
}

unsigned int KmerIntPair::GetVal() const {
  //uint8_t tmp = *reinterpret_cast<const uint8_t*>(this+KmerIntPair::IntOffset);
  return (uint8_t)this->v[KmerIntPair::IntOffset];
}

const Kmer& KmerIntPair::GetKey() const {
  return *reinterpret_cast<const Kmer*>(this + KmerIntPair::KmerOffset);
}

void KmerIntPair::SetKey(const Kmer& km) {
  memcpy(this, &km, sizeof(Kmer));
}

void SetKmerKey::operator()(KmerIntPair *value, const Kmer& km) {
  memcpy(value + KmerIntPair::KmerOffset, &km, sizeof(Kmer));
}


