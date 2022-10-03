#ifndef __BINOMIAL_TREE_CPP
#define __BINOMIAL_TREE_CPP

#include "binomial_tree.hpp"


BinomialTreeBase::BinomialTreeBase(){}
BinomialTreeBase::BinomialTreeBase(int num_steps, Option *pOption) : _num_steps(num_steps), _pOption(pOption) {}
BinomialTreeBase::~BinomialTreeBase(){}

#endif 