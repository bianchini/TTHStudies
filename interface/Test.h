#ifndef TEST_H
#define TEST_H

#include <string>
#include <map>
#include <stdio.h>


class Test {

 public:

  Test();
  ~Test();

  float  maxPt2jet(float , float);
  float  maxPt3jet(float , float,  float);
  float  maxPt4jet(float , float,  float, float);
  float  maxPt5jet(float , float,  float, float, float);
  float  maxPt6jet(float , float,  float, float, float, float);


};

#endif
