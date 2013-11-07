#include "./Test.h"
#include "TMath.h"


Test::Test(){}

Test::~Test(){}



float Test::maxPt2jet(float ptH1, float ptH2){

  return TMath::Max( ptH1, ptH2);

}


float Test::maxPt3jet(float ptH1, float ptH2, float ptA1){

  return  TMath::Max( ptA1 , this->maxPt2jet(ptH1, ptH2) );

}


float Test::maxPt4jet(float ptH1, float ptH2, float ptA1, float ptA2){

  return  TMath::Max( ptA2 , this->maxPt3jet(ptH1, ptH2, ptA1) );

}

float Test::maxPt5jet(float ptH1, float ptH2, float ptA1, float ptA2, float ptA3){

  return  TMath::Max( ptA3 , this->maxPt4jet(ptH1, ptH2, ptA1, ptA2) );

}

float Test::maxPt6jet(float ptH1, float ptH2, float ptA1, float ptA2, float ptA3, float ptA4){

  return  TMath::Max( ptA4 , this->maxPt5jet(ptH1, ptH2, ptA1, ptA2, ptA3) );

}
