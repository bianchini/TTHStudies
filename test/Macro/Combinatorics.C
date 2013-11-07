#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>

using namespace std; 


class Weight {

public:
  
  Weight(unsigned int,std::vector<float>*,std::vector<float>*,int);
  
  ~Weight();

  float probability();
  float atLeastOneT( unsigned int );
  float atLeastOneL( unsigned int );


 private:

  unsigned int nparticles_;
  unsigned int index_;
  int verbose_;
  std::vector<float>* probsL_;
  std::vector<float>* probsT_;

};



Weight::Weight(unsigned int nparticles, 
	       std::vector<float>* probsL,
	       std::vector<float>* probsT,
	       int verbose
	       ){
  
  nparticles_ = nparticles;
  index_      = -1;
  verbose_    = verbose;

  probsL_ = new vector<float>();
  probsT_ = new vector<float>();
  
  if( probsL->size() == probsT->size() && probsT->size() == nparticles){
    for(unsigned int k = 0; k < probsL->size(); k++){
      probsL_->push_back( (*probsL)[k] );
    }
    for(unsigned int k = 0; k < probsT->size(); k++){
      probsT_->push_back( (*probsT)[k] );
    }
  }
  
}


Weight::~Weight(){
  delete probsL_; 
  delete probsT_; 
}


float Weight::probability(){

  if( nparticles_< 2 ) return 0.;
  index_++;
  if(verbose_) cout << "probability(): Index " << index_ << endl;

  float prob = 0.;
  if( index_< nparticles_- 2 ){
    if(verbose_) cout << "Index " << index_ << ", I am in [1]" << endl;
    prob = ((*probsL_)[index_]-(*probsT_)[index_])* (this->Weight::atLeastOneT(index_+1)) + (*probsT_)[index_]*(this->Weight::atLeastOneL(index_+1)) +  (1-(*probsL_)[index_])*(this->Weight::probability());
  }
  else{
    prob = ((*probsL_)[nparticles_-2]-(*probsT_)[nparticles_-2]) * (*probsT_)[nparticles_-1] + (*probsT_)[nparticles_-2] * (*probsL_)[nparticles_-1];
    if(verbose_) cout << "Index " << index_ << ", I am in [2]" << endl;
  }

  return prob;

}

float Weight::atLeastOneT( unsigned int indexT ){

  if( indexT >= nparticles_) return 0.;

  if( indexT< nparticles_- 1 ){
    if(verbose_) cout << "indexT " << indexT  << ", I am in [3]" << endl;
    return ( (*probsT_)[indexT] + (1-(*probsT_)[indexT])* (this->Weight::atLeastOneT( indexT + 1 ) ) );
  }
  else return (*probsT_)[indexT];

}

float Weight::atLeastOneL( unsigned int indexL ){

  if( indexL >= nparticles_) return 0.;

  if( indexL< nparticles_- 1 ){
    if(verbose_) cout << "indexL " << indexL  << ", I am in [4]" << endl;
    return ( (*probsL_)[indexL] + (1-(*probsL_)[indexL])* (this->Weight::atLeastOneL( indexL + 1 ) ) );
  }
  else return (*probsL_)[indexL];

}


void myfunction(unsigned int nparticles = 3){

  std::vector<float>* probsL = new  std::vector<float>();
  std::vector<float>* probsT = new  std::vector<float>();

  for(unsigned int i = 0; i < nparticles; i++ ){
    probsL->push_back(0.10);
    probsT->push_back(0.02);
  }

  //probsL->push_back(0.20);
  //probsL->push_back(0.20);
  //probsL->push_back(0.20);
  //probsT->push_back(0.05);
  //probsT->push_back(0.05);
  //probsT->push_back(0.05);


  Weight* myWeight = new Weight( nparticles, probsL, probsT , 0);

  cout << "Starting now..." << endl;

  cout << myWeight->probability() << endl;

  delete myWeight;

}
