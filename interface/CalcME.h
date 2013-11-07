#ifndef CALCME_H
#define CALCME_H


#include "TLorentzVector.h"
#include "TMatrixT.h"

#include <string>
#include <map>
#include <limits>  
#include <utility> 

using namespace std;


class CalcME {

 public:
  
  CalcME(TLorentzVector,TLorentzVector,TLorentzVector,TLorentzVector,TLorentzVector, double, double);
  ~CalcME();

  pair<TMatrixD,TMatrixD> Gamma(int);
  pair<TMatrixD,TMatrixD> Slash(TLorentzVector);
  pair<TMatrixD,TMatrixD> Product4    (TMatrixD , TMatrixD , TMatrixD , TMatrixD, TMatrixD , TMatrixD , TMatrixD , TMatrixD);
  TMatrixD Product4Real(TMatrixD , TMatrixD , TMatrixD , TMatrixD);

  pair<TMatrixD,TMatrixD> A0s(int,int,int);
  pair<TMatrixD,TMatrixD> A0tu(int,int,TLorentzVector,TLorentzVector,int);
  double Trace(int, int , int , int );
  double meSquared();

  void debug();

 private:

  TLorentzVector q1_;
  TLorentzVector q2_;
  TLorentzVector t1_;
  TLorentzVector t2_;
  TLorentzVector h_;
  double Mt_;
  double Mt2_;
  double Mh_;
  double Mh2_;
};




CalcME::CalcME(TLorentzVector q1,TLorentzVector q2,TLorentzVector t1, TLorentzVector t2, TLorentzVector h, double Mt, double Mh){

  q1_ = q1;
  q2_ = q2;
  t1_ = t1;
  t2_ = t2;
  h_  = h;
  
  Mt_  = Mt;
  Mt2_ = Mt*Mt;
  Mh_  = Mh;
  Mh2_ = Mh*Mh;
}


CalcME::~CalcME(){
  cout << "Destructor" << endl;
}

pair<TMatrixD,TMatrixD> CalcME::Gamma(int index){

  TMatrixD outR(4,4);
  TMatrixD outI(4,4);
  
  const double elementsR0[16] = 
    { 1, 0,  0,  0,
      0, 1,  0,  0,
      0, 0, -1,  0,
      0, 0,  0, -1
    };
  const double elementsI0[16] = 
    { 0, 0,  0,  0,
      0, 0,  0,  0,
      0, 0,  0,  0,
      0, 0,  0,  0
    };

  const double elementsR1[16] = 
    { 0,  0,  0,  1,
      0,  0,  1,  0,
      0, -1,  0,  0,
      -1, 0,  0,  0 
    };
  const double elementsI1[16] = 
    { 0, 0,  0,  0,
      0, 0,  0,  0,
      0, 0,  0,  0,
      0, 0,  0,  0
    };

  const double elementsR2[16] = 
    { 0, 0,  0,  0,
      0, 0,  0,  0,
      0, 0,  0,  0,
      0, 0,  0,  0
    };
  const double elementsI2[16] = 
    { 0,  0,  0, -1,
      0,  0,  1,  0,
      0,  1,  0,  0,
      -1, 0,  0,  0
    };

  const double elementsR3[16] = 
    { 0, 0,  1,  0,
      0, 0,  0, -1,
     -1, 0,  0,  0,
      0, 1,  0,  0
    };
  const double elementsI3[16] = 
    { 0, 0,  0,  0,
      0, 0,  0,  0,
      0, 0,  0,  0,
      0, 0,  0,  0
    };

  const double elementsR4[16] = 
    { 1, 0,  0,  0,
      0, 1,  0,  0,
      0, 0,  1,  0,
      0, 0,  0,  1
    };
  const double elementsI4[16] = 
    { 0, 0,  0,  0,
      0, 0,  0,  0,
      0, 0,  0,  0,
      0, 0,  0,  0
    };

  const double elementsR5[16] = 
    { 1, 0,  0,  0,
      0,-1,  0,  0,
      0, 0, -1,  0,
      0, 0,  0, -1
    };
  const double elementsI5[16] = 
    { 1, 0,  0,  0,
      0,-1,  0,  0,
      0, 0, -1,  0,
      0, 0,  0, -1
    };

  switch(index){
  case 0:
    outR.SetMatrixArray( elementsR0, "F");
    outI.SetMatrixArray( elementsI0, "F");
    break;
  case 1:
    outR.SetMatrixArray( elementsR1, "F");
    outI.SetMatrixArray( elementsI1, "F");
    break;
  case 2:
    outR.SetMatrixArray( elementsR2, "F");
    outI.SetMatrixArray( elementsI2, "F");
    break;
  case 3:
    outR.SetMatrixArray( elementsR3, "F");
    outI.SetMatrixArray( elementsI3, "F");
    break;
  case 4:
    outR.SetMatrixArray( elementsR4, "F");
    outI.SetMatrixArray( elementsI4, "F");
    break;
  case 5:
    outR.SetMatrixArray( elementsR5, "F");
    outI.SetMatrixArray( elementsI5, "F");
    break;
  default:
    break;
  }

  return make_pair( outR, outI);
   
}


pair<TMatrixD,TMatrixD> CalcME::Slash(TLorentzVector v){

 
  TMatrixD Re( Gamma(0).first );
  Re          *= (v[3]);
  Re          += ( Gamma(1).first  *= (-v[0]) );
  Re          += ( Gamma(3).first  *= (-v[2]) );

  TMatrixD Im( Gamma(2).second );
  Im          *= (-v[1]);

  return make_pair(Re,Im);

}



pair<TMatrixD,TMatrixD> CalcME::A0s(int mu, int nu, int dagger){

  TMatrixD ReA( Gamma(4).second ); 
  ReA += Slash(t1_).first;
  ReA += Slash(h_).first;
  ReA += (Gamma(4).first *= Mt_);
  ReA *= (1./(2*t1_.Dot(h_) + Mh2_));

  //cout << "ImA : " << endl;

  TMatrixD ImA( Gamma(4).second ) ;
  ImA += Slash(t1_).second;
  ImA += Slash(h_).second;
  ImA *= (1./(2*t1_.Dot(h_) + Mh2_));

  //cout << "ReB : " << endl;

  TMatrixD ReB( Gamma(4).second  );
  ReB += (Slash(t2_).first *= (-1));
  ReB += (Slash(h_).first  *= (-1));
  ReB += (Gamma(4).first *= Mt_);
  ReB *= (1./(2*t2_.Dot(h_) + Mh2_));

  //cout << "ImB : " << endl;

  TMatrixD ImB( Gamma(4).second );
  ImB += (Slash(t2_).second *= (-1));
  ImB += (Slash(h_).second  *= (-1));
  ImB *= (1./(2*t2_.Dot(h_) + Mh2_));

  TMatrixD ReTot(Gamma(4).second) ; // null
  TMatrixD ImTot(Gamma(4).second) ; // null

  //cout << "Loop : " << endl;

  for(int alpha = 0; alpha<4; alpha++){

    //cout << "Alpha = " << alpha << endl;

    int alphaLV = alpha;
    int muLV    = mu;
    int nuLV    = nu;

    switch(alpha){
    case 0:
      alphaLV = 3;      
      break;
    case 1:
      alphaLV = 0;
      break;
    case 2:
      alphaLV = 1;
      break;
    case 3:
      alphaLV = 2;
      break;
    default:
      break;
    }

    switch(mu){
    case 0:
      muLV = 3;      
      break;
    case 1:
      muLV = 0;
      break;
    case 2:
      muLV = 1;
      break;
    case 3:
      muLV = 2;
      break;
    default:
      break;
    }

    switch(nu){
    case 0:
      nuLV = 3;      
      break;
    case 1:
      nuLV = 0;
      break;
    case 2:
      nuLV = 1;
      break;
    case 3:
      nuLV = 2;
      break;
    default:
      break;
    }

    //cout << "ReG,ImG : " << endl;
    double V = 
      (  q1_[alphaLV] -   q2_[alphaLV])*(Gamma(5).first)(mu,nu)    + 
      (  q1_[muLV]    + 2*q2_[muLV])   *(Gamma(5).first)(nu,alpha) -
      (2*q1_[nuLV]    +   q2_[nuLV])   *(Gamma(5).first)(mu,alpha) ;
    V /= (2*q1_.Dot(q2_));

    TMatrixD ReG( Gamma(alpha).first  * (Gamma(5).first)(alpha,alpha) );
    TMatrixD ImG( Gamma(alpha).second * (Gamma(5).first)(alpha,alpha) );

    //cout << "Fill ReTot... " << endl;

    TMatrixD ReTotTmp( Gamma(4).first ); 

    if(dagger==0)
      ReTotTmp.Mult( ReA, ReG );
    else
      ReTotTmp.Mult( ReG, ReA );
    ReTotTmp  *= V;
    ReTot    += ReTotTmp;
    if(dagger==0)
      ReTotTmp.Mult(ReG,ReB);
    else
      ReTotTmp.Mult(ReB,ReG);
    ReTotTmp  *= V;
    ReTot    += ReTotTmp;
    if(dagger==0)
      ReTotTmp.Mult(ImA,ImG);
    else
      ReTotTmp.Mult(ImG,ImA);
    ReTotTmp  *= V;
    ReTot    -= ReTotTmp; 
    if(dagger==0)
      ReTotTmp.Mult(ImG,ImB);
    else
      ReTotTmp.Mult(ImB,ImG);
    ReTotTmp  *= V;
    ReTot    -= ReTotTmp;

    //cout << "Fill ImTot... " << endl;

    TMatrixD ImTotTmp( Gamma(4).first ); 
    
    if(dagger==0)
      ImTotTmp.Mult( ReA, ImG);
    else
      ImTotTmp.Mult( ImG, ReA);
    ImTotTmp *= V;
    ImTot   +=  ImTotTmp;
    if(dagger==0)
      ImTotTmp.Mult(ImG,ReB);
    else
      ImTotTmp.Mult(ReB,ImG);
    ImTotTmp *= V;
    ImTot   +=  ImTotTmp;
    if(dagger==0)
      ImTotTmp.Mult(ImA,ReG);
    else
      ImTotTmp.Mult(ReG,ImA);
    ImTotTmp *= V;
    ImTot   +=  ImTotTmp ;
    if(dagger==0)
      ImTotTmp.Mult(ReG,ImB);
    else
      ImTotTmp.Mult(ImB,ReG);
    ImTotTmp *= V;
    ImTot   +=  ImTotTmp;

  }

  return make_pair(ReTot,ImTot);

}



pair<TMatrixD,TMatrixD> CalcME::Product4(TMatrixD ReA, TMatrixD ReB, TMatrixD ReC, TMatrixD ReD,
					 TMatrixD ImA, TMatrixD ImB, TMatrixD ImC, TMatrixD ImD){

  TMatrixD m1 ( Product4Real(ReA,ReB,ReC,ReD) ); // Re Re Re Re

  TMatrixD m2 ( Product4Real(ImA,ImB,ReC,ReD) ); // Im Im Re Re
  TMatrixD m3 ( Product4Real(ImA,ReB,ImC,ReD) );
  TMatrixD m4 ( Product4Real(ImA,ReB,ReC,ImD) );
  TMatrixD m5 ( Product4Real(ReA,ImB,ImC,ReD) );
  TMatrixD m6 ( Product4Real(ReA,ImB,ReC,ImD) );
  TMatrixD m7 ( Product4Real(ReA,ReB,ImC,ImD) );

  TMatrixD m8 ( Product4Real(ImA,ImB,ImC,ImD) ); // Im Im Im Im

  TMatrixD m9 ( Product4Real(ImA,ReB,ReC,ReD) ); // Im Re Re Re
  TMatrixD m10( Product4Real(ReA,ImB,ReC,ReD) );
  TMatrixD m11( Product4Real(ReA,ReB,ImC,ReD) );
  TMatrixD m12( Product4Real(ReA,ReB,ReC,ImD) );

  TMatrixD m13( Product4Real(ImA,ImB,ImC,ReD) ); // Im Im Im Re
  TMatrixD m14( Product4Real(ImA,ImB,ReC,ImD) );
  TMatrixD m15( Product4Real(ImA,ReB,ImC,ImD) );
  TMatrixD m16( Product4Real(ReA,ImB,ImC,ImD) );

  TMatrixD Re( m1 );
  Re -= m2;
  Re -= m3;
  Re -= m4;
  Re -= m5;
  Re -= m6;
  Re -= m7;
  Re += m8;

  TMatrixD Im( m9 );
  Im += m10;  
  Im += m11;  
  Im += m12;  
  Im -= m13;  
  Im -= m14;  
  Im -= m15;  
  Im -= m16;

  return make_pair( Re, Im);  
}


TMatrixD CalcME::Product4Real(TMatrixD A, TMatrixD B, TMatrixD C, TMatrixD D){

  TMatrixD tmp1( A );
  TMatrixD tmp2( A );
  TMatrixD tmp3( A );

  tmp1.Mult( A, B);
  tmp2.Mult( C, D);
  tmp3.Mult( tmp1, tmp2);

  return tmp3;

}


pair<TMatrixD,TMatrixD> CalcME::A0tu(int mu, int nu, TLorentzVector q1, TLorentzVector q2, int dagger){

  TMatrixD ReA( Gamma(4).second ); 
  ReA += Slash(t1_).first;
  ReA += Slash(h_).first;
  ReA += (Gamma(4).first *= Mt_);
  TMatrixD ImA( Gamma(4).second ) ;
  ImA += Slash(t1_).second;
  ImA += Slash(h_).second;

  TMatrixD ReB( Gamma(mu).first ); 
  TMatrixD ImB( Gamma(mu).second ); 

  TMatrixD ReC( Gamma(4).second ); 
  ReC += Slash(q2).first;
  ReC -= Slash(t2_).first;
  ReC += (Gamma(4).first *= Mt_);
  TMatrixD ImC( Gamma(4).second ); 
  ImC += Slash(q2).second;
  ImC -= Slash(t2_).second;

  TMatrixD ReD( Gamma(nu).first ); 
  TMatrixD ImD( Gamma(nu).second ); 

  TMatrixD ReE( Gamma(4).second ); 
  ReE +=  Slash(t1_).first;
  ReE -=  Slash(q1).first;
  ReE += (Gamma(4).first *= Mt_);
  TMatrixD ImE( Gamma(4).second ); 
  ImE +=  Slash(t1_).second;
  ImE -=  Slash(q1).second;

  TMatrixD ReF( Gamma(4).second ); 
  ReF -=  Slash(t2_).first;
  ReF -=  Slash(h_).first;
  ReF += (Gamma(4).first *= Mt_);
  TMatrixD ImF( Gamma(4).second ); 
  ImF -=  Slash(t2_).second;
  ImF -=  Slash(h_).second;


  TMatrixD Re1( Gamma(4).second  );
  TMatrixD Im1( Gamma(4).second  );
  if(dagger==0){
    Re1 = Product4(ReA,ReB,ReC,ReD, ImA, ImB, ImC,ImD).first ;
    Im1 = Product4(ReA,ReB,ReC,ReD, ImA, ImB, ImC,ImD).second;
  }
  else{
    Re1 = Product4(ReD, ReC, ReB, ReA, ImD, ImC, ImB, ImA).first ;
    Im1 = Product4(ReD, ReC, ReB, ReA, ImD, ImC, ImB, ImA).second;
  }


  TMatrixD Re2( Gamma(4).second  );
  TMatrixD Im2( Gamma(4).second  );
  if(dagger==0){
    Re2 = Product4(ReB,ReE,ReC,ReD, ImB, ImE, ImC,ImD).first ;
    Im2 = Product4(ReB,ReE,ReC,ReD, ImB, ImE, ImC,ImD).second ;
  }
  else{
    Re2 = Product4( ReD, ReC, ReE, ReB, ImD, ImC, ImE, ImB).first ;
    Im2 = Product4( ReD, ReC, ReE, ReB, ImD, ImC, ImE, ImB).second ;
  }

  TMatrixD Re3( Gamma(4).second  );
  TMatrixD Im3( Gamma(4).second  );
  if(dagger==0){
    Re3 = Product4(ReB,ReE,ReD,ReF, ImB, ImE, ImD,ImF).first ;
    Im3 = Product4(ReB,ReE,ReD,ReF, ImB, ImE, ImD,ImF).second ;
  }
  else{
    Re3 = Product4( ReF, ReD, ReE, ReB , ImF, ImD, ImE, ImB).first;
    Im3 = Product4( ReF, ReD, ReE, ReB , ImF, ImD, ImE, ImB).second;
  }



  Re1 *= (1./( 2*t1_.Dot(h_) + Mh2_)/(-2*q2.Dot(t2_)));
  Im1 *= (1./( 2*t1_.Dot(h_) + Mh2_)/(-2*q2.Dot(t2_)));


  Re2 *= (1./(-2*q1.Dot(t1_))/(-2*q2.Dot(t2_)));
  Im2 *= (1./(-2*q1.Dot(t1_))/(-2*q2.Dot(t2_)));


  Re3 *= (1./(-2*q1.Dot(t1_))/(2*t2_.Dot(h_) + Mh2_));
  Im3 *= (1./(-2*q1.Dot(t1_))/(2*t2_.Dot(h_) + Mh2_));

  Re1 += Re2;
  Re1 += Re3;

  Im1 += Im2;
  Im1 += Im3;

  return make_pair(Re1, Im1);

}


double CalcME::Trace(int mu, int nu, int alpha, int beta){

  TMatrixD ReA( Gamma(4).second ); 
  ReA += A0s(mu,nu,0).first;
  TMatrixD ImA( Gamma(4).second ); 
  ImA += A0s(mu,nu,0).second;

  TMatrixD ReB( Gamma(4).second ); 
  ReB += Slash(t2_).first;
  ReB += (Gamma(4).first *= Mt_);
  TMatrixD ImB( Gamma(4).second ); 
  ImB += Slash(t2_).second;

  TMatrixD ReC( Gamma(4).second ); 
  ReC += A0s (alpha,beta,1).first;
  ReC += A0tu(alpha,beta,q1_,q2_,1).first;
  ReC -= A0tu(beta,alpha,q2_,q1_,1).first;
  TMatrixD ImC( Gamma(4).second ); 
  ImC += A0s (alpha,beta,1).second;
  ImC += A0tu(alpha,beta,q1_,q2_,1).second;
  ImC -= A0tu(beta,alpha,q2_,q1_,1).second;

  TMatrixD ReD( Gamma(4).second ); 
  ReD += Slash(t1_).first;
  ReD += (Gamma(4).first *= Mt_);
  TMatrixD ImD( Gamma(4).second ); 
  ImD += Slash(t1_).second;

  TMatrixD ReE( Gamma(4).second ); 
  ReE += A0tu(mu,nu,q1_,q2_,0).first;
  ReE += A0tu(nu,mu,q2_,q1_,0).first;
  TMatrixD ImE( Gamma(4).second ); 
  ImE += A0tu(mu,nu,q1_,q2_,0).second;
  ImE += A0tu(nu,mu,q2_,q1_,0).second;

  TMatrixD ReF( Gamma(4).second ); 
  ReF += A0tu(alpha,beta,q1_,q2_,1).first;
  ReF += A0tu(beta,alpha,q2_,q1_,1).first;
  TMatrixD ImF( Gamma(4).second ); 
  ReF += A0tu(alpha,beta,q1_,q2_,1).second;
  ReF += A0tu(beta,alpha,q2_,q1_,1).second;

  TMatrixD ReG( Gamma(4).second ); 
  ReG += A0tu(mu,nu,q1_,q2_,0).first;
  TMatrixD ImG( Gamma(4).second ); 
  ImG += A0tu(mu,nu,q1_,q2_,0).second;

  TMatrixD ReH( Gamma(4).second ); 
  ReH += A0tu(beta,alpha,q2_,q1_,1).first;
  TMatrixD ImH( Gamma(4).second ); 
  ImH += A0tu(beta,alpha,q2_,q1_,1).second;




  TMatrixD Re1( Gamma(4).second ); 
  Re1 += Product4(ReA,ReB,ReC,ReD, ImA, ImB, ImC,ImD).first ; //second
  Re1.Mult( ReB, ReD);

  TMatrixD Re2( Gamma(4).second ); 
  Re2 += Product4(ReE,ReB,ReF,ReD, ImE, ImB, ImF,ImD).first ; //second

  TMatrixD Re3( Gamma(4).second ); 
  Re3 += Product4(ReG,ReB,ReH,ReD, ImG, ImB, ImH,ImD).first ; //second

  double trace1 = Re1(0,0)+Re1(1,1)+Re1(2,2)+Re1(3,3);
  //double trace2 = Re2(0,0)+Re2(1,1)+Re2(2,2)+Re2(3,3);
  //double trace3 = Re3(0,0)+Re3(1,1)+Re3(2,2)+Re3(3,3);


  /////////////////////////////////////////////// TEST ////////////////////////////////////

  TMatrixD ReDeb( Gamma(4).second ); 
  TMatrixD ReDebTmp( Gamma(4).first ); 
  ReB.Print("");
  ReD.Print("");
  ImB.Print("");
  ImD.Print("");
  ReDebTmp.Mult(ReB,ReD);
  ReDeb += ReDebTmp;
  ReDebTmp.Mult(ImB,ImD);
  ReDeb -= ReDebTmp;
  ReDeb.Print("");

  cout << "Trace of the Re part: " <<  ReDeb(0,0)+ReDeb(1,1)+ReDeb(2,2)+ReDeb(3,3) << endl;
  trace1 = ReDeb(0,0)+ReDeb(1,1)+ReDeb(2,2)+ReDeb(3,3);

  TMatrixD ImDeb( Gamma(4).second ); 
  TMatrixD ImDebTmp( Gamma(4).first ); 
  ImDebTmp.Mult(ReB,ImD);
  ImDeb += ImDebTmp;
  ImDebTmp.Mult(ImB,ReD);
  ImDeb += ImDebTmp;
  ImDeb.Print("");

  cout << "Trace of the Im part: " << (ImDeb(0,0)+ImDeb(1,1)+ImDeb(2,2)+ImDeb(3,3)) << endl;
  trace1 += (ImDeb(0,0)+ImDeb(1,1)+ImDeb(2,2)+ImDeb(3,3));
  ///////////////////////////////////////////////

  return //(12*trace1 + 2/3.*trace2 - 13*trace3);
    trace1;

}

double CalcME::meSquared(){


  return Trace(0,0,0,0);

  double me = 0.0;

  for(int mu = 0; mu<4; mu++){
    for(int nu = 0; nu<4; nu++){
      for(int alpha = 0; alpha<4; alpha++){
	for(int beta = 0; beta<4; beta++){

	  if(mu!=alpha || nu!=beta) continue;

	  me += ( (Gamma(5).first)(mu,alpha) * (Gamma(5).first)(nu,beta) * Trace(mu,nu,beta,alpha)  );  

	}
      }
    }
  }

  return me;
}



void CalcME::debug(){

  cout << "Element (0,0):" << endl;
  cout << "Real part: " << endl;
  (A0s(0,0,0).first).Print("");
  cout << "Imaginary part: " << endl;
  (A0s(0,0,0).second).Print("");
  
  cout << "Element (1,1):" << endl;
  cout << "Real part: " << endl;
  (A0s(1,1,0).first).Print("");
  cout << "Imaginary part: " << endl;
  (A0s(1,1,0).second).Print("");

  cout << "Element (2,2):" << endl;
  cout << "Real part: " << endl;
  (A0s(2,2,0).first).Print("");
  cout << "Imaginary part: " << endl;
  (A0s(2,2,0).second).Print("");

  cout << "Element (3,3):" << endl;
  cout << "Real part: " << endl;
  (A0s(3,3,0).first).Print("");
  cout << "Imaginary part: " << endl;
  (A0s(3,3,0).second).Print("");
}

//CalcME::


#endif
