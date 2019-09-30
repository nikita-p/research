#include <iostream>
#include "MC.C"

using namespace std;


void events(){
  MC a("../inputs/model/trees");
  a.Loop();
  //return;
}
