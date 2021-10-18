#include "matrix.h"
#include <iostream>
#include <initializer_list>
int main(){
matrix C;
matrix A({{0,1,2},{3,4,5},{6,7,8}});
matrix B({{0,1},{2,3}});
A.display();
std::cout<<A(1,1)<<std::endl;
B = A(B);
B.display();
//std::initializer_list list={1,2,3};
//for(int i =0,i<3,i++){
//  std::cout<<list[i]<<std::endl;
//}
//A.csvread("A.csv");
return 0;
}
