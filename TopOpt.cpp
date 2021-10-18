#include "matrix.h"
int main(){
//Discretization parameters.
int nelx = 9;//Number of x-elements.
int nely = 9;//Number of y-elements.
int nx = nelx+1;//Number of nodes along x-direction.
int ny = nely+1;//Number of nodes along y-direction.
int size = nx*ny*2;
//Elasticity constants.
double E = 1;//Young's Modulus
double nu =0.3;//Poisson's ratio.

//Lame's constants for the elements.
matrix k({{(3-nu)/6},
          {(1+nu)/8},
          {-(3+nu)/12},
          {(3*nu-1)/8},
          {(nu-3)/12},
          {-(nu+1)/8},
          {nu/6},
          {(1-3*nu)/8}});
k*=(E/(1-nu*nu));
k.display();
//Local stiffness matrix
matrix KE({{k(0),k(1),k(2),k(3),k(4),k(5),k(6),k(7)},
           {k(1),k(0),k(7),k(6),k(5),k(4),k(3),k(2)},
           {k(2),k(7),k(0),k(5),k(6),k(3),k(4),k(1)},
           {k(3),k(6),k(5),k(0),k(7),k(2),k(1),k(4)},
           {k(4),k(5),k(6),k(7),k(0),k(1),k(2),k(3)},
           {k(5),k(4),k(3),k(2),k(1),k(0),k(7),k(6)},
           {k(6),k(3),k(4),k(1),k(2),k(7),k(0),k(5)},
           {k(7),k(2),k(1),k(4),k(3),k(6),k(5),k(0)}});//hardcoded for all elements
KE.display();

//Displacement vector.
matrix U(size,1);//Preallocate memory for the displacement vector.
U(size*size-1) = 0;//Clamp vertical displacement at BR corner.
U(size*size-2) = 0;//Clamp horizontal displacement at BR corner.

//Force vector.
matrix F(size,1);//Preallocate memory for the force vectors.
F(1) = -1;//Specify the load at the middle.
//Global matrix builder
matrix K(size,size)
int n1, n2;
for(int i=0;i<nely;i++){
  for(int j=0;j<nely;j++){
    n1 = ny*();
    n2 = ;
  }//end "j" loop
}//end "i" loop


}// end main
