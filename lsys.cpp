#include<string>
#include <iomanip>
#include "matrix.h"
class lsys{
  //CLASS MEMBERS
  public:
    matrix A;//Coefficient matrix
    matrix RHS;//RHS vector(s)

  //CONSTRUCTORS
  public:

    lsys(matrix matA, matrix matB){
      if(matA.R != matA.C){//Check to see whether A matrix is square.
        std::cout<<"ERROR: Coefficient matrix is not square!"<<std::endl;
      }//end if
      else if(matA.C != matB.R){//Check whether columns of A match rows of B.
        std::cout<<"ERROR: Columns of A do not match rows of B!"<<std::endl;
      }//end else if
      else{
        this->A = matA;
        this->RHS = matB;
      }//end else
    }//end constructor

    //METHODS:

    void show(void){
      for(int i=0;i<A.R;i++){//For all rows of the linear system
        std::cout<<"[";
        for(int j=0;j<A.C;j++){
          std::cout<<setw(10)<<A.DATA[i*A.C+j]<<" ";
        }//end "j" loop
        std::cout<<"][x"<<i<<"][";
        for(int j=0;j<RHS.C;j++){
          std::cout<<setw(10)<<RHS.DATA[i*RHS.C+j];
        }//end "j" loop
        std::cout<<"]"<<std::endl;
      }//end "i" loop
    }

    matrix solA(std::string b){//Analytically solve
      felAB(A,RHS);//Forward eliminate Coefficient matrix and RHS.
      matrix x = bas(A,RHS);//produce answer using backward substitution.

      //if(b == "fpb"){//forward elimination + backward substitution
      //}//end if
      /*
      else if(b == "LU"){//LU factorization
        std::cout<<"WIP"<<std::endl;
      }//end else if
      else if(b =="QR"){//QR factorization
        std::cout<<"WIP"<<std::endl;
      }//end else if
      else if(b == "Cholesky"){
        std::cout<<"WIP"<<std::endl;
      }//end else if
      else{
        std::cout<<"ERROR: Invalid solution method specified!"<<std::endl;
      }//end else
      */
      return x;
    }//End solA

    void solI(std::string b){//Iteratively solve
      if(b=="Jacobi"){//Jacobi Iterations
        std::cout<<"WIP"<<std::endl;
      }//end if
      else if(b=="Gauss-Seidel"){//Gauss-Seidel
        std::cout<<"WIP"<<std::endl;
      }//end else if
      else{
        std::cout<<"ERROR: Invalid solution method specified!"<<std::endl;
      }
    }//end Isol

    //Backward substitution
    matrix bas(matrix matA, matrix matB){
      matrix x(matA.R,matB.C);
      int rA = matA.R;
      int cA = matA.C;
      x.DATA[rA-1] = matB.DATA[rA-1]/matA.DATA[rA*cA-1];//Assign the "last" element of the solution vector.
      double sum = 0.0;//Running sum variable (must reset at every row).
      double term = 0.0;//Contribution variable (adds on to sum).

      for(int i=1; i<rA; i++){//Scan all rows from bottom to top.
        sum = 0.0;//Reset running sum.
        for(int j=0; j<i; j++){//Scan all columns from right to left.
          term = matA.DATA[cA*rA-1-i*cA-j]*x.DATA[rA-1-j];//Compute new "contribution".
          sum = sum + term;//Update sum
        }//end "j" loop.
      x.DATA[rA-1-i] = (matB.DATA[rA-1-i]-sum)/(matA.DATA[rA*cA-1-cA*i-i]);//Compute the next element of x.
      }//end "i" loop.
      return x;//Output solution vector(s)
    }//end bas function

    //Naive forward elimination
    void felAB(matrix matA, matrix matB){
      int rA = matA.R;
      int cA = matA.C;
      int rB = matB.R;
      int cB = matB.C;
      double ratio;//Row elimination ratio
      int debug = 0;//Set to 1 to print matrix at each stage of elimination.
      int swap = 0;//Counter variable for number of partial pivots.
      for(int i=0; i<(matA.R-1); i++){//Scan all rows expcet the very last
        for(int k=i+1; k<rA; k++){

          //Conditional Partial pivoting.
          if(matA.DATA[i*cA+i] == 0){
            double tempA[cA];//temporary vector to store row with pivot  = 0;
            double tempB[cB];//temporary vector to store corresponding matB's rows.
            for(int l=0; l<cA;l++){
              tempA[l] = matA.DATA[i*cA+l]; //Copy elements of 0-pivot row to temporary vector
            }//end "l" loop
            for(int l=0;l<cB;l++){
              tempB[l] = matB.DATA[i*cB+l];//Do the same but for the RHS vector(s)
            }//end "l" loop

            for(int l=0; l<cA;l++){
              matA.DATA[i*cA+l] = matA.DATA[rA*cA-cA+l];//Overwrite current row with a copy of the last row.
            }//end "l" loop
            for(int l=0; l<cB;l++){
              matB.DATA[i*cB+l] = matB.DATA[cB*rB-cB+l];//Same on RHS vector(s)
            }//end "l" loop

            for(int l=0; l<cA;l++){
              matA.DATA[rA*cA-cA+l] = tempA[l];//Overwrite last row with data stored in the temporary vector.
            }//end "l" loop
            for(int l=0; l<cB;l++){
              matB.DATA[rB*cB-cB+l] = tempB[l];//Overwrite last row with data stored in the temporary vector.
            }//end "l" loop

            swap += 1;
            std::cout<<"NOTE: Partial pivoting deployed! (# swaps = "<<swap<<")"<<std::endl;//Notify user that partial pivoting was needed.
          }//end if statement

          ratio = matA.DATA[k*cA+i]/matA.DATA[i*cA+i];
          for(int m=1; m<cA; m++){//loop for eliminating the remainder of the row.
            matA.DATA[k*cA+m] = matA.DATA[k*cA+m] - ratio*matA.DATA[i*cA+m];
          }//end "m" loop
          for(int m=1;m<cB;m++){
            matB.DATA[k*cB+m] = matB.DATA[k*cB+m] - ratio*matB.DATA[i*cB+m];
          }//end "m" loop

          if(debug==1){
            matA.DATA[k*cA]=0.0;// (DEBUG): Replace eliminated element with 0.
            std::cout<<k+1<<"th row eliminated"<<std::endl;//(DEBUG): Label row that was eliminated.
            matA.display();//(DEBUG): output matrix at intermediate stages
            std::cout<<ratio<<std::endl;//(DEBUG): output elimination ratio
          }//end if statement
        }//end "k" loop
      }//end "i" loop
      double detA=1.0;//variable for the determinant (starts as a running product)
      for(int i=0; i<rA; i++){
        detA *= matA.DATA[i*cA+i];//Multiply diagonal elements.
      }//end "i" loop
      if(swap%2==1){//If an odd number of swaps
        detA *= -1;//Account for the sign change on the determinant
      }//end if
      std::cout<<"The determinant of A is: detA = "<<detA<<std::endl;
    }//end fel function

};//end of class

int main(){
  matrix A(3,3,"randi");//Create a random matrix
  matrix B(3,1,"randi");//Create a random RHS
  /*
  std::cout<<"matrix A is:"<<std::endl;
  A.display();
  std::cout<<"matrix B is:"<<std::endl;
  B.display();
  */
  lsys prob(A,B);//Instantiate the linear system
  prob.show();
  matrix x = prob.solA("fpb");
  std::cout<<"the solution vector is"<<std::endl;
  x.display();
  return 0;
}//end main()
