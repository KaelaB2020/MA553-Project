#include <iostream>//Standard Library for input/output
#include <iomanip>//Standrd Library for I/0 manipulation (need setw() function)
#include <cstdlib>//Allows use of the rand() function
#include <ctime>//Allows querying numbers from the clock.
#include <string>
#include <fstream>//Needed to enable matrix saving.
#include <limits>
#include <initializer_list>
using namespace std;
//TO-DO:
//Create a class specifically for solving linear systems.
//Make an LU factorization function?
//Cholesky decomposition?
//Version 1.1:
//  - csvwrite method enables saving matrices as .csv files.
//Version 1.0:
//	- Matrix-matrix addition and multiplication.
//	- Matrix-scalar addition and multiplication.
//	- Constructors for memory preallocation and random 1-digit integer matrices.
//	- Method foor finding a matrix' eigenvalues, trace, and determinant.
class matrix{
//CLASS MEMBERS
public:
  int R;//Rows of the matrix (stack)
  int C;//Columns of the matrix (stack)
  int N;//Number of elements in the matrix (stack)
  double* DATA;//Matrix coefficients (heap allocated)

//--------------------------CONSTRUCTORS-------------------------------
public:
  //Default Constructor
  matrix(){
    this->R=1;//Default number of rows is 1.
    this->C=1;//Default number of columns is 1.
    this->N=1;//Default number of elements is 1.
    this->DATA = new double[1];//Default data is a 1x1 0-matrix
  }//end constructor

  //Constructor for simple memory allocation
  matrix(int r, int c){
    this->R=r;
    this->C=c;
    this->N=r*c;
    this->DATA = new double[N];
  }//end constructor

  //Constructor for predefined values (i.e., MATLAB's zeros and ones)
  matrix(int r, int c, double val){
    this->R=r;
    this->C=c;
    this->N=r*c;
    this->DATA=new double[N];
    for(int i=0; i<N;i++){
      DATA[i]=val;//Assign the same value to all matrix entries.
    }//end "i" loop
  }//end constructor

  //Constructor for special matrices
  matrix(int r, int c, string a){
    this->R=r;
    this->C=c;
    this->N=r*c;
    this->DATA=new double[N];

    //Random integer matrix
    if(a=="randi"){
      srand(time(NULL));//Use clock time to generate random seeds.
      for(int i=0; i<N;i++){
        DATA[i] = rand() % 10;//Assign pseudo-random integers (0-9) to vector B
      }//end "i" loop
    }//end if

    //Identity matrix
    else if(a=="eye"){
      int pivot;
      for(int i=0;i<R;i++){
        pivot = i*(1+C);//Index of the pivot element
        *(DATA+pivot) = 1.0;//Set pivot elements equal to 1.
        for(int j =pivot+1;j<pivot+C;j++){//All elements until next pivot...
          *(DATA+j) = 0.0;//Set them to zero
        }//end "j" loop
      }//end "i" loop
    }//end else if
  }//end constructor

  //Constructor with initializer lists
  matrix(std::initializer_list<std::initializer_list<double>> init_list){
    int k=0;//Counting variable
    int size1=0;
    int size2=0;
    //ensure that each sublist is of the same size.
    for(auto row : init_list){
      size1 = row.size();//Get number of elements of sublist.
      if(k<1){
        size2=size1;//Make a copy of size1 for comparisson purposes.
        k+=1;//Update k so that no further copying happens.
      }//end if
      else if(size1 != size2){
        throw "ERROR: Rows with uneven elements input!";
      }
    }//end "row" loop
    k=0;//reset counter.
    this->R=init_list.size();//Assign the number of rows.
    this->C=size1;//Assign the number of columns.
    this->N=init_list.size()*size1;//Compute number of element entries.
    this->DATA=new double[N];//Allocate memory for the matrix coefficients.
    for(auto row : init_list){//re-scan the sublists
      for(auto col : row){//scan the elements inside the sublists.
        *(DATA+k) = col;//assign them to the DATA array contiguously.
        k+=1;//Update counter
      }//end "row" loop
    }//end "col" loop
  };//end constructor

  //Class destructor:
  ~matrix(){
    delete &(this->DATA);
    delete &(this->R);
    delete &(this->C);
    delete &(this->N);
  };//end destructor

//--------------------------METHODS-------------------------------
public:
  //Make a copy of a matrix
  matrix copymat(matrix A){
    matrix copy(A.R,A.C);//Declare a
    for(int i=0;i<A.N;i++){
      copy.DATA[i] = A.DATA[i];
    }//end "i" loop
    return copy;
  }

  //Compute the trace of a matrix
  double trace(void){
    if(R==C){
      double sum = 0;//Running sum variable for the trace.
      for(int i=0;i<C;i++){
        sum += DATA[i*(C+1)];//Reference diagonal elements
      }//end "i" loop
      return sum;
    }//end if
    else{
      cout<<"WARNING: Trace undefined for non-square matrix!"<<endl;
      return std::numeric_limits<double>::quiet_NaN();
    }//end else
  }//end trace method

  //Compute the eigenvalues of a square matrix.
  double eigen(void){
    if(R==C){
      cout<<"WIP"<<endl;
    }//end if
    else{
      cout<<"WARNING: Eigenvalues undefined for non-square matrix!"<<endl;
      return std::numeric_limits<double>::quiet_NaN();
    }//end else
  }//end eigenvalue method.

  //Naive forward elimination
  void fel(void){
    int swap = 0;//Counter variable for number of partial pivots.
    for(int i=0; i<(R-1); i++){//Scan all rows expcet the very last
      for(int k=i+1; k<R; k++){
        //Conditional Partial pivoting.
        if(*(DATA+i*C+i) == 0){
          double temp[C];//temporary vector to store row with pivot  = 0;
          for(int l=0; l<R;l++){
            temp[l] = *(DATA+i*C+l); //Copy elements of 0-pivot row to temporary vector
          }//end "l" loop
          for(int l=0; l<R;l++){
            *(DATA+i*C+l) = *(DATA+N-C+l);//Overwrite current row with a copy of the last row.
          }//end "l" loop
          for(int l=0; l<R;l++){//end constructor
            *(DATA+N-C+l) = temp[l];//Overwrite last row with data stored in the temporary vector.
          }//end "l" loop
          swap += 1;
        }//end if statement
        double ratio = *(DATA+k*C+i)/(*(DATA+i*C+i));
        for(int m=1; m<C; m++){//loop for eliminating the remainder of the row.
          *(DATA+k*C+m) = *(DATA+k*C+m) - ratio*(*(DATA+i*C+m));
          *(DATA+k*C+i) = 0.0;
        }//end "m" loop
      }//end "k" loop
    }//end "i" loop
  }//end naive forward elimination function.

  //Return number of elements in matrix
  int nel(void){
    return N;
  }//end nel() method.

  //Return number of rows
  int cols(void){
    return C;
  }//end cols() method

  //Return number of columns
  int rows(void){
    return R;
  }//end rows() method.

  //Print matrix to console
  void display(void){
    for(int i=0; i<R; i++){//Scan the rows.
      cout <<"[";//Announce beginning of row by using a square bracket.
      for(int j=0; j<C; j++){//scan the columns.
        cout<< setw(10);//This is a reasonable width for outputting matrix elements.
        cout<< *(DATA+i*C+j) <<' ';//Output values of the matrix without jumping lines.
      }//end "j" for loop
      cout<<"]"<<endl;//Jump to the next line.
    }//end "i" for loop
  }//end display() method

  //Write matrix contents to .csv file
  void csvwrite(std::string name){
  ofstream outfile;//Ouput stream file object.
  outfile.open(name);//Initilize the output stream.
  for(int i=0;i<R;i++){//Scan rows of matrix to save.
    for(int j=0;j<C;j++){//Scan columns of matrix to save.
      outfile<<*(DATA+i*C+j)<<",";//Write column entries to the stream.
    }//end "j" loop.
    outfile<<std::endl;//jump line in the stream.
  }//end "i" loop
  outfile.close();//Output file
  }//end csvwrite.

void csvread(std::string name){
  ifstream infile;//Input stream file object.
  char data;
  infile.open(name);
  infile>>data;
  std::cout<<data<<std::endl;
  infile.close();
}//end csvread


//--------------------------OPERATOR OVERLOADING-------------------------------
  //Single index matrix element accessor:
  double &operator()(int i){
    if(i>=N){
      throw "ERROR: Index out of bounds!";
    }//end if
    return *(DATA +i);
  }//end single index accessor

  //Subscript (double index) element accessor:
  double &operator()(int r,int c){
    if(r > (R-1)||r<0){
      throw "ERROR: Row index out of bounds!";
    }
    if(c > (C-1)||c<0){
      throw "ERROR: Column index out of bounds!";
    }
    return *(DATA + r*C+c);
  }//end double index accessor

  //Slice matrix element accessor
  matrix &operator()(const matrix& B){
    if(B.R > (this->R)){
      throw "ERROR: Index matrix rows exceed those of reference matrix!";
    }
    if(B.C > (this->C)){
      throw "ERROR: Index matrix columns exceed those of reference matrix!";
    }
    //if(*std::max_element(B.DATA,B.DATA+B.N) >= (this->N)){
    //  throw "ERROR: Index matrix has at least one index out of bounds!";
    //}
    matrix C(B.R,B.C);//Instantiate the subset matrix.
    for(int i=0;i<B.R;i++){//Scan rows of index matrix.
      for(int j=0;j<B.C;j++){//Scan columns of index matrix.
        C.DATA[(i*B.C)+j] = *(DATA+int(B.DATA[(i*B.C)+j]));
      }//end "j" loop
    }//end "i" loop
    return C;
  }//end matrix slice accessor

  //Copy a matrix to another variable declared as a matrix. (matA=matB)
  void operator=(matrix B){
    this->R = B.R;
    this->C = B.C;
    this->N = B.N;
    for(int i =0; i<N;i++){
      *(DATA + i) = B.DATA[i];
    }//end "i" loop
  }//end matA=matB

  //Matrix-constant addition (matA+B) output to new matrix instance.
  matrix operator+(double B){
    matrix ApB(this->R,this->C);
    for(int i=0;i<N;i++){
      ApB.DATA[i] = *(DATA+i) + B;
    }//end "i" loop
    return ApB;
  }//end matA + B (instantiate)

  //Matrix-constant addition (matA+B) update current matrix.
  void operator+=(double B){
    for(int i=0;i<N;i++){
      *(DATA+i) += B;
    }//end "i" loop
  }//end matA + B (update)

  //Matrix-constant multiplication (matA*B) output to new matrix instance.
  matrix operator*(double B){
    matrix AtB(this->R,this->C);
    for(int i=0;i<N;i++){
      AtB.DATA[i] = (*(DATA+i))*B;
    }//end "i" loop
    return AtB;
  }// matA*B (instantiate)

  //Matrix-constant multiplication (matA*B) update current matrix.
  void operator*=(double B){
    for(int i=0;i<N;i++){
      *(DATA+i) *= B;
    }//end "i" loop
  }//end matA + B (update)

  //Matrix-matrix addition (matA+matB)
  matrix operator+(const matrix& B){
    if (R == B.R && C == B.C){
      matrix ApB(B.R,B.C);
      for(int i=0; i<N; i++){
        ApB.DATA[i] = *(DATA+i) + B.DATA[i];//Add elements of both matrices together.
      }//end "i" lloop
      return ApB;
    }//end if
    else{
      cout<<"ERROR: matrices are of unequal size"<<endl;
      return matrix();
    }//end else
  }//end matA+matB

  //Matrix-matrix  multiplication (matA*matB).
  matrix operator*(const matrix& B){
    if(C==B.R){
      matrix AtB(R,B.C);//Allocate space for the output matrix
      double sum = 0.0;//Running sum variable for elementwise computation
      for(int i=0;i<R;i++){//For all rows of the left matrix
        for(int j=0;j<B.C;j++){//For all columns of the right matrix
          for(int k=0;k<B.R;k++){//For all elements of the jth column's elements
            sum += (*(DATA+i*C+k))*B.DATA[j+k*C];//
          }//end "k" loop
          AtB.DATA[j+i*C] = sum;//Assign element to the output matrix.
          sum = 0.0;//Reset running sum variable
        }//end "j" loop
      }//end "i" loop
      return AtB;
    }//end if
    else{
      cout<<"ERROR: matrices are of unequal size"<<endl;
      return matrix();
    }//end else
  }//end matA*matB


};//end matrix CLASS
