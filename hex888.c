#include <petsc.h>
static char help[] = "Creates the element stiffness matrix.\n"
"Option prefix = opt_.\n";
//./hex888 -ksp_view_solution :sol.txt -mat_view draw -draw_pause 1
int main(int argc, char **args){
  PetscErrorCode ierr;//Error tracking variable.
  Mat K;//Global Stifness matrix.
      //KE;//Element Stiffness matrix

  Vec U, //Global Displacement vector.
      F; //Force vector.

  KSP ksp; //Linear solver object.

  double t1,t2;//time keepers
  PetscInt nelx = 2,//#elements distributed along the x-direction.
           nely = 2,//#elements distributed along the y-direction.
           nelz = 1,//#elements distributed along the z-direction.
           nx,//#of nodes distributed along the x-direction.
           ny,//#of nodes distributed along the y-direction.
           nz,//#of nodes distributed along the z-direction.
           dof = 3,//Degrees of Freedom per node.
           size,//Size of the global stiffness matrix K.
           plane;//Elements in x-y plane.
  PetscScalar nu = 0.3,//Poisson's ratio.
            lame,
            //E = 1,//Young's modulus
            load = -1;//Force magnitude.

  PetscInt //index24[24] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23},
           edof[24];
  PetscInt p,q,r,o,//looping variables.
           n1, n2, n3, n4;//Linear index variables of top face corners of element_pqr.
  ierr = PetscInitialize(&argc,&args,NULL,help); if(ierr) return ierr;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"opt_","options for mesh","");CHKERRQ(ierr);
  ierr = PetscOptionsInt("-nelx","Number of x-elements","TopOpt.c",nelx,&nelx,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsInt("-nely","Number of y-elements","TopOpt.c",nely,&nely,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsInt("-nelz","Number of z-elements","TopOpt.c",nelz,&nelz,NULL); CHKERRQ(ierr);
  //ierr = PetscOptionsReal("-E","Young's Modulus","TopOpt.c",E,&E,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsReal("-nu","Poisson's Ratio","TopOpt.c",nu,&nu,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsReal("-load","Force applied.","TopOpt.c",load,&load,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

//***********************BEGIN (ELEMENT STIFFNESS MATRIX)***********************
lame = 1/((nu+1)*(1-2*nu));
  PetscReal k0  = (32-48*nu)*lame/144,//0
            k1  = (6)*lame/144,//1
            k2  = (-8)*lame/144,//2
            k3  = (6-24*nu)*lame/144,//3
            k4  = (-6+24*nu)*lame/144,//4
            k5  = (4)*lame/144,//5
            k6  = (3)*lame/144,//6
            k7  = (-6)*lame/144,//7
            k8  = (-10+12*nu)*lame/144,//8
            k9  = (3-12*nu)*lame/144,//9
            k10 = (-3)*lame/144,//10
            k11 = (-3+12*nu)*lame/144,//11
            k12 = (-4+12*nu)*lame/144,//12
            k13 = (-8+12*nu)*lame/144;//13

//Currently no better way (need to learn better syntax.)
PetscReal KEc[24][24] = {
{	k0 ,k1 ,k1 ,k2 ,k4 ,k4 ,k8 ,k7 ,k11,k5 ,k3 ,k6 ,k5 ,k6 ,k3 ,k8 ,k11,k7 ,k13,k10,k10,k12,k9 ,k9 },//0
{	k1 ,k0 ,k1 ,k3 ,k5 ,k6 ,k7 ,k8 ,k11,k4 ,k2 ,k4 ,k6 ,k5 ,k3 ,k9 ,k12,k9 ,k10,k13,k10,k11,k8 ,k7 },//1
{	k1 ,k1 ,k0 ,k3 ,k6 ,k5 ,k6 ,k9 ,k12,k6 ,k3 ,k5 ,k4 ,k4 ,k2 ,k7 ,k11,k8 ,k10,k10,k13,k11,k7 ,k8 },//2
{	k2 ,k3 ,k3 ,k0 ,k7 ,k7 ,k5 ,k4 ,k10,k8 ,k1 ,k9 ,k8 ,k9 ,k1 ,k5 ,k10,k4 ,k12,k11,k11,k13,k6 ,k6 },//3
{	k4 ,k5 ,k6 ,k7 ,k0 ,k1 ,k3 ,k2 ,k4 ,k1 ,k8 ,k11,k11,k12,k9 ,k10,k5 ,k3 ,k9 ,k8 ,k7 ,k6 ,k13,k10},//4
{	k4 ,k6 ,k5 ,k7 ,k1 ,k0 ,k10,k3 ,k5 ,k11,k9 ,k12,k1 ,k11,k8 ,k3 ,k4 ,k2 ,k9 ,k7 ,k8 ,k6 ,k10,k13},//5
{	k8 ,k7 ,k6 ,k5 ,k3 ,k10,k0 ,k1 ,k7 ,k2 ,k4 ,k3 ,k13,k10,k6 ,k12,k9 ,k11,k5 ,k6 ,k4 ,k8 ,k11,k1 },//6
{	k7 ,k8 ,k9 ,k4 ,k2 ,k3 ,k1 ,k0 ,k7 ,k3 ,k5 ,k10,k10,k13,k6 ,k11,k8 ,k1 ,k6 ,k5 ,k4 ,k9 ,k12,k11},//7
{	k11,k11,k12,k10,k4 ,k5 ,k7 ,k7 ,k0 ,k4 ,k10,k5 ,k6 ,k6 ,k13,k9 ,k1 ,k8 ,k3 ,k3 ,k2 ,k1 ,k9 ,k8 },//8
{	k5 ,k4 ,k6 ,k8 ,k1 ,k11,k2 ,k3 ,k4 ,k0 ,k7 ,k1 ,k12,k11,k9 ,k13,k6 ,k10,k8 ,k9 ,k7 ,k5 ,k10,k3 },//9
{	k3 ,k2 ,k3 ,k1 ,k8 ,k9 ,k4 ,k5 ,k10,k7 ,k0 ,k7 ,k9 ,k8 ,k1 ,k6 ,k13,k6 ,k11,k12,k11,k10,k5 ,k4 },//10
{	k6 ,k4 ,k5 ,k9 ,k11,k12,k3 ,k10,k5 ,k1 ,k7 ,k0 ,k11,k1 ,k8 ,k10,k6 ,k13,k7 ,k9 ,k8 ,k4 ,k3 ,k2 },//11
{	k5 ,k6 ,k4 ,k8 ,k11,k1 ,k13,k10,k6 ,k12,k9 ,k11,k0 ,k1 ,k7 ,k2 ,k4 ,k3 ,k8 ,k7 ,k6 ,k5 ,k3 ,k10},//12
{	k6 ,k5 ,k4 ,k9 ,k12,k11,k10,k13,k6 ,k11,k8 ,k1 ,k1 ,k0 ,k7 ,k3 ,k5 ,k10,k7 ,k8 ,k9 ,k4 ,k2 ,k3 },//13
{	k3 ,k3 ,k2 ,k1 ,k9 ,k8 ,k6 ,k6 ,k13,k9 ,k1 ,k8 ,k7 ,k7 ,k0 ,k4 ,k10,k5 ,k11,k11,k12,k10,k4 ,k5 },//14
{	k8 ,k9 ,k7 ,k5 ,k10,k3 ,k12,k11,k9 ,k13,k6 ,k10,k2 ,k3 ,k4 ,k0 ,k7 ,k1 ,k5 ,k4 ,k6 ,k8 ,k1 ,k11},//15
{	k11,k12,k11,k10,k5 ,k4 ,k9 ,k8 ,k1 ,k6 ,k13,k6 ,k4 ,k5 ,k10,k7 ,k0 ,k7 ,k3 ,k2 ,k3 ,k1 ,k8 ,k9 },//16
{	k7 ,k9 ,k8 ,k4 ,k3 ,k2 ,k11,k1 ,k8 ,k10,k6 ,k13,k3 ,k10,k5 ,k1 ,k7 ,k0 ,k6 ,k4 ,k5 ,k9 ,k11,k12},//17
{	k13,k10,k10,k12,k9 ,k9 ,k5 ,k6 ,k3 ,k8 ,k11,k7 ,k8 ,k7 ,k11,k5 ,k3 ,k6 ,k0 ,k1 ,k1 ,k2 ,k4 ,k4 },//18
{	k10,k13,k10,k11,k8 ,k7 ,k6 ,k5 ,k3 ,k9 ,k12,k9 ,k7 ,k8 ,k11,k4 ,k2 ,k4 ,k1 ,k0 ,k1 ,k3 ,k5 ,k6 },//19
{	k10,k10,k13,k11,k7 ,k8 ,k4 ,k4 ,k2 ,k7 ,k11,k8 ,k6 ,k9 ,k12,k6 ,k3 ,k5 ,k1 ,k1 ,k0 ,k3 ,k6 ,k5 },//20
{	k12,k11,k11,k13,k6 ,k6 ,k8 ,k9 ,k1 ,k5 ,k10,k4 ,k5 ,k4 ,k10,k8 ,k1 ,k9 ,k2 ,k3 ,k3 ,k0 ,k7 ,k7 },//21
{	k9 ,k8 ,k7 ,k6 ,k13,k10,k11,k12,k9 ,k10,k5 ,k3 ,k3 ,k2 ,k4 ,k1 ,k8 ,k11,k4 ,k5 ,k6 ,k7 ,k0 ,k1 },//22
{	k9 ,k7 ,k8 ,k6 ,k10,k13,k1 ,k11,k8 ,k3 ,k4 ,k2 ,k10,k3 ,k5 ,k11,k9 ,k12,k4 ,k6 ,k5 ,k7 ,k1 ,k0 } //23
};
//Local stiffness matrix [this is a PETSc matrix]
/*
ierr = MatCreate(PETSC_COMM_WORLD,&KE); CHKERRQ(ierr);
ierr = MatSetSizes(KE,PETSC_DECIDE,PETSC_DECIDE,24,24); CHKERRQ(ierr);
ierr = MatSetFromOptions(KE); CHKERRQ(ierr);
ierr = MatSetUp(KE); CHKERRQ(ierr);
for(o=0;o<24;o++){
  ierr = MatSetValues(KE,1,&o,24,index24,KEc[o],INSERT_VALUES); CHKERRQ(ierr);
}//end "i" loop.
ierr = MatAssemblyBegin(KE,MAT_FINAL_Aierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
SSEMBLY); CHKERRQ(ierr);
ierr = MatAssemblyEnd(KE,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
*/
//************************END (ELEMENT STIFFNESS MATRIX)************************

//***********************BEGIN (GLOBAL STIFFNESS MATRIX)***********************
t1 = MPI_Wtime();
nx = nelx + 1; //3
ny = nely + 1; //3
nz = nelz + 1; //2
size = nx*ny*nz*dof;//Size of global stifness matrix. //54
plane = nx*ny;//Stepper between planes (in stencil). //9
ierr = MatCreate(PETSC_COMM_WORLD,&K); CHKERRQ(ierr);
ierr = MatSetSizes(K,PETSC_DECIDE,PETSC_DECIDE,size,size); CHKERRQ(ierr);
ierr = MatSetFromOptions(K); CHKERRQ(ierr);
ierr = MatSetUp(K); CHKERRQ(ierr);
for(p=0;p<nelz;p++){//Scan planes (z_dir).
  for(q=0;q<nely;q++){//Scan rows of plane_p (y_dir).
    for(r=0;r<nelx;r++){//Scan columns of row_pq (x_dir).
      n1 = r + q*nx+ p*plane;//LI of Corner 1 in stencil (Top-Left).
      n2 = n1 + 1;//LI of Corner 2 in stencil (Top-Right).
      n3 = n1 + plane;//LI of Corner 3 in stencil (Bottom-Left).
      n4 = n3 + 1;//LI of Corner 4 in stencil (Bottom-Right).

      //Node 1.
      edof[0] = dof*n1;//X-displacement of node 1 of element_pqr
      edof[1] = edof[0] + 1;//Y-displacement of node 1 of element_pqr
      edof[2] = edof[1] + 1;//Z-displacement of node 1 of element_pqr

      //Node 2.
      edof[3]  = dof*n2;//X-displacement of node 2 of element_pqr
      edof[4]  = edof[3] + 1;//Y-displacement of node 2 of element_pqr
      edof[5]  = edof[4] + 1;//Z-displacement of node 2 of element_pqr

      //Node 3.
      edof[6]  = dof*n3;//X-displacement of node 3 of element_pqr
      edof[7]  = edof[6] + 1;//Y-displacement of node 3 of element_pqr
      edof[8]  = edof[7] + 1;//Z-displacement of node 3 of element_pqr

      //Node 4.
      edof[9]  = dof*n4;//X-displacement of node 4 of element_pqr
      edof[10] = edof[9] + 1;//Y-displacement of node 4 of element_pqr
      edof[11] = edof[10]+ 1;//Z-displacement of node 4 of element_pqr

      //Node 5 (below node 1).
      edof[12] = edof[0] + dof*nx;//X-displacement of node 5 of element_pqr
      edof[13] = edof[12]+ 1;//Y-displacement of node 5 of element_pqr
      edof[14] = edof[13]+ 1;//Z-displacement of node 5 of element_pqr

      //Node 6 (below node 2).
      edof[15] = edof[3] + dof*nx;//X-displacement of node 6 of element_pqr
      edof[16] = edof[15]+ 1;//Y-displacement of node 6 of element_pqr
      edof[17] = edof[16]+ 1;//Z-displacement of node 6 of element_pqr

      //Node 7 (below nodPETSC_VIEWER_DRAW_WORLDe 3).
      edof[18] = edof[6] + dof*nx;//X-displacement of node 7 of element_pqr
      edof[19] = edof[18]+ 1;//Y-displacement of node 7 of element_pqr
      edof[20] = edof[19]+ 1;//Z-displacement of node 7 of element_pqr

      //Node 8 (below node 4).
      edof[21] = edof[9] + dof*nx;//X-displacement of node 8 of element_pqr
      edof[22] = edof[21]+ 1;//Y-displacement of node 8 of element_pqr
      edof[23] = edof[22]+ 1;//Z-displacement of node 8 of element_pqr

      //Print LI (debugging purposes).
      //printf("%4d\t%4d\t%4d\t%4d\t%4d\t%4d\n",edof[0],edof[1],edof[2],edof[3],edof[4],edof[5]);
      //printf("%4d\t%4d\t%4d\t%4d\t%4d\t%4d\n",edof[6],edof[7],edof[8],edof[9],edof[10],edof[11]);
      //printf("%4d\t%4d\t%4d\t%4d\t%4d\t%4d\n",edof[12],edof[13],edof[14],edof[15],edof[16],edof[17]);
      //printf("%4d\t%4d\t%4d\t%4d\t%4d\t%4d\n",edof[18],edof[19],edof[20],edof[21],edof[22],edof[23]);

      for(o=0;o<24;o++){
        ierr = MatSetValues(K,1,&edof[o],24,edof,KEc[o],ADD_VALUES); CHKERRQ(ierr);
      }//end "o" loop.
    }//end "r" loop.
  }//end "q" loop.
}//end "p" loop.
t2 = MPI_Wtime();
printf("Global stiffness matrix (K) build time:%10.5f [s]\n",t2-t1);

//Enforce Dirichlet condition on X-displacements at symmetry plane.
ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
PetscInt rdex;//Row index.
PetscScalar diag = 1;
for(p=0;p<nz;p++){//Scan z-rows of nodes
  for(q=0;q<ny;q++){//Scan y-nodes of row z.
    rdex = dof*(nx*q + plane*p);//LI in K of clamped X-displacements.
    ierr = MatZeroRows(K,1,&rdex,diag,NULL,NULL); CHKERRQ(ierr);
  }//end "q" loop.
}//end "p" loop.

//Enforce Dirichlet condition on Z displacements at pinned boundary.
for(o=0;o<nz;o++){
  //rdex = dof*(nx*ny-1 + o*plane) +2;//LI in K of clamped Z-displacements.
  rdex = dof*(-1 + (1+o)*plane) +2;//LI in K of clamped Z-displacements.
  ierr = MatZeroRows(K,1,&rdex,diag,NULL,NULL); CHKERRQ(ierr);
}//end "o" loop.
MatView(K,PETSC_VIEWER_DRAW_WORLD);
//************************END (GLOBAL STIFFNESS MATRIX)************************


//***********************BEGIN (FORCE VECTOR)***********************
ierr = VecCreate(PETSC_COMM_WORLD,&F); CHKERRQ(ierr);
ierr = VecSetSizes(F,PETSC_DECIDE,size); CHKERRQ(ierr);
ierr = VecSetFromOptions(F); CHKERRQ(ierr);
for(o=0;o<nz;o++){
  PetscInt idx = dof*plane*o + 2;
  //printf("%4d\n",idx);
  ierr = VecSetValues(F,1,&idx,&load,INSERT_VALUES);
}//end "o" loop.
ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
//************************END (FORCE VECTOR)************************

//***********************BEGIN (DISPLACEMENT VECTOR)***********************
ierr = VecDuplicate(F,&U); CHKERRQ(ierr);
ierr = VecSet(U,0.0); CHKERRQ(ierr);
//***********************BEGIN (DISPLACEMENT VECTOR)***********************

t1 = MPI_Wtime();
//***********************BEGIN (SOLUTION)***********************
//Linear solver setup.
ierr = KSPCreate(PETSC_COMM_WORLD,&ksp); CHKERRQ(ierr);
ierr = KSPSetOperators(ksp,K,K); CHKERRQ(ierr);
ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
ierr = KSPSolve(ksp,F,U); CHKERRQ(ierr);
ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//***********************END (SOLUTION)***********************
t2 = MPI_Wtime();

printf("Solve KU = F time:%10.5f [s]\n",t2-t1);
PetscViewer mfile;
const char* name = "Displacement";
ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,name,FILE_MODE_WRITE,&mfile);
ierr = VecView(U,mfile);


//Cleanup
KSPDestroy(&ksp);
//MatDestroy(&KE);
MatDestroy(&K);
VecDestroy(&U); VecDestroy(&F);
return PetscFinalize();
}//end main program.
