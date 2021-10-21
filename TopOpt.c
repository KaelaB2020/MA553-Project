
static char help[] = "PETSc version of O.Sigmund's 1999 MATLAB code.\n"
"Option prefix = opt_.\n";
#include <petsc.h>
//./TopOpt -ksp_view_solution -mat_view draw -draw_pause 100

int main(int argc, char **args){
  //Variable declarations
  PetscErrorCode ierr;
  //PetscViewer viewer;
  Vec F, U;//Force vector (F) and displacement vector (U).
  Mat K;//Global stiffness matrix.
  KSP ksp;//Linear solver object.

  //Grid parameters.
  PetscInt nelx = 2;//Number of elements along x-direction.
  PetscInt nely = 2;//Number of elements along y-direction.

  //Grid variables.
  PetscInt i, j, w, l;//Looping variables
  PetscInt n1,n2;//Storage variables for the LI of the TL and TR nodes of element_ji;
  PetscInt edof[8];//Array to store LI of DOF's for the global stiffness matrix.
  PetscInt LIFy[1] = {1};//LI of the dof where Fy is to act.

  //Elasticity constants
  PetscReal E = 1;//Young's modulus.
  PetscReal nu = 0.3;//Poisson's ratio

  PetscReal Fy=-1;//Force applied on the beam.

  ierr = PetscInitialize(&argc,&args,NULL,help); if(ierr) return ierr;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"opt_","options for mesh","");CHKERRQ(ierr);
  ierr = PetscOptionsInt("-nelx","Number of x-elements","TopOpt.c",nelx,&nelx,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsInt("-nely","Number of y-elements","TopOpt.c",nely,&nely,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsReal("-E","Young's Modulus","TopOpt.c",E,&E,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsReal("-nu","Poisson's Ratio","TopOpt.c",nu,&nu,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsReal("-Fy","Force applied.","TopOpt.c",Fy,&Fy,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  PetscInt nx = nelx+1;//Number of X-nodes.
  PetscInt ny = nely+1;//Number of Y-nodes.
  PetscInt size = 2*nx*ny; //Size of the global stiffness matrix.

  //Coefficients of the local stiffness matrix.
  PetscReal mult = E/(1-nu*nu);//A common factor shared by the element stifnesses.
  PetscReal k[8] = {+(3-nu)*mult/6,
                    +(1+nu)*mult/8,
                    -(3+nu)*mult/12,
                    +(3*nu-1)*mult/8,
                    +(nu-3)*mult/12,
                    -(nu+1)*mult/8,
                    +nu*mult/6,
                    +(1-3*nu)*mult/8};

  //Local stiffness matrix (same for all elements).
  PetscReal KE[8][8] = {{k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7]},
                        {k[1],k[0],k[7],k[6],k[5],k[4],k[3],k[2]},
                        {k[2],k[7],k[0],k[5],k[6],k[3],k[4],k[1]},
                        {k[3],k[6],k[5],k[0],k[7],k[2],k[1],k[4]},
                        {k[4],k[5],k[6],k[7],k[0],k[1],k[2],k[3]},
                        {k[5],k[4],k[3],k[2],k[1],k[0],k[7],k[6]},
                        {k[6],k[3],k[4],k[1],k[2],k[7],k[0],k[5]},
                        {k[7],k[2],k[1],k[4],k[3],k[6],k[5],k[0]}};

  //Force vector (Right-Hand-Side vector RHS).
  ierr = VecCreate(PETSC_COMM_WORLD,&F); CHKERRQ(ierr);
  ierr = VecSetSizes(F,PETSC_DECIDE,size); CHKERRQ(ierr);
  ierr = VecSetFromOptions(F); CHKERRQ(ierr);
  ierr = VecSetValues(F,1,LIFy,&Fy,INSERT_VALUES); CHKERRQ(ierr);//Apply the point force on the beam.
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

  //Global stiffness matrix
  ierr = MatCreate(PETSC_COMM_WORLD,&K); CHKERRQ(ierr);
  ierr = MatSetSizes(K,PETSC_DECIDE,PETSC_DECIDE,size,size); CHKERRQ(ierr);
  ierr = MatSetFromOptions(K); CHKERRQ(ierr);
  ierr = MatSetUp(K); CHKERRQ(ierr);

  //Building the global stiffness matrix.
  for(i=0;i<nely;i++){//scan rows of elements.
    for(j=0;j<nelx;j++){//scan columns of elements.

      //LI of nodes in the mesh
      n1 = j + i*ny;//LI of TL corner of element_ij.
      n2 = n1 + ny;
      /*
      //LI of the corresponding dof's of element_ij in displacement vector U.
      edof[0] = 2*n1;//LI of X-displacement of TL corner in U.
      edof[1] = edof[0] + 1;//LI of Y-displacement of TL corner in U.
      edof[2] = edof[1] + 1;//LI X-displacement of TR corner in U.
      edof[3] = edof[2] + 1;//LI Y-displacement of TR corner in U.
      edof[4] = edof[0] + 2*ny;//LI X-displacement of BL corner in U.
      edof[5] = edof[4] + 1;//LI Y-displacement of BL corner in U.
      edof[6] = edof[5] + 1;//LI X-displacement of BR corner in U.
      edof[7] = edof[6] + 1;//LI Y-displacement of BR corner in U.
      */
      edof[0] = 2*n1;//LI of X-displacement of TL corner in U.
      edof[1] = 2*n1 + 1;//LI of Y-displacement of TL corner in U.
      edof[2] = 2*n2;//LI X-displacement of TR corner in U.
      edof[3] = 2*n2 + 1;//LI Y-displacement of TR corner in U.
      edof[4] = 2*n2 + 2;//LI X-displacement of BL corner in U.
      edof[5] = 2*n2 + 3;//LI Y-displacement of BL corner in U.
      edof[6] = 2*n1 + 2;//LI X-displacement of BR corner in U.
      edof[7] = 2*n1 + 3;//LI Y-displacement of BR corner in U.

      //This is the subset step that Sigmund uses in MATLAB.
      for(w=0;w<8;w++){//Scan rows.
        for(l=0;l<8;l++){//Scan columns.
          ierr = MatSetValues(K,1,&edof[w],1,&edof[l],&KE[w][l],ADD_VALUES); CHKERRQ(ierr);
        }//end "l" loop
      }//end "w" loop
    }//end "j" loop
  }//end "i" loop
  ierr = MatAssemblyBegin(K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
  PetscReal zero = 0;
  PetscReal one = 1;
  PetscInt last = size-1;
  //Cleanup steps - X-displacements clamped on left edge.
  for(w=0;w<ny;w++){//Scan along the y-nodes
    n1 = w*2;//"n1" now stores the row index of the clamped X-nodes.
    for(l=0;l<size;l++){//Scan along the K matrix
      ierr = MatSetValues(K,1,&n1,1,&l,&zero,INSERT_VALUES); CHKERRQ(ierr);//Clears redundant coefficients.
    }
    ierr = MatSetValues(K,1,&n1,1,&n1,&one,INSERT_VALUES); CHKERRQ(ierr);//Effectively imposes the roller boundary conditions.
  }//end "w" loop
  //Cleanup step: Y-displacement clamped at BR corner
  for(w=0;w<size;w++){
    ierr = MatSetValues(K,1,&last,1,&w,&zero,INSERT_VALUES); CHKERRQ(ierr);
  }//end "w" loop
  ierr = MatSetValues(K,1,&last,1,&last,&one,INSERT_VALUES); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);


  //Linear solver setup.
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,K,K); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
  ierr = VecDuplicate(F,&U); CHKERRQ(ierr);
  ierr = VecSet(U,0.0); CHKERRQ(ierr);
  ierr = KSPSolve(ksp,F,U); CHKERRQ(ierr);

  //Cleanup (free-up memory)
  KSPDestroy(&ksp); MatDestroy(&K);
  VecDestroy(&F); VecDestroy(&U);

  return PetscFinalize();
}//end main function


//
/*
//Lame's constants for the elements.
mult = E/(1-nu*nu);
k[0] = (3-nu)*mult/6;
k[1] = (1+nu)*mult/8;
k[2] = -(3+nu)*mult/12;
k[3] = (3*nu-1)*mult/8;
k[4] = (nu-3)*mult/12;
k[5] = -(nu+1)*mult/8;
k[6] = nu*mult/6;
k[7] = (1-3*nu)*mult/8;
*/
