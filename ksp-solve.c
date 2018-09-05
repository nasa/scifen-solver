#include "ksp-solve.h"
#include <stdio.h>
#include <ctype.h> 
#include <string.h>
int main(int argc, char **args){
  PetscInt        verbosity=1,maxits=1000,ndof,iters;
  KSP             ksp;// Krylov subspace solver
  PC              pc;// preconditioner
  PetscBool       flg,has_refu,has_u0,save_sol,save_csv;
  PetscErrorCode  ierr;
  char            kfile[PETSC_MAX_PATH_LEN];// input file names
  char            pfile[PETSC_MAX_PATH_LEN];
  char            ufile[PETSC_MAX_PATH_LEN];
  char            gfile[PETSC_MAX_PATH_LEN];
  char            sfile[PETSC_MAX_PATH_LEN];// output file name
  char            vfile[PETSC_MAX_PATH_LEN];// csv file name
  // -------------------- Get Started ---------------------
  ierr = PetscInitialize(&argc, &args, NULL,help); CHKERRQ(ierr);
  //
  PetscOptionsGetInt(PETSC_HACK_POGX,"-v",&verbosity,&flg);
  if(verbosity>2){PetscPrintf(PETSC_COMM_WORLD,"Parsing command line...\n");};
  //
  // -rmat and -rrhs are required
  PetscOptionsGetString(PETSC_HACK_POGX,"-rmat",kfile,PETSC_MAX_PATH_LEN,&flg);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Use -rmat <stiffness-filename>");
  //
  PetscOptionsGetString(PETSC_HACK_POGX,"-rrhs",pfile,PETSC_MAX_PATH_LEN,&flg);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Use -rrhs <loading-filename>");
  //
  // These are optional.
  PetscOptionsGetString(PETSC_HACK_POGX,"-rref",ufile,PETSC_MAX_PATH_LEN,&has_refu);
  PetscOptionsGetString(PETSC_HACK_POGX,"-ruin",gfile,PETSC_MAX_PATH_LEN,&has_u0  );
  PetscOptionsGetString(PETSC_HACK_POGX,"-wsol",sfile,PETSC_MAX_PATH_LEN,&save_sol);
  PetscOptionsGetString(PETSC_HACK_POGX,"-wcsv",vfile,PETSC_MAX_PATH_LEN,&save_csv);
  //
  // - Load the matrix and rhs, then destroy the viewers. -
  PetscViewer matsave, vecsave;
  if(verbosity>1){PetscPrintf(PETSC_COMM_WORLD,"Loading matrix: %s...\n",kfile);}
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,
    kfile, FILE_MODE_READ, &matsave); CHKERRQ(ierr);
  Mat k; Vec p,u;
  MatCreate(PETSC_COMM_WORLD,&k);
  MatSetType(k,MATSEQAIJ);
  ierr = MatLoad(k,matsave); CHKERRQ(ierr);
  PetscViewerDestroy(&matsave);
  //
  if(verbosity>1){PetscPrintf(PETSC_COMM_WORLD,"Loading RHS: %s...",pfile);}
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,
    pfile, FILE_MODE_READ, &vecsave); CHKERRQ(ierr);
  VecCreate(PETSC_COMM_WORLD,&p);
  ierr = VecLoad(p,vecsave); CHKERRQ(ierr);
  PetscViewerDestroy(&vecsave);
  ierr = VecGetSize(p,&ndof);CHKERRQ(ierr);
  maxits=ndof;
  if(verbosity>1){PetscPrintf(PETSC_COMM_WORLD,"%i DOF\n",ndof);}
  //
  // ----- Set Solver, Preconditioner, and Tolerances -----
  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
  ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
  //ierr = KSPSetTolerances(ksp, PETSC_DEFAULT,1e-16, PETSC_DEFAULT, 20000); 
  ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,maxits);
  //
  // Command-line options override above defaults
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
  ierr = PCSetFromOptions(pc); CHKERRQ(ierr);
  KSPType kt; PCType pt; PetscReal rtol,abstol,dtol;
  // Keep solver and preconditioner type strings after PETSc destruction
  char cpt[8],ckt[8];
  memset(cpt, '\0', sizeof(cpt));//FIXME There may be an off-by-one error here.
  memset(ckt, '\0', sizeof(ckt));
  if(verbosity>0){
    ierr = KSPGetType(ksp,&kt);
    ierr = PCGetType(pc,&pt);
    ierr = KSPGetTolerances(ksp,&rtol,&abstol,&dtol,&maxits);
    strncpy(ckt, kt, sizeof(ckt));
    strncpy(cpt, pt, sizeof(cpt));
  }
  if(verbosity>2){// Check settings
    PetscPrintf(PETSC_COMM_WORLD,"Solver: KSP-%s, ",kt);
    PetscPrintf(PETSC_COMM_WORLD,"Preconditioner: PC-%s\n",pt);
    PetscPrintf(PETSC_COMM_WORLD,"Using Preconditioned Norm Tolerances\n\
      Relative: %e, Absolute: %e,\n\
      Diverges: %e, Max Iter: %i\n" ,
      rtol,abstol,dtol,maxits);
  };
  // -------------------- Set up Solver -------------------
  ierr = KSPSetOperators(ksp,k,k); CHKERRQ(ierr);
  // Set these up now so solve timing is more accurate
  ierr = KSPSetUp(ksp); CHKERRQ(ierr);
  ierr = KSPSetUpOnBlocks(ksp); CHKERRQ(ierr);
  ierr = MatCreateVecs(k,&u,NULL); CHKERRQ(ierr);
  ierr = PCSetUp(pc); CHKERRQ(ierr);
  //
  // ------------------- Load Solutions -------------------
  if(has_u0){
    PetscViewer u0save;
    if(verbosity>1){PetscPrintf(PETSC_COMM_WORLD,"Loading initial guess: %s...\n",gfile);}
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,
      gfile, FILE_MODE_READ, &u0save); CHKERRQ(ierr);
    ierr = VecLoad(u,u0save); CHKERRQ(ierr);;
    PetscViewerDestroy(&u0save);
    KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
  }
  else{ VecSet(u,0.0);}
  Vec uref;
  if(has_refu){
    if(verbosity>1){PetscPrintf(PETSC_COMM_WORLD,"Loading reference: %s...\n",ufile);}
    PetscViewer refsave;
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,
      ufile, FILE_MODE_READ, &refsave); CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD,&uref); CHKERRQ(ierr);
    ierr = VecLoad(uref,refsave); CHKERRQ(ierr);
    PetscViewerDestroy(&refsave);
  }
  // ----------------------- Solve ------------------------
  if(verbosity>1){PetscPrintf(PETSC_COMM_WORLD,"Solving...\n");}
  ierr = KSPSolve(ksp,p,u); CHKERRQ(ierr);
  ierr=KSPGetIterationNumber(ksp, &iters);
  KSPConvergedReason reason;
  ierr=KSPGetConvergedReason(ksp,&reason);
  //
  // -------------------- Print Norm ----------------------
  PetscReal norm;
  ierr = KSPGetResidualNorm(ksp, &norm); CHKERRQ(ierr);
  if(verbosity>0){PetscPrintf(PETSC_COMM_WORLD,"Preconditioned Residual Norm: %e\n",norm);}
  if(has_refu){
    PetscScalar a=-1.0; PetscReal norm_ref;
    ierr = VecAXPY(uref,a,u);
    VecNorm(uref,NORM_2,&norm_ref);
    if(verbosity>0){PetscPrintf(PETSC_COMM_WORLD,"Reference Solution True Error: %e\n",norm_ref);}
  }
  // ------------------- Save Solution --------------------
  if(save_sol){
    PetscViewer solsave;
    if(verbosity>1){PetscPrintf(PETSC_COMM_WORLD,"Saving solution: %s...\n",sfile);}
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,
      sfile, FILE_MODE_WRITE, &solsave); CHKERRQ(ierr);
    VecView(u,solsave); CHKERRQ(ierr);
    PetscViewerDestroy(&solsave);
  }
  // --------------------- Finish Up ----------------------
  if(verbosity>2){PetscPrintf(PETSC_COMM_WORLD,"Destroying PETSc objects...\n");}
  MatDestroy(&k);
  VecDestroy(&p);
  VecDestroy(&u);
  if(has_refu){ VecDestroy(&uref);};
  KSPDestroy(&ksp);
  //
  PetscFinalize();
  //
  if(save_csv){
    FILE *fp;
    fp = fopen(vfile,"a");
    fprintf(fp,"%8s,%8s,%e,%e,%e,%e,%i,%i,%i\n",
	   ckt,cpt,rtol,abstol,dtol,norm,maxits,iters,reason);
    fclose(fp);
  };
  //
  if(verbosity>2){printf("Done.\n");}
  if(verbosity>2){//FIXME This just dumps a csv line into stdout
    printf("%8s,%8s,%e,%e,%e,%e,%i,%i,%i\n",
	   ckt,cpt,rtol,abstol,dtol,norm,maxits,iters,reason);
  }
  return ierr;
}
