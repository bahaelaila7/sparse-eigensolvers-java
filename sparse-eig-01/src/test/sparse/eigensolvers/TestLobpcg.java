package test.sparse.eigensolvers;

//Import MTJ classes
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.QR;
import no.uib.cipr.matrix.SVD;
import no.uib.cipr.matrix.Matrix.Norm;
import no.uib.cipr.matrix.sparse.CompRowMatrix;

//Import sparse eigensolver package
//see http://code.google.com/p/sparse-eigensolvers-java/
import sparse.eigenvolvers.java.*;

public class TestLobpcg {

	public static void main(String[] args) {
		int[] tests;
		// Enter test numbers as args
		if (args.length>0){
			tests= new int [args.length];
			int i=0;
			for (String s: args) tests[i++]=Integer.parseInt(s);
	    }
		// Select tests to run as default
		else {
			tests=new int[]{1,2,3,4};
//			tests=new int[]{4};
		}
		
		for(int i=0;i<tests.length;i++){
			switch (tests[i]) {
	        case 1:	testInputParameters();
	                break;
	        case 2:	testLaplacian();
            		break;
	        case 3:	testGeneralizeEigenvalues();
	        		break;
	        case 4:	testConstraints();
    				break;		
	        default: break;
			}
		}
	}
	
	
	// Test various combinations of input parameters
	static void testInputParameters(){
		int blockSize=20;
		int maxIterations=200;
		double tolerance=1e-8;
		int verbosityLevel=0;
		int nx=20,ny=20;
		DenseMatrix blockVectorX;
		Lobpcg lobpcg;
		Operator operA;
		
		// Problem setup
		operA=new Operator(Utilities.getLaplacian(nx,ny));
		DenseMatrix A=operA.getDenseOperatorMatrix();
		int n=A.numColumns();
		DenseMatrix B=Utilities.mult(A,A);  // Make B=A^2
		lobpcg=new Lobpcg();
		lobpcg.setInformationString("Running test: testInputParameters()");
		operA=new Operator(A);
		blockVectorX=(DenseMatrix) Matrices.random(n, blockSize);
		
		lobpcg.setVerbosityLevel(0);
		lobpcg.setMaxIterations(200);
		if (lobpcg.runLobpcg(A)>0)
			System.out.println("1 lobpcg.runLobpcg(A): converged");
		if (lobpcg.runLobpcg(5,A)>0)
			System.out.println("2 lobpcg.runLobpcg(5,A): converged");
		if (lobpcg.runLobpcg(25,A,B)>0)
			System.out.println("3 lobpcg.runLobpcg(25,A,B): converged");
		if (lobpcg.runLobpcg(A,B)>0)
			System.out.println("4 lobpcg.runLobpcg(A,B): converged");
		
		if (lobpcg.runLobpcg(A,B,tolerance,maxIterations,verbosityLevel)>0)
			System.out.println("5 lobpcg.runLobpcg(A,B,tolerance,maxIterations,verbosityLevel): converged");
				
		if (lobpcg.runLobpcg(A,B,tolerance,maxIterations,verbosityLevel)>0)
			System.out.println("6 lobpcg.runLobpcg(A,B,tolerance,maxIterations,verbosityLevel): converged");
		
		if (lobpcg.runLobpcg(blockVectorX,A,B,tolerance,maxIterations,verbosityLevel)>0)
			System.out.println("7 lobpcg.runLobpcg(blockVectorX,A,B,tolerance,maxIterations,verbosityLevel): converged");
		
		if (lobpcg.runLobpcg(blockVectorX,operA,B,tolerance,maxIterations,verbosityLevel)>0)
			System.out.println("8 lobpcg.runLobpcg(blockVectorX,operA,B,tolerance,maxIterations,verbosityLevel): converged");
		
		if (lobpcg.runLobpcg(blockVectorX,tolerance,maxIterations,A,B,verbosityLevel)>0)
			System.out.println("9 lobpcg.runLobpcg(blockVectorX,tolerance,maxIterations,A,B,verbosityLevel): converged");
		
		// Run lobpcg() with no in line parameters
		lobpcg=new Lobpcg();
		lobpcg.setInformationString("Running test: testInputParameters() (no inline)");
		blockVectorX=(DenseMatrix) Matrices.random(n, blockSize);
		lobpcg.setBlockVectorX(blockVectorX);
		lobpcg.setA(A);
		OperatorPrecDenseCholesky operT = new OperatorPrecDenseCholesky(A);
		lobpcg.setT(operT);
		lobpcg.setMaxIterations(50);
		lobpcg.setResidualTolerance(1e-8);
		lobpcg.setVerbosityLevel(0);
		if (lobpcg.runLobpcg()>0) 
			System.out.println("10 Run lobpcg() with no in line parameters has converged");
		
		// run with simplist parameters
		lobpcg=new Lobpcg();
		lobpcg.setVerbosityLevel(0);
		lobpcg.setA(A);
		lobpcg.setMaxIterations(100);
		if (lobpcg.runLobpcg()>0) 
			System.out.println("11 Run lobpcg with simplist parameters and default block size has converged");
		// run with simplist parameters and change block size
		lobpcg=new Lobpcg();
		lobpcg.setVerbosityLevel(0);
		lobpcg.setBlockSize(25);
		lobpcg.setMaxIterations(100);
		lobpcg.setA(A);
		if (lobpcg.runLobpcg()>0) 
			System.out.println("12 Run lobpcg with simplist parameters and block size 25 has converged");
		System.out.println("Testing input parameters complete!");
	}
	
	// Test Laplacian against theoretical values
	static void testLaplacian(){
		int verbosityLevel=1;
		int blockSize;
		int maxIterations;
		double tolerance;
		int innerIterations;
		int nx,ny,nz;
		double error1,error2,error3;
		DenseMatrix blockVectorX;
		Lobpcg lobpcg;
		
		// Problem setup no preconditioner, Laplacian 19 x 19
		lobpcg=new Lobpcg();
		blockSize=20; maxIterations=200; tolerance=1e-10;
		nx=19;ny=19;
		CompRowMatrix ACompRow=Utilities.getLaplacian(nx,ny);
		Operator operA=new Operator(ACompRow);
		int n=operA.getOperatorSize();
		blockVectorX=(DenseMatrix) Matrices.random(n,blockSize);
		// Run lobpcg
		lobpcg.runLobpcg(blockVectorX,operA,tolerance,maxIterations,verbosityLevel);
		error1=Utilities.sub(Utilities.getLaplacianEigenvalues(blockSize,nx,ny),
			lobpcg.getEigenvaluesMatrix()).norm(Norm.Frobenius);
		
		// Problem setup with preconditioner, Laplacian 19 x 19
		lobpcg=new Lobpcg();
		blockSize=20; maxIterations=75; tolerance=1e-10;
		innerIterations=4;
		nx=19;ny=19;
		ACompRow=Utilities.getLaplacian(nx,ny);
		operA=new Operator(ACompRow);
		OperatorPrecCG operT= new OperatorPrecCG(ACompRow);
		operT.setCGNumberIterations(innerIterations);
		n=operA.getOperatorSize();
		blockVectorX=(DenseMatrix) Matrices.random(n,blockSize);
		// Run lobpcg
		lobpcg.runLobpcg(blockVectorX,operA,null,operT,tolerance,maxIterations,verbosityLevel);
		error2=Utilities.sub(Utilities.getLaplacianEigenvalues(blockSize,nx,ny),
				lobpcg.getEigenvaluesMatrix()).norm(Norm.Frobenius);
		
		// Problem setup with CG preconditioner, Laplacian 20 x 20 x 20
		lobpcg=new Lobpcg();
		blockSize=20; maxIterations=100; tolerance=1e-10;
		innerIterations=4;
		nx=20;ny=20;nz=20;
		ACompRow=Utilities.getLaplacian(nx,ny,nz);
		operA=new Operator(ACompRow);
		operT= new OperatorPrecCG(ACompRow);
		operT.setCGNumberIterations(innerIterations);
		n=operA.getOperatorSize();
		blockVectorX=(DenseMatrix) Matrices.random(n,blockSize);
		// Run lobpcg
		lobpcg.runLobpcg(blockVectorX,operA,null,operT,tolerance,maxIterations,verbosityLevel);
		error3=Utilities.sub(Utilities.getLaplacianEigenvalues(blockSize,nx,ny,nz),
				lobpcg.getEigenvaluesMatrix()).norm(Norm.Frobenius);
		
		System.out.println("Test results against theoretically computed eigenvalues:");
		System.out.printf("Error= %e (Laplacian 19 x 19)\n",error1);
		System.out.printf("Error= %e (Laplacian 19 x 19 with preconditioner)\n",error2);
		System.out.printf("Error= %e (Laplacian 20 x 20 x 20 with preconditioner)\n",error3);
		
		/* Known theoretical results for 2D Laplacian - Utilities.getLaplacian(19,19)
		First 20 eigenvalues: 
		4.9246637619450060e+00
		1.2251028621941822e+01
		1.2251028621942096e+01
		1.9577393481938643e+01
		2.4261027043297744e+01
		2.4261027043299325e+01
		3.1587391903295610e+01
		3.1587391903295863e+01
		4.0658933005982500e+01
		4.0658933005982700e+01
		4.3597390324652710e+01
		4.7985297865978750e+01
		4.7985297865980440e+01
		5.9995296287336124e+01
		5.9995296287337425e+01
		6.1040975643662440e+01
		6.1040975643662820e+01
		6.8367340503659480e+01
		6.8367340503660340e+01
		7.6393202250021080e+01
		*/
	}

	//Test generalized eigenvalues using Laplacian 
	static void testGeneralizeEigenvalues(){
		// Problem setup
		int blockSize=30;
		int maxIterations=200;
		double tolerance=1e-10;
		int verbosityLevel=1;
		int nx=19,ny=19;
		Operator operA=new Operator(Utilities.getLaplacian(nx,ny));
		DenseMatrix A=operA.getDenseOperatorMatrix();
		DenseMatrix B=Utilities.mult(A,A);  // Make B=A^2 for testing purposes
		Lobpcg lobpcg=new Lobpcg();
		int n=A.numColumns();
		DenseMatrix blockVectorX=(DenseMatrix) Matrices.random(n, blockSize);
		
		// Run lobpcg
		if (lobpcg.runLobpcg(blockVectorX,operA,B,tolerance,maxIterations,verbosityLevel)>0)
			System.out.println("lobpcg.runLobpcg(blockVectorX,A,B,tolerance,maxIterations,verbosityLevel): converged");
		System.out.println("Testing generalized eigenvalues complete!");

		/* Matlab results, A=getLaplacian(19,19), B=A^2
		Final Eigenvalues lambda 1.2577424483213851e-003 
		Final Eigenvalues lambda 1.2694399311631438e-003 
		Final Eigenvalues lambda 1.2694399311631447e-003 
		Final Eigenvalues lambda 1.2813570386711460e-003 
		Final Eigenvalues lambda 1.2890934126830505e-003 
		Final Eigenvalues lambda 1.2890934126830525e-003 
		Final Eigenvalues lambda 1.3013841645270752e-003 
		Final Eigenvalues lambda 1.3013841645270758e-003 
		Final Eigenvalues lambda 1.3169312756370118e-003 
		Final Eigenvalues lambda 1.3169312756370125e-003
		
		Final Eigenvalues lambda 1.3220472632018111e-003 
		Final Eigenvalues lambda 1.3297612362660755e-003 
		Final Eigenvalues lambda 1.3297612362660764e-003 
		Final Eigenvalues lambda 1.3513427617188377e-003 
		Final Eigenvalues lambda 1.3513427617188381e-003 
		Final Eigenvalues lambda 1.3532550074356823e-003 
		Final Eigenvalues lambda 1.3532550074356825e-003 
		Final Eigenvalues lambda 1.3668061246587936e-003 
		Final Eigenvalues lambda 1.3668061246587945e-003 
		Final Eigenvalues lambda 1.3819660112501049e-003
		
		Final Eigenvalues lambda 1.3896171620084673e-003 
		Final Eigenvalues lambda 1.3896171620084695e-003 
		Final Eigenvalues lambda 1.3984161454712112e-003 
		Final Eigenvalues lambda 1.3984161454712121e-003 
		Final Eigenvalues lambda 1.4128916540949296e-003 
		Final Eigenvalues lambda 1.4128916540949311e-003 
		Final Eigenvalues lambda 1.4220204745782987e-003 
		Final Eigenvalues lambda 1.4220204745782998e-003 
		Final Eigenvalues lambda 1.4372806196115639e-003 
		Final Eigenvalues lambda 1.4372806196115652e-003 
		*/
	}
	
	// Test lobpcg constraints
	static void testConstraints(){
		int nn;
		DenseMatrix Y;
		DenseMatrix blockVectorX;
		Lobpcg lobpcg;
		DenseMatrix eigVal1,eigVal2,Temp;

		// Problem setup
		int blockSize=20;
		int maxIterations=100;
		double tolerance=1e-8;
		int verbosityLevel=0;
		int sizeLaplacian=20;
		Operator operA,operB;
		
		operA=new Operator(Utilities.getLaplacian(sizeLaplacian,sizeLaplacian));
		DenseMatrix A=operA.getDenseOperatorMatrix();
		int n=A.numRows();
		DenseMatrix B=Utilities.mult(A,A); // Make B=A^2
		operA=new Operator(A);
		operB=new Operator(B);
		OperatorPrecDenseCholesky operT = new OperatorPrecDenseCholesky(A);
		lobpcg=new Lobpcg();
	
		// Test eigenvalue problem with random constraints
		System.out.println("Test lobpcg eigenvalue problem with random constraints");
		int yColumns=100;
		Y=(DenseMatrix) Matrices.random(n,yColumns);
		// Orthonormalize Y
		QR qr=QR.factorize(Y);
		Y=qr.getQ();
		// Setup and execute SVD
		DenseMatrix U;
		try {
			SVD svd=SVD.factorize(Utilities.mult(Y,Utilities.trans(Y)));
			U=svd.getU();
			DenseMatrix Ycomp=Utilities.getMatrix(U, yColumns, n-1);
			// Get singular vectors for orthogonal complement of Y (null space of Y*Y^T)
			// Compute RITZ values for this subspace and store first blockSize of these
			nn=Ycomp.numColumns();
			double[] eigValues=new double[nn];
			Temp=Utilities.mult(Utilities.trans(Ycomp),Utilities.mult(A,Ycomp));
			Utilities.eig(Temp,eigValues);
			eigVal1=new DenseMatrix(blockSize,1);
			for (int i=0; i<blockSize; ++i) eigVal1.set(i,0,eigValues[i]);
			// Compute Ritz values using constraints
			// This operator has Ritz values plus yColumns zeros
			Temp=Utilities.mult(Ycomp,Utilities.trans(Ycomp));
			operA=new Operator(Utilities.mult(Temp,Utilities.mult(A,Temp)));
			blockVectorX=(DenseMatrix)Matrices.random(n,blockSize);
			// The constraints should remove the zeros and only keep the Ritz values
			lobpcg.runLobpcg(blockVectorX,operA,null,null,Y,tolerance,maxIterations,0);
			eigVal2=lobpcg.getEigenvaluesMatrix();
			System.out.printf("Check eigenvalues: %e\n\n",(Utilities.sub(eigVal2,eigVal1)).norm(Norm.Frobenius));
		} catch (NotConvergedException e){}
		
		// Test eigenvalue problem with a preconditioner and constraints
		operA=new Operator(A); // Restore operator A
		System.out.println("Test eigenvalue problem with a preconditioner and constraints");
		blockVectorX=(DenseMatrix) Matrices.random(n, 2*blockSize);
		lobpcg.runLobpcg(blockVectorX,operA,null,operT,tolerance,maxIterations,verbosityLevel);
		Y=Utilities.getMatrix(lobpcg.getEigenvectors(),0,blockSize-1);
		eigVal1=Utilities.getMatrix(lobpcg.getEigenvaluesMatrix(),blockSize,2*blockSize-1,0,0);
		nn=Y.numColumns();
		System.out.printf("Check normalization of eigenvectors: %e\n",
			(Utilities.sub(Matrices.identity(nn),
			Utilities.mult(Utilities.trans(Y),Y))).norm(Norm.Frobenius));
		// Now solve constrained problem
		blockVectorX=(DenseMatrix) Matrices.random(n, blockSize);
		lobpcg.runLobpcg(blockVectorX,operA,null,operT,Y,tolerance,maxIterations,verbosityLevel);
		eigVal2=lobpcg.getEigenvaluesMatrix();
		System.out.printf("Check eigenvalues: %e\n\n",(Utilities.sub(eigVal2,eigVal1)).norm(Norm.Frobenius));
		
		// Test generalized eigenvalue problem with constraints 
		System.out.println("Test generalized eigenvalue problem with constraints");
		blockVectorX=(DenseMatrix) Matrices.random(n, 2*blockSize);
		lobpcg.runLobpcg(blockVectorX,operA,operB,tolerance,maxIterations,verbosityLevel);
		Y=Utilities.getMatrix(lobpcg.getEigenvectors(),0,n-1,0,blockSize-1);
		eigVal1=Utilities.getMatrix(lobpcg.getEigenvaluesMatrix(),blockSize,2*blockSize-1,0,0);
		nn=Y.numColumns();
		Temp=Utilities.mult(Utilities.trans(Y),Utilities.mult(B,Y));
		System.out.printf("Check normalization of generalize eigenvectors: %e\n",
				(Utilities.sub(Matrices.identity(nn),Temp)).norm(Norm.Frobenius));
		// Now solve constrained problem
		blockVectorX=(DenseMatrix) Matrices.random(n, blockSize);
		lobpcg.runLobpcg(blockVectorX,operA,operB,null,Y,tolerance,maxIterations,verbosityLevel);
		eigVal2=lobpcg.getEigenvaluesMatrix();
		System.out.printf("Check generalized eigenvalues: %e\n\n",
			(Utilities.sub(eigVal2,eigVal1)).norm(Norm.Frobenius));

		// Test generalized eigenvalue problem with preconditioner and constraints 
		System.out.println("Test generalized eigenvalue problem with preconditioner and constraints");
		blockVectorX=(DenseMatrix) Matrices.random(n, 2*blockSize);
		lobpcg.runLobpcg(blockVectorX,operA,operB,operT,tolerance,maxIterations,verbosityLevel);
		Y=Utilities.getMatrix(lobpcg.getEigenvectors(),0,n-1,0,blockSize-1);
		eigVal1=Utilities.getMatrix(lobpcg.getEigenvaluesMatrix(),blockSize,2*blockSize-1,0,0);
		nn=Y.numColumns();
		Temp=Utilities.mult(Utilities.trans(Y),Utilities.mult(B,Y));
		System.out.printf("Check normalization of generalize eigenvectors: %e\n",
				(Utilities.sub(Matrices.identity(nn),Temp)).norm(Norm.Frobenius));
		// Now solve constrained problem
		blockVectorX=(DenseMatrix) Matrices.random(n, blockSize);
		lobpcg.runLobpcg(blockVectorX,operA,operB,operT,Y,tolerance,maxIterations,verbosityLevel);
		eigVal2=lobpcg.getEigenvaluesMatrix();
		System.out.printf("Check generalized eigenvalues: %e\n\n",
				(Utilities.sub(eigVal2,eigVal1)).norm(Norm.Frobenius));
		System.out.println("Test constraints complete!");
	}
}
