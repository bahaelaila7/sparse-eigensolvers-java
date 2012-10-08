package sparse.eigenvolvers.java;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix.Norm;
import no.uib.cipr.matrix.sparse.CompRowMatrix;

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
			tests=new int[]{1,2,3,4,5,6,7,8,9,10,11,12};
//			tests=new int[]{1};
			tests=new int[]{2};
		}
		
		for(int i=0;i<tests.length;i++){
			switch (tests[i]) {
	        case 1:	testLaplacian(0);
	                break;
	        case 2:	testLaplacian(1);
	        		break;
//	        case 3:	testLobpcgPerformance();
//	                break;
//	        case 4:	testLobpcgSparse();
//	        		break;
//	        case 5: testGeneralizeEigenvalues();
//					break;	
//	        case 6: testLaplacian();
//	        		break;
	        default: break;
			}
		}
	}
	
	// Test Laplacian 
	static void testLaplacian(int verbosityLevel){
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
		
		// Problem setup with preconditioner, Laplacian 20 x 20 x 20
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

}
