package test.sparse.eigensolvers;

// Import MTJ classes
// see http://code.google.com/p/matrix-toolkits-java/
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;

// Import sparse eigensolver package
//see http://code.google.com/p/sparse-eigensolvers-java/
import sparse.eigenvolvers.java.*;

// This example illustrates how the user can define their own operator in a totally flexible way
public class ExampleExtendOperator {
	public static void main(String[] args) {
		// Problem setup
		int blockSize=10;
		int maxIterations=200;
		int verbosityLevel=1;
		int n=500;
		double tolerance=1e-6;
		Lobpcg lobpcg=new Lobpcg();
//		MyOperator operA=new MyOperator(n);
		MyOperator2 operA=new MyOperator2(n);
		DenseMatrix blockVectorX=(DenseMatrix) Matrices.random(n,blockSize);  //Random matrix of vectors to define starting subspace
		
		// Run lobpcg eigensolver to obtain the smallest blocksize eigenvalues
		// Requires more than 150 iterations to converge
		lobpcg.runLobpcg(blockVectorX,operA,tolerance,maxIterations,verbosityLevel);
		// Now eigenvectors are stored in blockVectorX
		// Can obtain eigenvalues for further processing in dense matrix
		Utilities.print(lobpcg.getEigenvaluesMatrix());
		// Or get eigenvalues in double array
		Utilities.print(lobpcg.getEigenvalues());
		
		// Now run problem with preconditioner, which converges much faster
		// Requires less than 30 iterations to converge
		MyPreconditioner T=new MyPreconditioner(n);
		blockVectorX=(DenseMatrix) Matrices.random(n,blockSize); 
		maxIterations=50;
		lobpcg.runLobpcg(blockVectorX,operA,null,T,tolerance,maxIterations,verbosityLevel);
	}
}

// User defined matrix vector multiply
// User defined operatorAction returns same size dense matrix as input matrix. Otherwise user has total flexibility.
// Operator may exist as a dense or sparse matrix or as some user defined algorithm.
// Size of problem must be defined when constructor executed.
class MyOperator extends Operator {
	DenseMatrix myOperatorMatrix;
	
	// Constructor creates diagonal matrix with 1:n on diagonal in this example
	// There must be some kind of constructor at least to define the size of the problem
	public MyOperator(int n){
		myOperatorMatrix=(DenseMatrix) Matrices.identity(n);
		for (int i=0;i<n;i++) myOperatorMatrix.set(i, i, (double) (i+1));
		setOperatorSize(n); // Must specify operator size ( n x n)
	}
	
	// Override operator action
	// Multiply operator matrix times block matrix X in this case
	public DenseMatrix operatorAction(DenseMatrix X){
		DenseMatrix Temp=new DenseMatrix(X.numRows(),X.numColumns()); 
		myOperatorMatrix.mult(X, Temp);
		return Temp;
	}
}


// User defined matrix vector multiply
// In this case the user defined function for a matrix vector multiply does not store matrix
class MyOperator2 extends Operator {
	double[] Diag;
	
	// Constructor saves diagonal values
	// There must be some kind of constructor at least to define the size of the problem
	public MyOperator2(int n){
		Diag=new double[n];
		for (int i=0;i<n;i++) Diag[i]=(double)(i+1);
		setOperatorSize(n); // Must specify operator size ( n x n)
	}
	
	// Override operator action
	// Use 1D double array to implement equivalent of matrix vector multiply
	public DenseMatrix operatorAction(DenseMatrix X){
		DenseMatrix Temp=new DenseMatrix(X.numRows(),X.numColumns()); 
		for (int i=0;i<X.numRows();i++){
			for (int j=0;j<X.numColumns();j++) Temp.set(i,j,X.get(i,j)*Diag[i]);
		}
		return Temp;
	}
}

//User defined preconditioner which is inverse of operator action
class MyPreconditioner extends Operator {
	double[] Diag;
	
	// Constructor saves diagonal values
	// There must be some kind of constructor at least to define the size of the problem
	public MyPreconditioner(int n){
		Diag=new double[n];
		for (int i=0;i<n;i++) Diag[i]=(double)(i+1);
		setOperatorSize(n); // Must specify operator size ( n x n)
	}
	
	// Override operator action
	// Use 1D double array to implement equivalent of inverse of operator
	public DenseMatrix operatorAction(DenseMatrix X){
		DenseMatrix Temp=new DenseMatrix(X.numRows(),X.numColumns()); 
		for (int i=0;i<X.numRows();i++){
			for (int j=0;j<X.numColumns();j++) Temp.set(i,j,X.get(i,j)/Diag[i]);
		}
		return Temp;
	}
}

