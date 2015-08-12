# Lobpcg Java Examples #

**These java source files are also provided in the download file.**


`ExampleDenseMatrix.java`
```
//Import MTJ classes
//see http://code.google.com/p/matrix-toolkits-java/
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;

//Import sparse eigensolver package
//see http://code.google.com/p/sparse-eigensolvers-java/
import sparse.eigenvolvers.java.*;

// This example illustrates how to find eigenvalues for a symmetric dense matrix
// Loppcg is really designed to handle sparse matrices, but this example is included as an illustration.
// There are many java programs in existence to compute eigenvalues for dense matrices.
public class ExampleDenseMatrix {
    public static void main(String[] args) {
        // Setup problem
	int n=500;
	int maxIterations=200;
	DenseMatrix A=(DenseMatrix) Matrices.random(n,n);
	//Make it symmetric
	DenseMatrix Temp=A.copy();
	A.add(Temp.transpose());
	Lobpcg lobpcg=new Lobpcg();
	
	// Run lobpcg with mostly default parameters
	lobpcg.runLobpcg(A,maxIterations);
	
	// Now eigenvectors are transferred to blockVectorX
	DenseMatrix blockVectorX=lobpcg.getEigenvectors();
	// Can obtain eigenvalues for further processing in dense matrix
	Utilities.print(lobpcg.getEigenvaluesMatrix());
	// Or get eigenvalues in double array
	Utilities.print(lobpcg.getEigenvalues());
	
	// Now run lobpcg again to compute largest eigenvalues
	int verbosityLevel=0;
	A.scale(-1);
	lobpcg.runLobpcg(A,maxIterations,verbosityLevel); // no output
	// Print negative of largest eigenvalues
	Utilities.print(lobpcg.getEigenvaluesMatrix());
    }
}
```


`ExampleLaplacianOperator.java`
```
//Import MTJ classes
//see http://code.google.com/p/matrix-toolkits-java/
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix.Norm;
import no.uib.cipr.matrix.sparse.CompRowMatrix;

//Import sparse eigensolver package
//see http://code.google.com/p/sparse-eigensolvers-java/
import sparse.eigenvolvers.java.*;

public class ExampleLaplacianOperator {
	
// This example illustrates how to find eigenvalues for a sparse operator stored
// in compressed row format.
// We generate a 3D Laplacian for testing purposes and check the computed eigenvalues
// against those computed theoretically.
public static void main(String[] args) {
	// Setup problem
	int blockSize=20;
	int maxIterations=100;
	double tolerance=1e-8;
	int verbosityLevel=1;
	int innerIterations=5;
	int nx=20,ny=20,nz=20;
	double error;
	Lobpcg lobpcg=new Lobpcg();
	CompRowMatrix ACompRow=Utilities.getLaplacian(nx,ny,nz); //Laplacian 20 x 20 x 20
	int n=ACompRow.numColumns();
	DenseMatrix blockVectorX=(DenseMatrix) Matrices.random(n,blockSize);
				
	// Run without preconditioner which takes many more iterations
	maxIterations=300;
	lobpcg.runLobpcg(blockVectorX,ACompRow,tolerance,maxIterations,verbosityLevel);
		
	// Run lobpcg with preconditioner which takes many fewer iterations
	// We need place holder "null" since we are not solving
	// a generalized eigenvalue problem and there is no operB
	// Use conjugate gradient preconditioner (operT approximates inverse of ACompRow)
	OperatorPrecCG operT= new OperatorPrecCG(ACompRow); 
	operT.setCGNumberIterations(innerIterations); // Need to try different values to see best convergence
	maxIterations=80;
	blockVectorX=(DenseMatrix) Matrices.random(n,blockSize);
	lobpcg.runLobpcg(blockVectorX,ACompRow,null,operT,tolerance,maxIterations,verbosityLevel);
	error=Utilities.sub(Utilities.getLaplacianEigenvalues(blockSize,nx,ny,nz),
	lobpcg.getEigenvaluesMatrix()).norm(Norm.Frobenius);
		
	System.out.println("Test results against theoretically computed eigenvalues:");
	System.out.printf("Error= %e (Laplacian 20 x 20 x 20 with preconditioner)\n",error);
	}
}
```

`ExampleMatrixMarket.java`
```
//Import MTJ classes
//see http://code.google.com/p/matrix-toolkits-java/
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.sparse.CompRowMatrix;

//Import sparse eigensolver package
//see http://code.google.com/p/sparse-eigensolvers-java/
import sparse.eigenvolvers.java.*;

// This example illustrates how user can find several of the largest eigenvalues 
// after loading a matrix market file
// Note that string file needs to be defined based on where user stores matrix market files
// See http://math.nist.gov/MatrixMarket/
public class ExampleMatrixMarket {
	public static void main(String[] args) {
		// Problem setup
		int blockSize=10;
		int maxIterations=100;
		int verbosityLevel=1;
		double tolerance=1e-6;
		String file="C:\\MatrixMarket\\plat1919.mtx";
		Lobpcg lobpcg=new Lobpcg();
		CompRowMatrix A=Utilities.readMatrixMarketFile(file);
		A.scale(-1); // Multiply operator times -1 to compute the largest eigenvalues
		
		// Run lobpcg
		// Note that if the first parameter is an integer, a random set of starting vectors
		// defining the initial subspace is generated inside lobpcg
		lobpcg.runLobpcg(blockSize,A,tolerance,maxIterations,verbosityLevel);
		
		// Now eigenvectors are transferred to blockVectorX
		DenseMatrix blockVectorX=lobpcg.getEigenvectors();
		// Can obtain eigenvalues for further processing in dense matrix
		Utilities.print(lobpcg.getEigenvaluesMatrix());
		// Or get eigenvalues in double array
		Utilities.print(lobpcg.getEigenvalues());
	}
}
```

`ExampleExtendOperator.java`
```
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
```