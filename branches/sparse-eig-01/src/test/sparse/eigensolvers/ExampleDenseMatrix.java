package test.sparse.eigensolvers;

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