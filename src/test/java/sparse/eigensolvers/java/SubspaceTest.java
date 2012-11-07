/*
 * Copyright (C) 2012 Rico Argentati
 * 
 * This file is part of SEJ (Sparse Eigensolvers for Java).
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 * 
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

package sparse.eigensolvers.java;

import static org.junit.Assert.*;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;

import org.junit.Test;

import sparse.eigenvolvers.java.Subspace;
import sparse.eigenvolvers.java.Utilities;

public class SubspaceTest  extends TestCase {

	@Test
	public void testAccuarcy() {
		int n=100;
		int p=10;
		int q=10;
		
//		DenseMatrix X=(DenseMatrix) Matrices.random(n, p);
//		DenseMatrix Y=(DenseMatrix) Matrices.random(n, q);
		
		DenseMatrix X=new DenseMatrix(2,1);
		DenseMatrix Y=new DenseMatrix(2,1);
		X.set(0, 0, 1);
		Y.set(0, 0, 1);
//		Y.set(1, 0, 1);
//		Y.set(1, 0, 1e-4);
//		Y.set(1, 0, 1e-6);
//		Y.set(1, 0, 1e-8);
//		Y.set(1, 0, 1e-10);
//		Y.set(1, 0, 1e-16);
		Y.set(1, 0, 1e-20);
//		Y.set(1, 0, 1e-30);
		
		double[] s=Subspace.getPrincipleAngles(X, Y);
		System.out.println("Check accuracy");
		Utilities.print(s);
		fail("Not yet implemented");
	}
	
	public void testAllPrincipalAngles() {
		int n=100;
		int p=10;
		int q=10;
		
		DenseMatrix X=(DenseMatrix) Matrices.random(n, p);
		DenseMatrix Y=(DenseMatrix) Matrices.random(n, q);
		
		double[] s=Subspace.getPrincipleAngles(X, Y);
		Utilities.print(s);
		System.out.println(Subspace.getPrincipleAnglesLargest(X, Y));
		s=Subspace.getPrincipleAnglesSines(X, Y);
		Utilities.print(s);
	}
	
//	public void testMisc() {
//		System.out.println("running testMisc");
//		int n=5;
//		int m=2;
//		DenseMatrix X=(DenseMatrix) Matrices.random(n, m);
//		DenseMatrix Y=Utilities.blockMatrix(1,2,X,X);
//		DenseMatrix Z=new DenseMatrix(0,0);
////		Y=Utilities.orth(X);
//		Utilities.print(Utilities.orth(X));
//		Utilities.print(Utilities.orth(Y));
////		Utilities.print(Y);
//		
//	}

}
