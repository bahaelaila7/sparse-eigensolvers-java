rem 10/11/2012

set CLASSPATH=.;D:\MyMainFiles\MyDocuments\Java\Lobpcg Project\ReleasesLobpcg\01\sparse-eig-01\lib\arpack_combined_all.jar
set CLASSPATH=%CLASSPATH%;D:\MyMainFiles\MyDocuments\Java\Lobpcg Project\ReleasesLobpcg\01\sparse-eig-01\lib\netlib-java-0.9.3.jar
set CLASSPATH=%CLASSPATH%;D:\MyMainFiles\MyDocuments\Java\Lobpcg Project\ReleasesLobpcg\01\sparse-eig-01\lib\mtj-0.9.14.jar
set CLASSPATH=%CLASSPATH%;D:\MyMainFiles\MyDocuments\Java\Lobpcg Project\ReleasesLobpcg\01\sparse-eig-01\sparse-eig-01.jar

java test/sparse/eigensolvers/ExampleDenseMatrix
java test/sparse/eigensolvers/ExampleExtendOperator 
java test/sparse/eigensolvers/ExampleLaplacianOperator
java test/sparse/eigensolvers/ExampleMatrixMarket
java test/sparse/eigensolvers/TestLobpcg

pause

