#############################################################################
##
#A  testall.tst           FrancyMonoids package       Pedro A. Garcia-Sanchez
#W                                                    Andres Herrera-Poyatos
##                                                    Manuel Martins
##
##  (based on the cooresponding file of the 'example' package,
##   by Alexander Konovalov)
##
##  To create a test file, place GAP prompts, input and output exactly as
##  they must appear in the GAP session. Do not remove lines containing 
##  START_TEST and STOP_TEST statements.
##
##  The first line starts the test. START_TEST reinitializes the caches and 
##  the global random number generator, in order to be independent of the 
##  reading order of several test files. Furthermore, the assertion level 
##  is set to 2 by START_TEST and set back to the previous value in the 
##  subsequent STOP_TEST call.
##
##  The argument of STOP_TEST may be an arbitrary identifier string.
## 
gap> START_TEST("FrancyMonoids package: testall.tst");

# Note that you may use comments in the test file
# and also separate parts of the test by empty lines

# First load the package without banner (the banner must be suppressed to 
# avoid reporting discrepancies in the case when the package is already 
# loaded)
gap> LoadPackage("FrancyMonoids",false);
true

# Check that the data are consistent  

gap> s:=NumericalSemigroup(3,5,7);
<Numerical semigroup with 3 generators>
gap> DrawTreeOfSonsOfNumericalSemigroup(s,3,MinimalGenerators);
<jupyter renderable>
gap> DrawTreeOfGluingsOfNumericalSemigroup(s);
<jupyter renderable>
gap> DrawRosalesGraph(10,s);
<jupyter renderable>
gap> DrawEliahouGraph(FactorizationsIntegerWRTList(10,[3,5,7]));
<jupyter renderable>
gap> DrawFactorizationGraph(FactorizationsIntegerWRTList(10,[3,5,7]));
<jupyter renderable>
gap> DrawOverSemigroupsNumericalSemigroup(s);
<jupyter renderable>
gap> DrawHasseDiagramOfNumericalSemigroup(s,[1,2,3]);
<jupyter renderable>



## Each test file should finish with the call of STOP_TEST.
## The first argument of STOP_TEST should be the name of the test file.
## The second argument is redundant and is used for backwards compatibility.
gap> STOP_TEST( "testall.tst", 10000 );

#############################################################################
##
#E
