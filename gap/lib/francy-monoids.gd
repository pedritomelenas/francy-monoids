#! @Arguments f
#! @Returns a drawing
#! @Description
#! f is a set of factorizations 
#! Draws the graph of factorizations associated to f: a complete graph 
#! whose vertices are the elements of f. Edges are labelled with
#! distances between nodes they join. Kruskal algorithm is used to 
#! draw in red a spannin tree with minimal distances. Thus the catenary
#! degree is reached in the edges of the tree.
DeclareOperation("DrawFactorizationGraph",[IsRectangularTable]);

#! @Arguments f
#! @Returns a drawing
#! @Description
#! f is a set of factorizations 
#! Draws the Eliahou graph of factorizations associated to f: a graph 
#! whose vertices are the elements of f, and there is an edge between
#! two vertices if they have common support. Edges are labelled with
#! distances between nodes they join.
DeclareOperation("DrawEliahouGraph",[IsRectangularTable]);

#! @Arguments n,s 
#! @Returns a drawing
#! @Description
#! s is either a numerical semigroup or an affine semigroup, and n is an
#! element of s
#! Draws the graph associated to n in s.
DeclareOperation("DrawRosalesGraph",[IsHomogeneousList,IsAffineSemigroup]);
DeclareOperation("DrawRosalesGraph",[IsInt,IsNumericalSemigroup]);


#! @Arguments s 
#! @Returns a drawing
#! @Description
#!  Draws the Hasse diagram of 
#!  oversemigroupstree of the numerical semigroup s.
#!  Irreducible numerical semigroups are highlighted.
DeclareGlobalFunction("DrawOverSemigroupsNumericalSemigroup");


#! @Arguments s,l,generators 
#! @Returns a drawing
#! @Description
#!  Draws the tree of sons of numerical semigroups up to level l with 
#!  respect to the minimal system of generators given by the function
#!  generators. 
DeclareGlobalFunction("DrawTreeOfSonsOfNumericalSemigroup");


#! @Arguments s,A 
#! @Returns a drawing
#! @Description
#! plots a graph whose set of vertices is A, which is a finite set of 
#! integers, and whose edges are provided by the order of the numerical 
#! semigroup s.
DeclareGlobalFunction("DrawHasseDiagramOfNumericalSemigroup");



#! @Arguments s,[expand] 
#! @Returns a drawing
#! @Description
#! Returns a Francy canvas with the tree of gluings of the numerical 
#! semigroup s. If the optional argument expand is provided, then
#! the tree is drawn fully expanded.
DeclareGlobalFunction("DrawTreeOfGluingsOfNumericalSemigroup"); 
