############################################################################
##
#O DrawFactorizationGraph(f)
##
## f is a set of factorizations 
## Draws the graph of factorizations associated to f: a complete graph 
## whose vertices are the elements of f. Edges are labelled with
## distances between nodes they join. Kruskal algorithm is used to 
## draw in red a spannin tree with minimal distances. Thus the catenary
## degree is reached in the edges of the tree.
##
#############################################################################
InstallMethod(DrawFactorizationGraph,"draws the factorization graph", [IsRectangularTable],
function(f)
    local graph, canvas, fs, c, nf, i, p, ln, distance, Kruskal, tv;

    Kruskal := function(V, E)
        local trees, needed, v, e, i,j, nv;

        trees := List(V, v-> [v]);
        needed := [];
        nv:=Length(V);
        for e in E do
          i:=First([1..Length(trees)], k-> e[1] in trees[k]);
          j:=First([1..Length(trees)], k-> e[2] in trees[k]);
          if i<>j then
            trees[i]:=Union(trees[i], trees[j]);
            trees[j]:=[];
            Add(needed,e);
          fi;
          if Length(needed)=nv-1 then
            break;
          fi;
        od;
        return needed;
    end;
 
    distance := function(a,b)
        local   k,  gcd,  i;

        k := Length(a);
        if k <> Length(b) then
            Error("The lengths of a and b are different.\n");
        fi;


        gcd := [];
        for i in [1..k] do
            Add(gcd, Minimum(a[i],b[i]));
        od;
        return(Maximum(Sum(a-gcd),Sum(b-gcd)));

    end;

    graph:=Graph(GraphType.UNDIRECTED);
    SetSimulation(graph,true);
    SetDrag(graph,true);
    SetShowNeighbours(graph,true);
    nf:=Length(f);
    fs:=[];
    for i in [1..nf] do 
        fs[i]:=Shape(ShapeType!.CIRCLE, Concatenation("(",JoinStringsWithSeparator(f[i],","),")"));
        SetLayer(fs[i],Sum(f[i]));
        SetSize(fs[i],1);
        Add(graph,fs[i]);
    od;
    c:=Cartesian([1..nf],[1..nf]);
    c:=Filtered(c,p->p[1]<p[2] and f[p[1]]*f[p[2]]<>0);
    Sort(c,function(e,ee) return distance(f[e[1]],f[e[2]])<distance(f[ee[1]],f[ee[2]]); end);
    tv:=Kruskal(f,List(c,p->[f[p[1]],f[p[2]]]));
    for p in c do 
        ln:=Link(fs[p[1]],fs[p[2]]);
        #SetWeight(ln, distance(f[p[1]],f[p[2]]));
        SetTitle(ln, String(distance(f[p[1]],f[p[2]])));
        if [f[p[1]],f[p[2]]] in tv then 
            SetColor(ln,"red");
        fi;
        Add(graph,ln);
    od;
    canvas:=Canvas("Factorizations graph");
    Add(canvas,graph);
    return Draw(canvas);
end);


############################################################################
##
#O DotEliahouGraph(f)
##
## f is a set of factorizations 
## Draws the Eliahou graph of factorizations associated to f: a graph 
## whose vertices are the elements of f, and there is an edge between
## two vertices if they have common support. Edges are labelled with
## distances between nodes they join.
##
#############################################################################

InstallMethod(DrawEliahouGraph, "draws Eliahou's graph for element in numerical semigroup", [IsRectangularTable], 
function(f)
    local graph, canvas, fs, c, nf, i, p;
    
    graph:=Graph(GraphType.UNDIRECTED);
    SetShowNeighbours(graph,true);
    SetSimulation(graph,true);
    SetDrag(graph,true);
    nf:=Length(f);
    fs:=[];
    for i in [1..nf] do 
        fs[i]:=Shape(ShapeType!.CIRCLE, Concatenation("(",JoinStringsWithSeparator(f[i],","),")"));
        SetLayer(fs[i],Sum(f[i]));
        SetSize(fs[i],1);
        Add(graph,fs[i]);
    od;
    c:=Cartesian([1..nf],[1..nf]);
    c:=Filtered(c,p->p[1]<p[2] and f[p[1]]*f[p[2]]<>0);
    for p in c do 
        Add(graph,Link(fs[p[1]],fs[p[2]]));
    od;
    canvas:=Canvas("Eliahou graph");
    Add(canvas,graph);
    return Draw(canvas);
end);

############################################################################
##
#O DrawRosalesGraph(n,s)
## s is either a numerical semigroup or an affine semigroup, and n is an
## element of s
## Draws the graph associated to n in s.
##
#############################################################################
InstallMethod(DrawRosalesGraph, "for numerical semigroups",  [IsInt,IsNumericalSemigroup], 
function(n,s)
    local graph, canvas, msg, msgs, c, nv, i, p;
    
    msg:=Filtered(MinimalGenerators(s), g->n-g in s);
    graph:=Graph(GraphType.UNDIRECTED);
    SetSimulation(graph,true);
    SetDrag(graph,true);
    nv:=Length(msg);
    msgs:=[];
    for i in [1..nv] do 
        msgs[i]:=Shape(ShapeType!.CIRCLE, String(msg[i]));
        SetSize(msgs[i],1);
        Add(graph,msgs[i]);
    od;
    c:=Cartesian([1..nv],[1..nv]);
    c:=Filtered(c,p->p[1]<p[2] and n-(msg[p[1]]+msg[p[2]]) in s);
    for p in c do 
        Add(graph,Link(msgs[p[1]],msgs[p[2]]));
    od;
    canvas:=Canvas("Rosales graph");
    Add(canvas,graph);
    return Draw(canvas);
end);

InstallMethod(DrawRosalesGraph, "for affine semigroups",  [IsHomogeneousList,IsAffineSemigroup], 
function(n,s)
    local graph, canvas, msg, msgs, c, nv, i, p;
    
    msg:=Filtered(MinimalGenerators(s), g->n-g in s);
    graph:=Graph(GraphType.UNDIRECTED);
    SetSimulation(graph,true);
    SetDrag(graph,true);
    nv:=Length(msg);
    msgs:=[];
    for i in [1..nv] do 
        msgs[i]:=Shape(ShapeType!.CIRCLE, String(msg[i]));
        SetSize(msgs[i],1);
        Add(graph,msgs[i]);
    od;
    c:=Cartesian([1..nv],[1..nv]);
    c:=Filtered(c,p->p[1]<p[2] and n-(msg[p[1]]+msg[p[2]]) in s);
    for p in c do 
        Add(graph,Link(msgs[p[1]],msgs[p[2]]));
    od;
    canvas:=Canvas("Rosales graph");
    Add(canvas,graph);
    return Draw(canvas);
end);


############################################################################
##
#F DrawOverSemigroupsNumericalSemigroup(s)
##  Draws the Hasse diagram of 
##  oversemigroupstree of the numerical semigroup s.
##  Irreducible numerical semigroups are highlighted.
##
############################################################################
InstallGlobalFunction(DrawOverSemigroupsNumericalSemigroup, function(s)
    local ov, graphHasse, canvas,c,i,r,ovs,n,hasse,lbl;
    
    hasse:=function(rel)
      local dom, out;
      dom:=Flat(rel);
      out:=Filtered(rel, p-> ForAny(dom, x->([p[1],x] in rel) and ([x,p[2]] in rel)));
      return Difference(rel,out);
    end;

    ov:=OverSemigroupsNumericalSemigroup(s);
    n:=Length(ov);
    graphHasse := Graph(GraphType.UNDIRECTED);
    SetSimulation(graphHasse,true);
    SetDrag(graphHasse,true);
    c:=Cartesian([1..n],[1..n]);
    c:=Filtered(c, p-> p[2]<>p[1]);
    c:=Filtered(c, p-> IsSubset(ov[p[1]],ov[p[2]]));
    c:=hasse(c);
    ovs:=[];
    for i in [1..n] do
    lbl:=Concatenation("$\\langle",JoinStringsWithSeparator(MinimalGenerators(ov[i]),","),"\\rangle$");


        if IsIrreducible(ov[i]) then
            ovs[i]:=Shape(ShapeType!.DIAMOND, lbl);
        else
            ovs[i]:=Shape(ShapeType!.CIRCLE, lbl);
        fi;
        SetLayer(ovs[i],Genus(ov[i]));
    SetSize(ovs[i],2);
    Add(graphHasse,ovs[i]);
    od;
    for r in c do
        Add(graphHasse,Link(ovs[r[1]],ovs[r[2]]));
    od;
    canvas:=Canvas("Oversemigroups");
    SetTexTypesetting(canvas, true);
    Add(canvas,graphHasse);
    return Draw(canvas);    
end);

############################################################################
##
#F DrawTreeOfSonsOfNumericalSemigroup(s, l, generators)
##  Draws the tree of sons of numerical semigroups up to level l with 
##  respect to the minimal system of generators given by the function
##  generators. 
##
############################################################################
InstallGlobalFunction(DrawTreeOfSonsOfNumericalSemigroup,
function(s,l,generators)
    local gens, frb, desc, graphTreee, d, shpr, shp, canvas, sonsf, lbl;


    sonsf:=function(s,n,lv)
        local gens, frb, desc, d, shp;
        if lv=0 then
            return ;
        fi;
        frb:=FrobeniusNumber(s);
        gens:=Filtered(generators(s), x-> x>frb);
        desc:=List(gens, g->RemoveMinimalGeneratorFromNumericalSemigroup(g,s));
        for d in desc do
            lbl:=Concatenation("$\\{",JoinStringsWithSeparator(generators(d),","),"\\}$");
            shp:=Shape(ShapeType!.CIRCLE, lbl);
            SetSize(shp,5);
            Add(graphTreee,shp);
            SetParentNode(shp,n);
            sonsf(d,shp,lv-1);
        od;
        if desc<>[] then
            return ;
        fi;
        #Add(canvas, FrancyMessage(FrancyMessageType.WARNING, "This semigroup is a leaf"));
        return ;
    end;

    frb:=FrobeniusNumber(s);
    gens:=Filtered(generators(s), x-> x>frb);
    desc:=List(gens, g->RemoveMinimalGeneratorFromNumericalSemigroup(g,s));

    graphTreee := Graph(GraphType.TREE);
    SetCollapsed(graphTreee,false);
    shpr:=Shape(ShapeType!.CIRCLE, "S");
    SetSize(shpr,5);
    Add(shpr,FrancyMessage(String(generators(s))));
    Add(graphTreee,shpr);
    canvas:=Canvas("Sons of a numerical semigroup");
    SetTexTypesetting(canvas, true);
    Add(canvas,graphTreee);
    sonsf(s,shpr,l);
    return Draw(canvas);
end);



########################################################################
##
#F DrawHasseDiagramOfNumericalSemigroup(s, A)
## plots a graph whose set of vertices is A, which is a finite set of 
## integers, and whose edges are provided by the order of the numerical 
## semigroup s.
##
#########################################################################
InstallGlobalFunction(DrawHasseDiagramOfNumericalSemigroup, 
function(s, A, t...)
    local hasse, order, layers, _layers, showfacts, graphHasse, rel, V, l, e, i, canvas, message, title, graph;    
    
    if not IsNumericalSemigroup(s) then
        Error("The argument must be a numerical semigroup.\n");
    fi;
    
    if Length(t) > 0 then
        title := t[1];
    else
        title := "Hasse diagram of numerical semigroup";
    fi;

    # rel is a set of lists with two elements representing a binary relation
    # hasse(rel) removes from rel the pairs [x,y] such that there exists
    # z with [x,z], [z,y] in rel
    hasse := function(rel)
        local dom, out;
        dom := Set(Flat(rel));
        out := Filtered(rel, p -> ForAny(dom, x -> ([p[1], x] in rel) and ([x, p[2]] in rel)));
        return Difference(rel, out);
    end;
    
    # determine the layer depth for every vertex of the graph
    layers := function(edges)
        local graph, leaves, depths;        
        graph := ListWithIdenticalEntries(Length(A), 0);
        Apply(graph, p -> []);
        depths := ListWithIdenticalEntries(Length(A), 0);
        leaves := ListWithIdenticalEntries(Length(A), true);
        for e in edges do
            Add(graph[e[1]], e[2]);
            leaves[e[2]] := false;            
        od;        
        
        for i in [1..Length(A)] do
            if leaves[i] = true then
                _layers(i, graph, depths);
            fi;
        od;        
                
        return depths;
    end;
    # Auxiliar function to compute the layers depth
    _layers := function(node, graph, depths)
        local son;        
        for son in graph[node] do
            depths[son] := Minimum(depths[son], depths[node]-1);
            _layers(son, graph, depths);            
        od;            
    end;
    
    # Function which shows the information regarding each vertex of the graph
    showfacts := function(x)
        if x in s then
            message := FrancyMessage(Concatenation(String(x), " factors as "), 
                                     String(FactorizationsElementWRTNumericalSemigroup(x, s)));
        else
            message := FrancyMessage(Concatenation(String(x)," is not an element of the numerical semigroup"));
        fi;
            
        SetId(message, Concatenation("message-for-", String(x)));
        Add(canvas, message);
        return Draw(canvas);
    end;
    
    # Initialize the graph
    graphHasse := Graph(GraphType.DIRECTED);
    SetSimulation(graphHasse, true);
    SetDrag(graphHasse, true);
    
    # Build the binary relation and apply the hasse function
    A := Set(A);
    rel := Cartesian([1..Length(A)],[1..Length(A)]);
    rel := Filtered(rel, p -> A[p[2]] <> A[p[1]]);
    rel := Filtered(rel, p -> A[p[2]] - A[p[1]] in s);
    rel := hasse(Set(rel));
    
    # Add the graph vertices
    V := [];
    l := layers(rel);    
    for i in [1..Length(A)] do
        V[i] := Shape(ShapeType!.CIRCLE, String(A[i]));
        SetLayer(V[i], l[i]);
        Add(V[i], Callback(showfacts, [A[i]]));
        if V[i] in s then
            message := String(FactorizationsElementWRTNumericalSemigroup(A[i],s));
            Add(V[i], FrancyMessage(message));
        fi;      
        Add(graphHasse, V[i]);
    od;
    
    # Add the graph vertices and edges
    for e in rel do
        Add(graphHasse, Link(V[e[1]], V[e[2]]));
    od;
    
    # Return the canvas
    canvas := Canvas(title);
    Add(canvas, graphHasse);
    return Draw(canvas);    
end);


########################################################################
##
#F DrawTreeOfGluingsOfNumericalSemigroup(s, expand...)
## Returns a Francy canvas with the tree of gluings of the numerical 
## semigroup s. If the optional argument expand is provided, then
## the tree is drawn fully expanded.
##
#########################################################################
InstallGlobalFunction(DrawTreeOfGluingsOfNumericalSemigroup, 
function(s, expand...)
    local SystemOfGeneratorsToString, rgluings, tree, canvas, root, expand_tree;

    SystemOfGeneratorsToString := function(sg)
        return Concatenation("〈", JoinStringsWithSeparator(sg, ","), "〉");
    end;
    
    if Length(expand) = 0 then
        expand_tree := false;
    else
        expand_tree := true;
    fi;

    # Recursively plot the gluings tree 
    rgluings := function(s, node)
        local lg, label, shape, son1, son2, gen1, gen2, p;
        
        # For each possible gluing plot the gluing and the numerical semigroups associated.
        lg := AsGluingOfNumericalSemigroups(s);        
        for p in lg do
            # Plot the gluing information
            label := Concatenation(SystemOfGeneratorsToString(p[1])," + ", SystemOfGeneratorsToString(p[2]));
            shape := Shape(ShapeType!.SQUARE, label);
            SetSize(shape,1);            
            Add(tree, shape);
            SetParentNode(shape, node);            
            
            # Plot the two numerical semigroups involved
            gen1 := p[1] / Gcd(p[1]);
            gen2 := p[2] / Gcd(p[2]);            
            son1 := Shape(ShapeType!.CIRCLE, SystemOfGeneratorsToString(gen1));
            son2 := Shape(ShapeType!.CIRCLE, SystemOfGeneratorsToString(gen2));
            SetSize(son1, 1);            
            SetSize(son2, 1);            
            
            Add(tree, son1);
            Add(tree, son2);
            SetParentNode(son1, shape);
            SetParentNode(son2, shape);
            
            rgluings(NumericalSemigroup(gen1), son1);            
            rgluings(NumericalSemigroup(gen2), son2);            
        od;
        
        return ;
    end;

    tree := Graph(GraphType.TREE);
    root := Shape(ShapeType!.CIRCLE, SystemOfGeneratorsToString(MinimalGenerators(s)));
    SetSize(root, 1);           
    Add(tree, root);

    if expand_tree then 	
	SetCollapsed(tree, false);
    fi;

    canvas := Canvas("Gluings of a numerical semigroup");
    Add(canvas, tree);
    rgluings(s, root);    
    return Draw(canvas);
end);
