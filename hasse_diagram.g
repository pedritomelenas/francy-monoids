#######################################################################################
# Functions to draw hasse diagrams of the integers with the order of a 
# numerical semigroup. The library Francis is used.
# See https://github.com/mcmartins/francy
#######################################################################################

########################################################################
##
#F HasseDiagramOfNumericalSemigroup(s, A)
## plots a graph whose set of vertices is A, which is a finite set of 
## integers, and whose edges are provided by the order of the numerical 
## semigroup s.
##
#########################################################################
HasseDiagramOfNumericalSemigroup := function(s, A, t...)
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
    graphHasse := Graph(GraphType.HASSE);
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
    return canvas;    
end;

########################################################################
##
#F BettiGraphOfNumericalSemigroup(s)
## plots the graph of the Betti elements of the numerical semigroup s.
##
#########################################################################
BettiHasseDiagramOfNumericalSemigroup := function(s)
    if not IsNumericalSemigroup(s) then
        Error("The argument must be a numerical semigroup.\n");
    fi;
    
    return HasseDiagramOfNumericalSemigroup(s, BettiElementsOfNumericalSemigroup(s), "Betti elements");    
end;

########################################################################
##
#F AperyHasseOfNumericalSemigroup(args)
## plots the graph of an Apèry set of a numerical semigroup, Ap(s, n).
## This function takes two arguments:
## - args[1] : the numerical semigroup s.
## - args[2] : an element n of the numerical semigroup (optional). The 
##             multiplicty of s is used by default.
##
#########################################################################
AperyHasseDiagramOfNumericalSemigroup := function(arg)
    local s, n;
    
    if Length(arg) = 1 then
        s := arg[1];
        n := MultiplicityOfNumericalSemigroup(s);
    fi;
    if Length(arg) = 2 then
        s := arg[1];
        n := arg[2];
    fi;
    if Length(arg) > 2 then
        Error("The number of arguments must be one or two");
    fi;
    
    if not IsNumericalSemigroup(s) then
        Error("The argument must be a numerical semigroup.\n");
    fi;
    
    return HasseDiagramOfNumericalSemigroup(s, AperyListOfNumericalSemigroup(s, n), "Apery set");    
end;
