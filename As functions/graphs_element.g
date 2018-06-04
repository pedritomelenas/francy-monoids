#######################################################################################
# Functions to draw the graphs of an element of a numerical semigroup. 
# There are two different graphs. One concerns the factorizations of the element
# and it is due to Eliahou. The other one was published by Rosales and its vertices
# are the minimal generators of the semigroup.
# 
# The library Francis is used.
# See https://github.com/mcmartins/francy
#######################################################################################

EliahouGraph:=function(ns,s)
    local facts, graphs, graph, canvas, fs, f, c, n, i, p;
    
    facts:=List(ns, n->FactorizationsElementWRTNumericalSemigroup(n,s));
    graphs:=[];
    for f in facts do
        graph:=Graph(GraphType.UNDIRECTED);
        SetSimulation(graph,true);
        SetDrag(graph,true);
        n:=Length(f);
        fs:=[];
        for i in [1..n] do 
            fs[i]:=Shape(ShapeType!.CIRCLE, Concatenation("(",JoinStringsWithSeparator(f[i],","),")"));
            SetLayer(fs[i],Sum(f[i]));
            SetSize(fs[i],1);
            Add(graph,fs[i]);
        od;
        c:=Cartesian([1..n],[1..n]);
        c:=Filtered(c,p->p[1]<p[2] and f[p[1]]*f[p[2]]<>0);
        for p in c do 
            Add(graph,Link(fs[p[1]],fs[p[2]]));
        od;
        Add(graphs,graph);
    od;
    canvas:=Canvas("Eliahour graph");
    for graph in graphs do
        Add(canvas,graph);
    od;
    return Draw(canvas);
end;

DrawEliahouGraph:=function(n,s)
    local graph, canvas, f, fs, c, nf, i, p;
    
    f:=FactorizationsElementWRTNumericalSemigroup(n,s);
    graph:=Graph(GraphType.UNDIRECTED);
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
end;

DrawRosalesGraph:=function(n,s)
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
end;

DrawFactorizationGraph:=function(n,s)
    local graph, canvas, f, fs, c, nf, i, p, ln, distance, Kruskal, tv;

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

    f:=FactorizationsElementWRTNumericalSemigroup(n,s);
    graph:=Graph(GraphType.UNDIRECTED);
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
end;


DrawFactorizationGraph:=function(f)
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
end;
