## TODO
# Integrate Arf functions
# Graphs of elements; add catenary "tree"
# Tree of complete intersections


apery:=function(arg)
    local ap,c,hasse, s, n, r, graphHasse, aps, es, canvas, i, order, showfacts, message;
    # rel is a list of lists with two elements representin a binary relation
    # hasse(rel) removes from rel the pairs [x,y] such that there exists
    # z with [x,z],[z,y] in rel
    hasse:=function(rel)
      local dom, out;
      dom:=Flat(rel);
      out:=Filtered(rel, p-> ForAny(dom, x->([p[1],x] in rel) and ([x,p[2]] in rel)));
      return Difference(rel,out);
    end;

    order:=function(x)
        return Maximum(LengthsOfFactorizationsElementWRTNumericalSemigroup(x,s));
    end;

    showfacts:=function(x)
        message := FrancyMessage(Concatenation(String(x), " factors as "), 
                    String(FactorizationsElementWRTNumericalSemigroup(x,s)));
        SetId(message, Concatenation("message-for-", String(x)));
        Add(canvas, message);
        return Draw(canvas);
    end;
    if Length(arg)=1 then
        s:=arg[1];
        n:=MultiplicityOfNumericalSemigroup(s);
    fi;
    if Length(arg)=2 then
        s:=arg[1];
        n:=arg[2];
    fi;
    if Length(arg)>2 then
        Error("The number of arguments must be one or two");
    fi;
  
    graphHasse := Graph(GraphType.HASSE);
    SetSimulation(graphHasse,true);
    SetDrag(graphHasse,true);
    ap:=AperyList(s,n);
    c:=Cartesian([1..n],[1..n]);
    c:=Filtered(c, p-> ap[p[2]]<>ap[p[1]]);
    c:=Filtered(c, p-> ap[p[1]]-ap[p[2]] in s);
    c:=hasse(c);
    aps:=[];
    for i in [1..n] do
        aps[i]:=Shape(ShapeType!.CIRCLE, String(ap[i]));
        SetLayer(aps[i],-order(ap[i]));
        Add(aps[i],Callback(showfacts,[ap[i]]));
        Add(aps[i],FrancyMessage(String(FactorizationsElementWRTNumericalSemigroup(ap[i],s))));
        Add(graphHasse,aps[i]);
    od;
    for r in c do
        Add(graphHasse,Link(aps[r[1]],aps[r[2]]));
    od;
    canvas:=Canvas("Apery");
    Add(canvas,graphHasse);
    return Draw(canvas);    
end;

sons:=function(s)
    local gens, frb, desc, graphHasse, d, shpr, shp, canvas, sonsf, i, gn, lbl;


    sonsf:=function(s,n)
        local gens, frb, desc, d, shp, i, lbl, gn;

        frb:=FrobeniusNumber(s);
        gens:=Filtered(MinimalGenerators(s), x-> x>frb);
        desc:=List(gens, g->RemoveMinimalGeneratorFromNumericalSemigroup(g,s));
        gn:=Genus(s);
        i:=0;
        for d in desc do
            i:=i+1;
            lbl:=Concatenation("S",String(gn),"-",String(i));
            shp:=Shape(ShapeType!.CIRCLE, lbl);
            shp!.layer:=Genus(d);
            Add(shp,Callback(sonsf,[d,shp]));
            Add(shp,HintMessage(String(MinimalGenerators(d))));
            Add(graphHasse,shp);
            Add(graphHasse,Link(n,shp));
            Print("adding new nodes");
        od;
        if desc<>[] then
            return Draw(canvas);
        fi;
    end;

    frb:=FrobeniusNumber(s);
    gens:=Filtered(MinimalGenerators(s), x-> x>frb);
    desc:=List(gens, g->RemoveMinimalGeneratorFromNumericalSemigroup(g,s));
    gn:=Genus(s);

    graphHasse := Graph(GraphType.HASSE);
    shpr:=Shape(ShapeType!.CIRCLE, String(Genus(s)));
    Add(shpr,HintMessage(String(MinimalGenerators(s))));
    shpr!.layer:=Genus(s);
    Add(graphHasse,shpr);
    i:=0;
    for d in desc do
        i:=i+1;
        lbl:=Concatenation("S",String(gn),"-",String(i));
        shp:=Shape(ShapeType!.CIRCLE, lbl);
        shp!.layer:=Genus(d);
        Add(shp,Callback(sonsf,[d,shp]));
        Add(shp,HintMessage(String(MinimalGenerators(d))));
        Add(graphHasse,shp);
        Add(graphHasse,Link(shpr,shp));
    od;
	canvas:=Canvas("Sons of a numerical semigroup");
    Add(canvas,graphHasse);
    return Draw(canvas);
end;



sonshasse:=function(s)
    local gens, frb, desc, graphTree, d, shpr, shp, canvas, sonsf, i, gn, lbl;


    sonsf:=function(s,n)
        local gens, frb, desc, d, shp, i, lbl, gn;

        frb:=FrobeniusNumber(s);
        gens:=Filtered(MinimalGenerators(s), x-> x>frb);
        desc:=List(gens, g->RemoveMinimalGeneratorFromNumericalSemigroup(g,s));
        i:=0;
        for d in desc do
            i:=i+1;
        lbl:=Concatenation(lb,":",String(i));
            gn:=String(MinimalGenerators(d));
            shp:=Shape(ShapeType!.CIRCLE, gn);
            Add(shp,Callback(sonsf,[d,shp]));
            Add(graphTree,shp);
            Add(n,shp);
        od;
        if desc<>[] then
            return Draw(canvas);
        fi;
        Add(canvas, HintMessage(MessageType.WARNING, "This semigroup is a leaf"));
        return Draw(canvas);
    end;

    frb:=FrobeniusNumber(s);
    gens:=Filtered(MinimalGenerators(s), x-> x>frb);
    desc:=List(gens, g->RemoveMinimalGeneratorFromNumericalSemigroup(g,s));
    gn:=Genus(s);

    graphTree := Graph(GraphType.TREE);
    graphTree!.collapse := false;
    shpr:=Shape(ShapeType!.CIRCLE, String(MinimalGenerators(s)));
    Add(graphTree,shpr);
    i:=0;
    for d in desc do
        i:=i+1;
        shp:=Shape(ShapeType!.CIRCLE, String(MinimalGenerators(d)));
        Add(shp,Callback(sonsf,[d,shp]));
        Add(graphTree,shp);
        Add(shpr,shp);
    od;
    canvas:=Canvas("Sons of a numerical semigroup");
    Add(canvas,graphTree);
    return Draw(canvas);
end;

sonstree:=function(s,l,generators)
    local gens, frb, desc, graphTreee, d, shpr, shp, canvas, sonsf;


    sonsf:=function(s,n,lv)
        local gens, frb, desc, d, shp;
        if lv=0 then
            return ;
        fi;
        frb:=FrobeniusNumber(s);
        gens:=Filtered(generators(s), x-> x>frb);
        desc:=List(gens, g->RemoveMinimalGeneratorFromNumericalSemigroup(g,s));
        for d in desc do
        shp:=Shape(ShapeType!.CIRCLE, String(generators(d)));
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
    Add(canvas,graphTreee);
    sonsf(s,shpr,l);
    return Draw(canvas);
end;


oversemigroups:=function(s)
    local ov, graphHasse, canvas,c,i,r,ovs,n,hasse;
    
    hasse:=function(rel)
      local dom, out;
      dom:=Flat(rel);
      out:=Filtered(rel, p-> ForAny(dom, x->([p[1],x] in rel) and ([x,p[2]] in rel)));
      return Difference(rel,out);
    end;

    ov:=OverSemigroupsNumericalSemigroup(s);
    n:=Length(ov);
    graphHasse := Graph(GraphType.HASSE);
    SetSimulation(graphHasse,true);
    SetDrag(graphHasse,true);
    c:=Cartesian([1..n],[1..n]);
    c:=Filtered(c, p-> p[2]<>p[1]);
    c:=Filtered(c, p-> IsSubset(ov[p[1]],ov[p[2]]));
    c:=hasse(c);
    ovs:=[];
    for i in [1..n] do
        if IsIrreducible(ov[i]) then
            ovs[i]:=Shape(ShapeType!.DIAMOND, String(MinimalGenerators(ov[i])));
        else
            ovs[i]:=Shape(ShapeType!.CIRCLE, String(MinimalGenerators(ov[i])));
        fi;
        SetLayer(ovs[i],Genus(ov[i]));
    SetSize(ovs[i],2);
    Add(graphHasse,ovs[i]);
    od;
    for r in c do
        Add(graphHasse,Link(ovs[r[1]],ovs[r[2]]));
    od;
    canvas:=Canvas("Oversemigroups");
    Add(canvas,graphHasse);
    return Draw(canvas);    
end;


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