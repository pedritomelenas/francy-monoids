#######################################################################################
# Functions to draw the tree of the sons of a numerical semigroup. 
# The library Francis is used. See https://github.com/mcmartins/francy
#######################################################################################

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
