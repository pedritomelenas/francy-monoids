## TODO
# Integrate Arf functions
# Graphs of elements; add catenary "tree"
# Tree of complete intersections


########################################################################
##
#F apery(arg)
## plots the hasse diagram (Apery(s), <=_s).
##
#########################################################################
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
    c:=Filtered(c, p-> ap[p[2]] <> ap[p[1]]);
    c:=Filtered(c, p-> ap[p[1]] - ap[p[2]] in s);
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
