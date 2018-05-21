#######################################################################################
# Functions to draw the gluings of a numerical semigroup with Francy.
# See https://github.com/mcmartins/francy
#######################################################################################


########################################################################
##
#F GluingsHasse(s)
## Returns a Francy canvas with a hasse diagram of the gluings of the 
## numerical semigroup s. When clicking on a numerical semigroup, its 
## gluings are computed.
##
#########################################################################
GluingsHasse := function(s)
    local SystemOfGeneratorsToString, pgluings, tree, canvas, root;
    
    # Return a string with the list sg between < and >
    SystemOfGeneratorsToString := function(sg)
        return Concatenation("〈", JoinStringsWithSeparator(sg, ","), "〉");
    end;
    
    # Draw the gluings of the semigroup s when the node is clicked 
    pgluings := function(s, layer, node)
        local lg, label, shape, son1, son2, gen1, gen2, p;
        
        # Each element of the list of gluings lg is drawn.
        lg := AsGluingOfNumericalSemigroups(s);        
        for p in lg do
            # Draw the gluing node
            label := Concatenation(SystemOfGeneratorsToString(p[1])," + ", SystemOfGeneratorsToString(p[2]));
            shape := Shape(ShapeType!.SQUARE, label);
            SetId(shape, label);
            SetLayer(shape, layer + 1);
            SetSize(shape,1);
            
            Add(tree, shape);
            Add(tree, Link(node, shape));
            
            # Draw the semigroups involved in the gluing
            gen1 := p[1] / Gcd(p[1]);
            gen2 := p[2] / Gcd(p[2]);            
            son1 := Shape(ShapeType!.CIRCLE, SystemOfGeneratorsToString(gen1));
            son2 := Shape(ShapeType!.CIRCLE, SystemOfGeneratorsToString(gen2));
            SetId(son1, Concatenation(SystemOfGeneratorsToString(gen1), String(layer+1)));
            SetId(son2, Concatenation(SystemOfGeneratorsToString(gen2), String(layer+1)));
            SetLayer(son1, layer + 2);
            SetLayer(son2, layer + 2);
            SetSize(son1, 1);            
            SetSize(son2, 1);            
            
            Add(tree, son1);
            Add(tree, son2);
            Add(tree, Link(shape, son1));            
            Add(tree, Link(shape, son2));
            
            # Add callback functionality to the semigroups' nodes.
            Add(son1, Callback(pgluings, [NumericalSemigroup(gen1), layer+2, son1]));
            Add(son2, Callback(pgluings, [NumericalSemigroup(gen2), layer+2, son2]));
        od;
        
        # Draw the new canvas
        if lg <> [] then
            return Draw(canvas);
        fi;
        Add(canvas, FrancyMessage(FrancyMessageType.WARNING, "This semigroup is not a gluing"));
        return Draw(canvas);
    end;
    
    # Build a hasse diagram with the numerical semigroup s.
    tree := Graph(GraphType.HASSE);
    root := Shape(ShapeType!.CIRCLE, SystemOfGeneratorsToString(MinimalGenerators(s)));
    SetLayer(root, 0);
    SetSize(root, 1);
    Add(tree, root);
    
    # Add callback functionality to the first node.
    Add(root, Callback(pgluings, [s, 0, root]));

    canvas := Canvas("Gluings of a numerical semigroup");
    Add(canvas, tree);
    return canvas;
end;

########################################################################
##
#F GluingsTree(s)
## Returns a Francy canvas with the tree of gluings of the numerical 
## semigroup s. 
##
#########################################################################
GluingsTree := function(s)
    local SystemOfGeneratorsToString, rgluings, tree, canvas, root;

    SystemOfGeneratorsToString := function(sg)
        return Concatenation("〈", JoinStringsWithSeparator(sg, ","), "〉");
    end;
    
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

    canvas := Canvas("Gluings of a numerical semigroup");
    Add(canvas, tree);
    rgluings(s, root);    
    return canvas;
end;
