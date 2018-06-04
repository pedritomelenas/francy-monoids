##  this creates the documentation, needs: GAPDoc and AutoDoc packages, pdflatex
##  
##  Call this with GAP from within the package directory.
##
## Thanks Sebastian!

if fail = LoadPackage("AutoDoc", ">= 2016.01.21") then
    Error("AutoDoc 2016.01.21 or newer is required");
fi;

AutoDoc( rec( scaffold := true, autodoc := true ) );
PrintTo("version", GAPInfo.PackageInfoCurrent.Version);
