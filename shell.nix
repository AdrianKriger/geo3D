with import <nixpkgs> { };
let
  # For packages pinned to a specific version
  pinnedHash = "nixos-23.11";
  pinnedPkgs = import (fetchTarball "https://github.com/NixOS/nixpkgs/archive/${pinnedHash}.tar.gz") { };
  pythonPackages = python3Packages;
in pkgs.mkShell rec {
  name = "impurePythonEnv";
  venvDir = "./.venv";
  buildInputs = [
    # A Python interpreter including the 'venv' module is required to bootstrap
    # the environment.
    pythonPackages.python

    # This executes some shell code to initialize a venv in $venvDir before
    # dropping into the shell
    pythonPackages.venvShellHook
    pinnedPkgs.virtualenv
    # Those are dependencies that we would like to use from nixpkgs, which will
    # add them to PYTHONPATH and thus make them accessible from within the venv.
    pythonPackages.numpy
    pinnedPkgs.git
    pythonPackages.pip
    pythonPackages.numpy
    pythonPackages.overpy
    pythonPackages.osmpythontools
    pythonPackages.setuptools
    pythonPackages.gdal
    pythonPackages.geopandas
    pythonPackages.mapbox-earcut
    pythonPackages.pydeck
    pythonPackages.overpy
    pythonPackages.jupyter
    pythonPackages.ipython
    pythonPackages.ipympl
    pinnedPkgs.vim
    pinnedPkgs.git
    pinnedPkgs.wget
    pinnedPkgs.gotop
    # For printing from jupyter
    # The list after scheme-small and latex are all .sty latex
    # modules that are needed for jupyter printing to work. 
    # I obtained the list from this issue:
    # https://github.com/jupyter/nbconvert/issues/1328#issue-659661022
    # Scheme-small is a small footprint latext install. The
    # latex schemes and the sytax for the entry below are 
    # described here:
    # https://nixos.wiki/wiki/TexLive
    # To actually generate a pdf in jupyter, do
    # File -> Save and export notebook as -> PDF
    (pinnedPkgs.texlive.combine { inherit (texlive)
        scheme-small latex adjustbox caption collectbox enumitem environ eurosym jknapltx
        parskip pgf rsfs tcolorbox titling trimspaces ucs ulem upquote 
        lastpage titlesec advdate pdfcol soul
        collection-langgerman collection-langenglish
    ;})
  ];
  # Run this command, only after creating the virtual environment
  PROJECT_ROOT = builtins.getEnv "PWD";

  postVenvCreation = ''
    unset SOURCE_DATE_EPOCH
    pip install -r requirements.txt
    echo "-----------------------"
    echo "🌈 Your Dev Environment is prepared."
    echo "Run ./impute.py"
    echo ""
    echo "📒 Note:"
    echo "-----------------------"
    echo "start vscode like this:"
    echo ""
    echo "./vscode.sh"
    echo "-----------------------"
  '';

  # Now we can execute any commands within the virtual environment.
  # This is optional and can be left out to run pip manually.
  postShellHook = ''
    # allow pip to install wheels
    unset SOURCE_DATE_EPOCH
  '';


}
