#
#
#   key        class    group      dir      name              alias             rwfunc        mkfunc       size    flag      default                                              
#   str        str      str        str      str               str               str           str          str     str       str                                                  
#################################################################################################################################################################################
    icgen      icgen    icgen      main     icgen             None              None          None         small   dir       None
    stat_ic    icgen    icgen      icgen    static_bgic       None              rw_snap       mk_ic        small   shpfix    None
    stat_mg    icgen    icgen      icgen    static_bgmg       None              rw_snap       merge_ic     small   shpfix    None
    stat_pm    icgen    icgen      icgen    static_bgpm       None              rw_param      mk_param     small   shpfix    /home/langmm/code/fortran/buildgal/buildgal.in
    makelog    icgen    icgen      icgen    static_bgmklog    None              rw_makelog    mk_makelog   small   shpfix    None
    input      input    input      main     input             None              None          None         small   dir       None
    exec       input    exec       input    buildgal          buildgal.stall    None          mk_exec      small   None      /home/langmm/code/fortran/buildgal/buildgal.stall
    makefile   input    exec       input    bgMakefile        None              rw_makefile   mk_exec      small   None      /home/langmm/code/fortran/buildgal/buildgal.compile
    wrap       input    input      input    bgwrap            None              rw_wrapper    mk_wrapper   small   None      None
    param      input    longin     input    bgpar{:d}         buildgal.in       rw_param      mk_param     large   selfdir   /home/langmm/code/fortran/buildgal/buildgal.in
    pbs        input    input      input    bgpbs             None              rw_pbs        mk_pbs       small   None      None
    sbatch     input    input      input    bgsbatch          None              rw_sbatch     mk_sbatch    small   None      None
    units      input    input      input    bgunits           None              rw_units      mk_units     small   None      None
    output     output   output     main     output            None              None          None         small   dir       None
    runout     output   output     output   output            None              None          None         small   None      None
    procout    output   longout    output   bgout{:d}         bgout             None          None         large   selfdir   None
    ic         output   ic         output   ic_buildgal{:d}   ic_buildgal.out   rw_snap       None         large   selfdir   None
    DDAT       output   longout    output   DDAT{:d}          DDATA             None          None         large   selfdir   None
    STAT       output   longout    output   STAT{:d}          STATS             None          None         large   selfdir   None
    TREEBI     output   longout    output   TREEBI{:d}        TREEBI            None          None         large   selfdir   None
    rotc       output   longout    output   rotc{:d}          rotcurve.dat      None          None         large   selfdir   None
    sigr       output   longout    output   sigr{:d}          sigratio.dat      None          None         large   selfdir   None
    sigt       output   longout    output   sigt{:d}          sigt.dat          None          None         large   selfdir   None
    sigz       output   longout    output   sigz{:d}          sigzero.dat       None          None         large   selfdir   None
    toom       output   longout    output   toom{:d}          toomreqave.dat    None          None         large   selfdir   None
    scflog     output   longout    output   scflog{:d}        scflog            None          None         large   selfdir   None
    ticklog    output   longout    output   ticklog{:d}       tickerlog         None          None         large   selfdir   None
