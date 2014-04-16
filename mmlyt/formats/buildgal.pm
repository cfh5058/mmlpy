#
#
#   key        class    group      dir      name                    rwfunc        mkfunc       size    flag      default                                              
#   str        str      str        str      str                     str           str          str     str       str                                                  
######################################################################################################################################################################
    icgen      icgen    icgen      main     icgen                   None          None         small   dir       None
    stat_ic    icgen    icgen      icgen    static_bgic             rw_snap       mk_ic        small   shpfix    None
    stat_pm    icgen    icgen      icgen    static_bgpm             rw_param      mk_param     small   shpfix    /home/langmm/bin/buildgalquiet/buildgal.in
    makelog    icgen    icgen      icgen    static_bgmklog          rw_makelog    mk_makelog   small   shpfix    None
    input      input    input      main     input                   None          None         small   dir       None
    exec       input    exec       input    buildgal                None          mk_exec      small   None      /home/langmm/bin/buildgalquiet/buildgal.stall
    makefile   input    exec       input    bgMakefile              rw_makefile   mk_exec      small   None      /home/langmm/bin/buildgalquiet/buildgal.compfile
    wrap       input    input      input    bgwrap                  rw_wrapper    mk_wrapper   small   None      None
    param      input    input      input    bgpar                   rw_param      mk_param     large   selfdir   /home/langmm/bin/buildgalquiet/buildgal.in
    pbs        input    input      input    bgpbs                   rw_pbs        mk_pbs       small   None      None
    units      input    input      input    bgunits                 rw_units      mk_units     small   None      None
    output     output   output     main     output                  None          None         small   dir       None
    runout     output   output     output   output                  None          None         small   None      None
    procout    output   output     output   bgout                   None          None         large   selfdir   None
    ic         output   ic         output   ic_buildgal             rw_snap       None         large   selfdir   None
    DDAT       output   output     output   DDAT                    None          None         large   selfdir   None
    STAT       output   output     output   STAT                    None          None         large   selfdir   None
    TREEBI     output   output     output   TREEBI                  None          None         large   selfdir   None
    rotc       output   output     output   rotc                    None          None         large   selfdir   None
    sigr       output   output     output   sigr                    None          None         large   selfdir   None
    sigt       output   output     output   sigt                    None          None         large   selfdir   None
    sigz       output   output     output   sigz                    None          None         large   selfdir   None
    toom       output   output     output   toom                    None          None         large   selfdir   None
    scflog     output   output     output   scflog                  None          None         large   selfdir   None
    ticklog    output   output     output   ticklog                 None          None         large   selfdir   None
