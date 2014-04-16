#
#
#   key        class    group      dir      name                  alias             rwfunc        mkfunc       size    flag      default
#   str        str      str        str      str                   str               str           str          str     str       str   
##################################################################################################################################################################
    icgen      icgen    icgen      main     icgen                 None              None          None         small   dir       None
    stat_pm    icgen    icgen      icgen    static_scfpm          None              rw_param      mk_param     small   shpfix    /home/langmm/code/fortran/scf/scfpar
    stat_md    icgen    icgen      icgen    static_scfmd          None              rw_mods       mk_mods      small   shpfix    /home/langmm/code/fortran/scf/scfmod
    makelog    icgen    icgen      icgen    static_scfmklog       None              rw_makelog    mk_makelog   small   shpfix    None
    input      input    input      main     input                 None              None          None         small   dir       None
    exec       input    exec       input    mpiscf                mpiscf            None          mk_exec      small   None      /home/langmm/code/fortran/scf/mpiscf
    wrap       input    input      input    scfwrap               None              rw_wrapper    mk_wrapper   small   None      None
    param      input    input      input    scfpar                scfpar            rw_param      mk_param     small   None      /home/langmm/code/fortran/scf/scfpar
    mods       input    input      input    scfmod                scfmod            rw_mods       mk_mods      small   None      /home/langmm/code/fortran/scf/scfmod
    pbs        input    input      input    scfpbs                None              rw_pbs        mk_pbs       small   None      None
    sbatch     input    input      input    scfsbatch             None              rw_sbatch     mk_sbatch    small   None      None
    units      input    input      input    scfunits              None              rw_units      mk_units     small   None      None
    isnap      input    snapin     input    scfisnap{:03d}_{:s}   scfbi             rw_isnap      mk_snap      large   selfdir   None
    icoef      input    coefin     input    scficoef{:03d}_{:s}   scficoef          rw_icoef      None         large   selfdir   None
    output     output   output     main     output                None              None          None         small   dir       None
    runout     output   output     output   output                None              None          None         small   None      None
    procout    output   output     output   scfout{:03d}_{:s}     scfout            None          None         large   selfdir   None
    log        output   output     output   scflog{:03d}_{:s}     scflog            None          None         large   selfdir   None
    osnap      output   output     output   scfosnap{:03d}_{:s}   snap001           rw_osnap      None         large   selfdir   None
    ocoef      output   output     output   scfocoef{:03d}_{:s}   scfocoef          rw_ocoef      None         large   selfdir   None
    el         output   output     output   scfel{:03d}_{:s}      scfel             None          None         large   selfdir   None
    chk        output   output     output   scfchk{:03d}_{:s}     scfchkpt          None          None         large   selfdir   None
    olil       output   output     output   scfolil{:03d}_{:s}    slil              None          None         large   selfdir   None
