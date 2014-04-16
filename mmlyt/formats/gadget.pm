#
#
#   key        class    group      dir      name                    rwfunc        mkfunc       size    flag     default                                              
#   str        str      str        str      str                     str           str          str     str      str                                                  
#####################################################################################################################################################################
    icgen      icgen    icgen      main     icgen                   None          None         small   dir      None
    stat_ic    icgen    icgen      icgen    static_gdic             rw_snap       mk_ic        small   shpfix   None
    stat_pm    icgen    icgen      icgen    static_gdpm             rw_param      mk_param     small   shpfix   /home/langmm/bin/Gadget/default.param
    makelog    icgen    icgen      icgen    static_gdmklog          rw_makelog    mk_makelog   small   shpfix   None
    input      input    input      main     input                   None          None         small   dir      None
    exec       input    exec       input    Gadget2                 None          mk_exec      small   None     /home/langmm/bin/Gadget/Gadget-2.0.6/Gadget2/Gadget2
    makefile   input    exec       input    Makefile                rw_makefile   mk_exec      small   None     /home/langmm/bin/Gadget/Gadget-2.0.6/Gadget2/Makefile
    ic         input    ic         input    ic_gadget               rw_snap       mk_ic        small   None     None
    param      input    input      input    gdpar                   rw_param      mk_param     small   None     /home/langmm/bin/Gadget/default.param
    pbs        input    input      input    gdpbs                   rw_pbs        mk_pbs       small   None     None
    resub      input    input      input    resub                   rw_sub        mk_sub       small   None     None
    outlist    input    input      input    output_times.txt        rw_outlist    mk_outlist   small   None     None
    output     output   output     main     output                  None          None         small   dir      None
    upar       output   output     output   parameters-usedvalues   rw_param      None         small   nopfix   None
    uparam     output   output     input    gdpar-usedvalues        rw_param      None         small   None     None
    runout     output   output     output   output                  None          None         small   None     None
    snapbase   output   snapshot   output   snapshot                rw_snap       None         large   None     None
    restbase   output   restart    output   restart                 None          None         large   None     None
    energy     output   output     output   energy                  rw_energy     None         small   None     None
    cpu        output   output     output   cpu                     rw_cpu        None         small   None     None
    info       output   output     output   info                    rw_info       None         small   None     None
    timings    output   output     output   timings                 rw_timings    None         small   None     None
