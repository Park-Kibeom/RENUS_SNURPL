! file name = thop.h - nominal t/h operation data
      real ::        powfa     , &  !assembly power in watt
                     powfa0    , &  !assembly power in Mw
                     powlin    , &  !linear power density
                     cfrperfa  , &  !coolant mass flow rate per assy in kg/sec
                     fracdc    , &  !fraction of heat deposited directly in coolant
                     fracdf    , &  !fraction of heat in fuel
                     tin       , &  !inlet temperature
                     din       , &  !inlet density
                     tdopin    , &  !inlet doppler temperature
                     hgap      , &  !gap conductance in w/m^2-C
                     tcnom     , &  !nominal coolant temperature
                     tfnom     , &  !nominal fuel temperature
                     toutavg   , &  !average outlet temperature
                     tfmax     , &  !maximum fuel centerline temperature
                     tmmax     , &  !maximum moderator temperature   ! 2014_10_06 . scb
                     fracdvb   , &  !fraction of heat deposited directly into vessel bypass
                     fracdwr   , &  !fraction of heat deposited directly into water rods
                     byp_dsat  , &  !saturation density of coolant
                     byp_a_frac, &  !bypass area fraction
                     wr_a_frac      !water rod area fraction
 

      common /thop/  powfa, powfa0, powlin, cfrperfa, fracdc, fracdf, tin, din, tdopin, hgap, &
                     tcnom, tfnom, toutavg, tfmax, tmmax, fracdvb, fracdwr, byp_dsat, byp_a_frac, wr_a_frac
 
