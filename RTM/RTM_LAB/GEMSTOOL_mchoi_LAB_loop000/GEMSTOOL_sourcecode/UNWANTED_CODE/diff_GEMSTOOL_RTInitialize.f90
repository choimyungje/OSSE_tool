6a7
> USE VLIDORT_LIN_IO_DEFS
52,54c53,65
<    VLIDORT_FixIn%Bool%TS_DO_SSCORR_TRUNCATION = .FALSE.
<    VLIDORT_FixIn%Bool%TS_DO_SSFULL            = .FALSE.
<    VLIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL    = .FALSE.
---
> !  Set MS-mode operations in VLIDORT, depenidng on first-order choices
> !    if using FO code, then VLIDORT operates only in Multiple-scatter mode.....
> 
>    if ( GEMSTOOL_INPUTS%RTcontrol%do_firstorder_option ) then
>       VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE     = .FALSE.
>       VLIDORT_FixIn%Bool%TS_DO_SS_EXTERNAL      = .TRUE.
>    else
>       VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE     = .TRUE.
>       VLIDORT_FixIn%Bool%TS_DO_SS_EXTERNAL      = .FALSE.
>    endif
> 
>    VLIDORT_FixIn%Bool%TS_DO_SSCORR_TRUNCATION   = .FALSE.
>    VLIDORT_FixIn%Bool%TS_DO_SSFULL              = .FALSE.
58a70,71
>    VLIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL      = .FALSE.
> 
62,64d74
<    VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE  = .TRUE.
< !      LAMBERTIAN_ALBEDO   = Will be set in Wavelength loop
< 
72a83,85
>    VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE  = .TRUE.
> !      LAMBERTIAN_ALBEDO   = Will be set in Wavelength loop
> 
85,95d97
< !  Set MS-mode operations in VLIDORT, depenidng on first-order choices
< !    if using FO code, then VLIDORT operates only in Multiple-scatter mode.....
< 
<    if ( GEMSTOOL_INPUTS%RTcontrol%do_firstorder_option ) then
<       VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE     = .false.
<       VLIDORT_FixIn%Bool%TS_DO_SS_EXTERNAL      = .true.
<    else
<       VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE     = .true.
<       VLIDORT_FixIn%Bool%TS_DO_SS_EXTERNAL      = .false.
<    endif
< 
106,134d107
<    VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES       = .TRUE.
<    
<    VLIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY = .FALSE.
<    VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION    = .TRUE.
< 
< !  Need to turn on Mean-valu (flux) output for Spherical Albedo calculation
< 
<    if ( GEMSTOOL_INPUTS%RTcontrol%do_SphericalAlbedo ) then
<       VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT    = .TRUE.
<    else
<       VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT    = .FALSE.
<    endif
< 
< !  Not doing Flux-alon calculation
< 
<    VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY          = .FALSE.
< 
< !  No thermal
<    
<    VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY   = .FALSE.
< 
< !  Alwasy need this flag
< 
<    VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES       = .TRUE. 
< 
< !  New flag (Observational Geometry). Default - Use Observational geometry.
< 
<    VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY = .TRUE.
< 
138,139c111,112
<       VLIDORT_ModIn%MBool%TS_DO_SSCORR_NADIR    = .false.
<       VLIDORT_ModIn%MBool%TS_DO_SSCORR_OUTGOING = .false.
---
>       VLIDORT_ModIn%MBool%TS_DO_SSCORR_NADIR    = .FALSE.
>       VLIDORT_ModIn%MBool%TS_DO_SSCORR_OUTGOING = .FALSE.
149a123,127
>    VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES       = .TRUE.
> 
>    VLIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY = .FALSE.
>    VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION    = .TRUE.
> 
164c142
<       VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST = .false. 
---
>       VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST = .false.
166,167c144,145
<       VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING = .false.  
<       VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING = .false. 
---
>       VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING = .false.
>       VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING = .false.
174c152
<       VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST = .true. 
---
>       VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST = .true.
176,177c154,155
<       VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING = .false.  
<       VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING = .false. 
---
>       VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING = .false.
>       VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING = .false.
179a158,181
> !  Alwasy need this flag
> 
>    VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES       = .TRUE.
> 
> !  Need to turn on Mean-valu (flux) output for Spherical Albedo calculation
> 
>    if ( GEMSTOOL_INPUTS%RTcontrol%do_SphericalAlbedo ) then
>       VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT    = .TRUE.
>    else
>       VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT    = .FALSE.
>    endif
> 
> !  Not doing Flux-alon calculation
> 
>    VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY          = .FALSE.
> 
> !  No thermal
> 
>    VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY   = .FALSE.
> 
> !  New flag (Observational Geometry). Default - Use Observational geometry.
> 
>    VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY = .TRUE.
> 
206a209,213
> !    -- not used
> 
>    VLIDORT_FixIn%Cont%TS_NLAYERS_NOMS     = 0
>    VLIDORT_FixIn%Cont%TS_NLAYERS_CUTOFF   = 0
> 
249a257,265
> !    -- One Level only, Level will be set later.
> 
>    VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS = 1           ! One level
>    VLIDORT_ModIn%MUserVal%TS_USER_LEVELS  = zero        ! Zero for now
> 
> !    -- Geometry specification height = Bottom of height grid (set later, zeroed here)
> 
>    VLIDORT_ModIn%MUserval%TS_GEOMETRY_SPECHEIGHT = zero
> 
259,267d274
< !    -- One Level only, Level will be set later.
<    
<    VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS = 1           ! One level
<    VLIDORT_ModIn%MUserVal%TS_USER_LEVELS  = zero        ! Zero for now
< 
< !    -- Geometry specification height = Bottom of height grid (set later, zeroed here)
< 
<    VLIDORT_ModIn%MUserval%TS_GEOMETRY_SPECHEIGHT = zero
< 
295d301
<    VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO     = zero
299,300c305,315
<    VLIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT     = zero
<    VLIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT     = zero
---
>    VLIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT      = zero
>    VLIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT      = zero
> 
> !    -- This is set later, zeroed here
> 
>    VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO     = zero
> 
> !    -- not used
> 
>    VLIDORT_FixIn%Optical%TS_LTE_DELTAU_VERT_INPUT = zero
>    VLIDORT_FixIn%Optical%TS_LTE_THERMAL_BB_INPUT  = zero
310c325
<       
---
> 
316c331
<       
---
> 
327,328c342,343
< !   VLIDORT_LinFixIn%Cont%TS_DO_LTE_LINEARIZATION    = .false.
< !   VLIDORT_LinFixIn%Cont%TS_DO_SURFBB_LINEARIZATION = .false.
---
>    !VLIDORT_LinFixIn%Cont%TS_DO_SURFBB_LINEARIZATION = .false.
>    !VLIDORT_LinFixIn%Cont%TS_DO_LTE_LINEARIZATION    = .false.
330,340c345,358
< !   VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS  = 0
< !   VLIDORT_LinFixIn%Cont%TS_COLUMNWF_NAMES     = ' '   
<    
< !  These quantities are se, just initialized here
<    
< !   VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG    = .false.
< !   VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER  = 0
< !   VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS = 0
< !   VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS      = 0
< !   VLIDORT_LinFixIn%Cont%TS_DO_SIMULATION_ONLY = .false.
< !   VLIDORT_LinFixIn%Cont%TS_PROFILEWF_NAMES    = ' '
---
>    !VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS  = 0
>    !VLIDORT_LinFixIn%Cont%TS_COLUMNWF_NAMES     = ' '
> 
>    !VLIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS     = .false.
>    !VLIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS       = 0
> 
> !  These quantities are just initialized here
> 
>    !VLIDORT_LinFixIn%Cont%TS_DO_SIMULATION_ONLY = .false.
>    !VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG    = .false.
>    !VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER  = 0
>    !VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS = 0
>    !VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS      = 0
>    !VLIDORT_LinFixIn%Cont%TS_PROFILEWF_NAMES    = ' '
345c363
< !  These things are set here 
---
> !  These things are set here
397,401c415,419
< !   VLIDORT_LinSup%BRDF%TS_LS_EXACTDB_BRDFUNC = ZERO
< !   VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0        = ZERO
< !   VLIDORT_LinSup%BRDF%TS_LS_BRDF_F          = ZERO
< !   VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0   = ZERO
< !   VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F     = ZERO
---
>    !VLIDORT_LinSup%BRDF%TS_LS_EXACTDB_BRDFUNC = ZERO
>    !VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0        = ZERO
>    !VLIDORT_LinSup%BRDF%TS_LS_BRDF_F          = ZERO
>    !VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0   = ZERO
>    !VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F     = ZERO
413,414c431,432
< !   VLIDORT_LinSup%BRDF%TS_LS_EMISSIVITY      = ZERO
< !   VLIDORT_LinSup%BRDF%TS_LS_USER_EMISSIVITY = ZERO
---
>    !VLIDORT_LinSup%BRDF%TS_LS_EMISSIVITY      = ZERO
>    !VLIDORT_LinSup%BRDF%TS_LS_USER_EMISSIVITY = ZERO
