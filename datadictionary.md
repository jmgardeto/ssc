<a name="toc"></a>
- [Parametros](#parameters)
  * [calc_design_pipe_vals](#calc_design_pipe_vals)
  * [custom_sf_pipe_sizes](#custom_sf_pipe_sizes)
  * [custom_sgs_pipe_sizes](#custom_sgs_pipe_sizes)
  * [custom_tes_p_loss](#custom_tes_p_loss)
  * [DP_SGS](#dp_sgs)
  * [has_hot_tank_bypass](#has_hot_tank_bypass)
  * [k_tes_loss_coeffs](#k_tes_loss_coeffs)
  * [L_rnr_pb](#l_rnr_pb)
  * [L_rnr_per_xpan](#l_rnr_per_xpan)
  * [L_xpan_hdr](#l_xpan_hdr)
  * [L_xpan_rnr](#l_xpan_rnr)
  * [Min_rnr_xpans](#min_rnr_xpans)
  * [N_hdr_per_xpan](#n_hdr_per_xpan)
  * [N_max_hdr_diams](#N_max_hdr_diams)
  * [northsouth_field_sep](#northsouth_field_sep)
  * [offset_xpan_hdr](#offset_xpan_hdr)
  * [sf_hdr_diams](#sf_hdr_diams)
  * [sf_hdr_lengths](#sf_hdr_lengths)
  * [sf_hdr_wallthicks](#sf_hdr_wallthicks)
  * [sf_rnr_diams](#sf_rnr_diams)
  * [sf_rnr_lengths](#sf_rnr_lengths)
  * [sf_rnr_wallthicks](#sf_rnr_wallthicks)
  * [sgs_diams](#sgs_diams)
  * [sgs_lengths](#sgs_lengths)
  * [sgs_wallthicks](#sgs_wallthicks)
  * [tanks_in_parallel](#tanks_in_parallel)
  * [T_tank_hot_inlet_min](#t_tank_hot_inlet_min)
  * [V_hdr_cold_max](#v_hdr_cold_max)
  * [V_hdr_cold_min](#v_hdr_cold_min)
  * [V_hdr_hot_max](#v_hdr_hot_max)
  * [V_hdr_hot_min](#v_hdr_hot_min)
  * [V_tes_des](#v_tes_des)
- [Outputs](#outputs)
  * [pipe_header_diams](#pipe_header_diams)
  * [pipe_header_expansions](#pipe_header_expansions)
  * [pipe_header_lengths](#pipe_header_lengths)
  * [pipe_header_mdot_dsn](#pipe_header_mdot_dsn)
  * [pipe_header_P_dsn](#pipe_header_P_dsn)
  * [pipe_header_T_dsn](#pipe_header_T_dsn)
  * [pipe_header_vel_dsn](#pipe_header_vel_dsn)
  * [pipe_header_wallthk](#pipe_header_wallthk)
  * [pipe_loop_P_dsn](#pipe_loop_P_dsn)
  * [pipe_loop_T_dsn](#pipe_loop_T_dsn)
  * [pipe_runner_diams](#pipe_runner_diams)
  * [pipe_runner_expansions](#pipe_runner_expansions)
  * [pipe_runner_lengths](#pipe_runner_lengths)
  * [pipe_runner_mdot_dsn](#pipe_runner_mdot_dsn)
  * [pipe_runner_P_dsn](#pipe_runner_P_dsn)
  * [pipe_runner_T_dsn](#pipe_runner_T_dsn)
  * [pipe_runner_vel_dsn](#pipe_runner_vel_dsn)
  * [pipe_runner_wallthk](#pipe_runner_wallthk)
  * [pipe_sgs_diams](#pipe_sgs_diams)
  * [pipe_sgs_mdot_dsn](#pipe_sgs_mdot_dsn)
  * [pipe_sgs_P_dsn](#pipe_sgs_p_dsn)
  * [pipe_sgs_T_dsn](#pipe_sgs_t_dsn)
  * [pipe_sgs_vel_dsn](#pipe_sgs_vel_dsn)
  * [pipe_sgs_wallthk](#pipe_sgs_wallthk)
 
 
<!-- toc -->

## Parametros
### calc_design_pipe_vals
true si las temperaturas y presiones htf en condiciones de diseño en los canales, el cabezal más lejano y el bucle más lejano deben calcularse y emitirse. Predeterminado = true. [^](#toc)

### custom_sf_pipe_sizes
true si deben usarse los parámetros de diámetros de canal y cabezal, espesores de pared y longitudes en lugar de calcularlos. Tenga en cuenta que cambiar las longitudes no afecta el diseño del campo. [^](#toc)

### custom_sgs_pipe_sizes
true si los parámetros de diámetros y espesores de pared de SGS deben usarse en lugar de calcularlos. (Tenga en cuenta que las longitudes de SGS siempre son entradas). [^](#toc)

### custom_tes_p_loss
true si las pérdidas de la tubería TES deben calcularse utilizando las longitudes de las tuberías TES y los coeficientes de pérdida menor (k_tes_loss_coeffs) o false si se utilizan los parámetros de potencia de bombeo en la página de parásitos. Predeterminado = false. [^](#toc)

### DP_SGS
la caída de presión en bar dentro del sistema generador de vapor (SGS) Predeterminado = 0. [^](#toc)

### has_hot_tank_bypass
true si la válvula de derivación del campo solar hace que el campo htf evite solo el tanque caliente (y el bloque de alimentación y la caldera auxiliar) y entre en el tanque frío antes de regresar al campo. El valor es false si la válvula de derivación pasa por alto tanto el tanque frío como el caliente. Predeterminado = false. [^](#toc)

### k_tes_loss_coeffs
coeficientes de pérdida menor combinados de los accesorios y válvulas en la colección (incluida la derivación) y los bucles de generación en la tubería TES que se utilizarán en la ecuación DP = K * U^2 * rho / 2. Predeterminado = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} [^](#toc)

### L_rnr_pb
longitud de tramo de tubería en metros, ya sea para las líneas frías o calientes. Esta longitud se compartió anteriormente con el otro conjunto idéntico de tramos  para la otra mitad del campo solar, pero este ya no es el caso. Por defecto = 25 m. [^](#toc)

### L_rnr_per_xpan
longitud umbral de tramo de tubería recto sin bucle de expansión. Una vez que se ha alcanzado esta longitud, se agrega un bucle de expansión (sin aumentar la distancia lineal). Por defecto = 70 m. [^](#toc)

### L_xpan_hdr
longitud combinada en metros de los dos segmentos perpendiculares de un bucle de expansión del cabezal. Esta es la longitud de tubería adicional para cada bucle de expansión. Predeterminado = 20 m [^](#toc)

### L_xpan_rnr
longitud combinada en metros de los dos segmentos perpendiculares de un bucle de expansión de tramo. Esta es la longitud de tubería adicional para cada bucle de expansión. Predeterminado = 20 m [^](#toc)

### Min_rnr_xpans
número mínimo de bucles de expansión por sección de tramo de un solo diámetro. Predeterminado = 1 [^](#toc)

### N_hdr_per_xpan
número de bucles de colector por bucles de expansión de cabecera. Predeterminado = 2. Valor = 1 significa que hay bucles de expansión entre cada bucle de colector. [^](#toc)

### N_max_hdr_diams
número máximo de diámetros permitidos en cada uno de los cabezales fríos y calientes. El número máximo de diámetros en los cabezales fríos y calientes es 2*N_max_hdr_diams. Predeterminado = 10. [^](#toc)

### northsouth_field_sep
separación norte / sur entre subcampos, en metros, definida como la distancia más corta entre SCA en los diferentes campos. Por defecto = 20 m. Valor = 0 significa que los SCA se están tocando. [^](#toc)

### offset_xpan_hdr
ubicación del primer bucle de expansión del cabezal. Predeterminado = 1, lo que significa que el primer bucle de expansión está después del primer bucle colector más cercano al tramo. [^](#toc)

### sf_hdr_diams
custom inner diameters for the header piping as read from the modified output files. Utilized if custom_sf_pipe_sizes is true. Do not change the number of values (sections) as this will result in unpredictable model behavior. [^](#toc)
 
### sf_hdr_lengths
custom lengths for the header piping as read from the modified output files. Utilized if custom_sf_pipe_sizes is true. Changing the lengths does not affect the field layout. Do not change the number of values (sections) as this will result in unpredictable model behavior. [^](#toc)
 
### sf_hdr_wallthicks
custom wall thicknesses for the header piping as read from the modified output files. Utilized if custom_sf_pipe_sizes is true. Do not change the number of values (sections) as this will result in unpredictable model behavior. [^](#toc)
 
### sf_rnr_diams
custom inner diameters for the runner piping as read from the modified output files. Utilized if custom_sf_pipe_sizes is true. Do not change the number of values (sections) as this will result in unpredictable model behavior. [^](#toc)
 
### sf_rnr_lengths
custom lengths for the runner piping as read from the modified output files. Utilized if custom_sf_pipe_sizes is true. Changing the lengths does not affect the field layout. Do not change the number of values (sections) as this will result in unpredictable model behavior. [^](#toc)
 
### sf_rnr_wallthicks
custom wall thicknesses for the runner piping as read from the modified output files. Utilized if custom_sf_pipe_sizes is true. Do not change the number of values (sections) as this will result in unpredictable model behavior. [^](#toc)
 
### sgs_diams
custom inner diameters for the SGS piping as read from the modified output files. Utilized if custom_sgs_pipe_sizes is true. Do not change the number of values (sections) as this will result in unpredictable model behavior. [^](#toc)

Collection Sections:
- 0: &nbsp;&nbsp;&nbsp; Solar field (SF) pump suction header to individual SF pump inlet
- 1: &nbsp;&nbsp;&nbsp; Individual SF pump discharge to SF pump discharge header
- 2: &nbsp;&nbsp;&nbsp; SF pump discharge header to collection field section headers (i.e., runners)
- 3: &nbsp;&nbsp;&nbsp; Collector field section outlet headers (i.e., runners) to expansion vessel (indirect storage) or hot thermal storage tank (direct storage)
- 4: &nbsp;&nbsp;&nbsp; Bypass branch - Collector field section outlet headers (i.e., runners) to pump suction header (indirect) or cold thermal storage tank (direct)

Generation Sections:
- 5: &nbsp;&nbsp;&nbsp; SGS pump suction header to individual SGS pump inlet (applicable only for storage in series with SF)
- 6: &nbsp;&nbsp;&nbsp; Individual SGS pump discharge to SGS pump discharge header (only for series storage)
- 7: &nbsp;&nbsp;&nbsp; SGS pump discharge header to steam generator supply header (only for series storage)
- 8: &nbsp;&nbsp;&nbsp; Steam generator supply header to inter-steam generator piping
- 9: &nbsp;&nbsp;&nbsp; Inter-steam generator piping to steam generator outlet header
- 10: &nbsp;&nbsp;&nbsp; Steam generator outlet header to SF pump suction header (indirect) or cold thermal storage tank (direct)

### sgs_lengths
length of piping in the SGS collection flow loop followed by the generation flow loop [m]. These are not read from the modified output files. Defaults = {0, 90, 100, 120, 0, 0, 0, 0, 80, 120, 80}. Lengths at indices 0, 1, 5 and 6 are the summed lengths of the multiple individual pump sections. Do not change the number of values (sections) as this will result in unpredictable model behavior. [^](#toc)

### sgs_wallthicks
custom wall thicknesses for the SGS piping as read from the modified output files. Utilized if custom_sgs_pipe_sizes is true. Do not change the number of values (sections) as this will result in unpredictable model behavior. [^](#toc)

### tanks_in_parallel
true if the hot and cold storage tank branch is in parallel with the solar field (traditional case), or false if the tanks are in series with the solar field (only applicable for direct storage). Default = true. [^](#toc)

### T_tank_hot_inlet_min
the minimum field htf temperature that may enter the hot tank [C]. If below this temperature the bypass valve is opened and the field recirculates. Default = 400 C. [^](#toc)
				
### V_hdr_cold_max
maximum allowed velocity in the cold header at design conditions. This value can be exceeded if the minimum would also be exceeded, but only if this puts it less out of range. [^](#toc)

### V_hdr_cold_min
minimum allowed velocity in the cold header at design conditions. This value can be exceeded if the maximum would also be exceeded, but only if this puts it less out of range. [^](#toc)

### V_hdr_hot_max
maximum allowed velocity in the hot header at design conditions. This value can be exceeded if the minimum would also be exceeded, but only if this puts it less out of range. [^](#toc)

### V_hdr_hot_min
minimum allowed velocity in the hot header at design conditions. This value can be exceeded if the maximum would also be exceeded, but only if this puts it less out of range. [^](#toc)

### V_tes_des
design-point velocity for sizing the diameters of the TES piping [m/s]. Default = 1.85 m/s. [^](#toc)


## Outputs
### pipe_header_diams
inner diameters in meters of all of the header sections in the cold and hot headers in one subfield. The first diameter is that before the first set of loops in the cold header and the last diameter is that after the last set of loops in the hot header. [^](#toc)

### pipe_header_expansions
number of expansions or contractions in the given header section [^](#toc)

### pipe_header_lengths
lengths in meters of the all of the header sections, including the added lengths of any expansion loops. The first length is that before the first set of loops in the cold header and last length is that after the last set of loops in the hot header. [^](#toc)

### pipe_header_mdot_dsn
mass flow rate in kg/s of the heat transfer fluid in each header section at design conditions. The first value is in the section before the first set of loops in the cold header and the last value is in the section after the last set of loops in the hot header. The mass flow for the cold header sections is the same as that entering the section, and the mass flow for the hot header sections is the same as that leaving the section. [^](#toc)

### pipe_header_P_dsn
gauge pressure in bar of the heat transfer fluid entering each section of the farthest header at design conditions. The first value is for the section before the first set of loops in the cold header and the last value is for the section after the last set of loops in the hot header. [^](#toc)

### pipe_header_T_dsn
temperature in Celsius of the heat transfer fluid entering each section of the farthest header at design conditions. The first value is for the section before the first set of loops in the cold header and the last value is for the section after the last set of loops in the hot header. [^](#toc)

### pipe_header_vel_dsn
velocity in m/s of the heat transfer fluid in each header section at design conditions. The first value is in the section before the first set of loops in the cold header and the last value is in the section after the last set of loops in the hot header. The velocity for the cold header sections is the same as that entering the section, and the velocity for the hot header sections is the same as that leaving the section. [^](#toc)

### pipe_header_wallthk
wall thickness of header pipe sections in [m] [^](#toc)

### pipe_loop_P_dsn
gauge pressure in bar of the heat transfer fluid entering each node in the farthest loop at design conditions. The values correspond to: [^](#toc)
- 0: &nbsp;&nbsp;&nbsp; the inlet interconnect carrying twice the loop mass flow rate
- 1: &nbsp;&nbsp;&nbsp; the interconnect before the first SCA
- 2: &nbsp;&nbsp;&nbsp; the first SCA
- 3: &nbsp;&nbsp;&nbsp; the interconnect between the first and second SCA
- 4: &nbsp;&nbsp;&nbsp; the second SCA
- ...
- n-3: &nbsp;&nbsp;&nbsp; the last SCA
- n-2: &nbsp;&nbsp;&nbsp; the interconnect after the last SCA
- n-1: &nbsp;&nbsp;&nbsp; the outlet interconnect carrying twice the loop mass flow rate

### pipe_loop_T_dsn
temperature in Celsius of the heat transfer fluid entering each node in the farthest loop at design conditions. The values correspond to: [^](#toc)
- 0: &nbsp;&nbsp;&nbsp; the inlet interconnect carrying twice the loop mass flow rate
- 1: &nbsp;&nbsp;&nbsp; the interconnect before the first SCA
- 2: &nbsp;&nbsp;&nbsp; the first SCA
- 3: &nbsp;&nbsp;&nbsp; the interconnect between the first and second SCA
- 4: &nbsp;&nbsp;&nbsp; the second SCA
- ...
- n-3: &nbsp;&nbsp;&nbsp; the last SCA
- n-2: &nbsp;&nbsp;&nbsp; the interconnect after the last SCA
- n-1: &nbsp;&nbsp;&nbsp; the outlet interconnect carrying twice the loop mass flow rate

### pipe_runner_diams
inner diameters in meters of the runners listed in L_runner. The first diameter is for the runner that carries half the total mass flow. Example diameters are: [^](#toc)
* 2 field sections = {x1}
* 4 field sections = {x1, x1}
* 6 field sections = {x1, x2}
* 8 field sections = {x1, x1, x3}
* 10 field sections = {x1, x4, x5}

### pipe_runner_expansions
number of expansions or contractions in the given runner section [^](#toc)

### pipe_runner_lengths
lengths in meters of the different diameter runners that extend away from the power block in one direction. L_runner[0] is currently defaulted to 25, which is for the runner piping in and around the power block before it heads out to the field in the main runners. L_runner[0] was previously shared with the other identical set of runners for the other half of the solar field, but this is no longer the case. The runner lengths include expansion loops, except for L_runner[0]. For a given row spacing, SCA length, gap between SCAs, and number of SCA's, example values are: [^](#toc)
* 2 field sections = {L_rnr_pb}
* 4 field sections = {L_rnr_pb, x}
* 6 field sections = {L_rnr_pb, 2x}
* 8 field sections = {L_rnr_pb, x, 2x}
* 10 field sections = {L_rnr_pb, 2x, 2x}

### pipe_runner_mdot_dsn
mass flow rate in kg/s of the heat transfer fluid in each runner section at design conditions. The first value is in the section in and around the power block before it heads out to the field in the main runners. The last value is in the section in and around the power block after it comes back from the field. The mass flow for the cold runner sections is the same as that entering the section, and the mass flow for the hot runner sections is the same as that leaving the section. [^](#toc)

### pipe_runner_P_dsn
pressure in bar of the heat transfer fluid entering each runner section at design conditions. The first value is for the section in and around the power block before it heads out to the field in the main runners. The last value is in the section in and around the power block after it comes back from the field. [^](#toc)

### pipe_runner_T_dsn
temperature in Celsius of the heat transfer fluid entering each runner section at design conditions. The first value is for the section in and around the power block before it heads out to the field in the main runners. The last value is in the section in and around the power block after it comes back from the field. [^](#toc)

### pipe_runner_vel_dsn
velocity in m/s of the heat transfer fluid in each runner section at design conditions. The first value is in the section in and around the power block before it heads out to the field in the main runners. The last value is in the section in and around the power block after it comes back from the field. The velocity for the cold runner sections is the same as that entering the section, and the velocity for the hot runner sections is the same as that leaving the section. [^](#toc)

### pipe_runner_wallthk
wall thickness of runner pipe sections in [m] [^](#toc)

### pipe_sgs_diams
SGS pipe inner diameters in [m] [^](#toc)
 
### pipe_sgs_mdot_dsn
SGS mass flow in each pipe section in [kg/s] [^](#toc)

### pipe_sgs_P_dsn
SGS pressure in each pipe section in [bar] [^](#toc)

### pipe_sgs_T_dsn
SGS temperature in each pipe section in [C] [^](#toc)
 
### pipe_sgs_vel_dsn
SGS velocity in each pipe section in [m/s] [^](#toc)
 
### pipe_sgs_wallthk
SGS wall thickness of each pipe section in [m] [^](#toc)
