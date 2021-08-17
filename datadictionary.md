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
diámetros internos personalizados para la tubería del cabezal según se lee en los archivos de salida modificados. Se utiliza si custom_sf_pipe_sizes es true. No cambie el número de valores (secciones) ya que esto resultará en un comportamiento impredecible del modelo.  [^](#toc)
 
### sf_hdr_lengths
longitudes personalizadas para la tubería del cabezal según se leen de los archivos de salida modificados. Se utiliza si custom_sf_pipe_sizes es verdadero. Cambiar las longitudes no afecta el diseño del campo. No cambie el número de valores (secciones) ya que esto resultará en un comportamiento impredecible del modelo. [^](#toc)
 
### sf_hdr_wallthicks
espesores de pared personalizados para la tubería del cabezal como se lee en los archivos de salida modificados. Se utiliza si custom_sf_pipe_sizes es verdadero. No cambie el número de valores (secciones) ya que esto resultará en un comportamiento impredecible del modelo. [^](#toc)
 
### sf_rnr_diams
Diámetros internos personalizados para la tubería del tramo como se lee en los archivos de salida modificados. Se utiliza si custom_sf_pipe_sizes es verdadero. No cambie el número de valores (secciones) ya que esto resultará en un comportamiento impredecible del modelo. [^](#toc)
 
### sf_rnr_lengths
longitudes personalizadas para tramo de tubería como se leen de los archivos de salida modificados. Se utiliza si custom_sf_pipe_sizes es verdadero. Cambiar las longitudes no afecta el diseño del campo. No cambie el número de valores (secciones) ya que esto resultará en un comportamiento impredecible del modelo. [^](#toc)
 
### sf_rnr_wallthicks
espesores de pared personalizados para tramo de tubería como se lee en los archivos de salida modificados. Se utiliza si custom_sf_pipe_sizes es verdadero. No cambie el número de valores (secciones) ya que esto resultará en un comportamiento impredecible del modelo. [^](#toc)
 
### sgs_diams
Diámetros internos personalizados para la tubería SGS como se lee en los archivos de salida modificados. Se utiliza si custom_sgs_pipe_sizes es verdadero. No cambie el número de valores (secciones) ya que esto resultará en un comportamiento impredecible del modelo. [^](#toc)

Secciones Collection:
- 0: &nbsp;&nbsp;&nbsp; Cabezal de succión de la bomba de campo solar (SF) a la entrada de la bomba SF individual
- 1: &nbsp;&nbsp;&nbsp; Descarga de bomba SF individual al cabezal de descarga de la bomba SF
- 2: &nbsp;&nbsp;&nbsp; Cabezal de descarga de la bomba SF a los cabezales de la sección del campo de recolección (es decir, canales)
- 3: &nbsp;&nbsp;&nbsp; Cabezales de salida de la sección del campo del colector (es decir, canales) al vaso de expansión (almacenamiento indirecto) o al tanque de almacenamiento térmico caliente (almacenamiento directo)
- 4: &nbsp;&nbsp;&nbsp; Bypass branch - Ramal de derivación: cabezales de salida de la sección del campo del colector (es decir, canales) al cabezal de succión de la bomba (indirecto) o al tanque de almacenamiento térmico frío (directo)

Secciones Generation:
- 5: &nbsp;&nbsp;&nbsp; Cabezal de succión de bomba SGS a entrada de bomba SGS individual (aplicable solo para almacenamiento en serie con SF)
- 6: &nbsp;&nbsp;&nbsp; Descarga de bomba SGS individual al cabezal de descarga de bomba SGS (solo para almacenamiento en serie)
- 7: &nbsp;&nbsp;&nbsp; Cabezal de descarga de la bomba SGS al cabezal de suministro del generador de vapor (solo para almacenamiento en serie)
- 8: &nbsp;&nbsp;&nbsp; Cabezal de suministro del generador de vapor a la tubería entre generadores de vapor
- 9: &nbsp;&nbsp;&nbsp; Tubería entre generadores de vapor al cabezal de salida del generador de vapor
- 10: &nbsp;&nbsp;&nbsp; Cabezal de salida del generador de vapor al cabezal de succión de la bomba SF (indirecto) o al tanque de almacenamiento térmico frío (directo)

### sgs_lengths
longitud de la tubería en el circuito de flujo de recolección SGS seguida por el circuito de flujo de generación [m]. Estos no se leen de los archivos de salida modificados. Valores predeterminados = {0, 90, 100, 120, 0, 0, 0, 0, 80, 120, 80}. Las longitudes en los índices 0, 1, 5 y 6 son las longitudes sumadas de las múltiples secciones individuales de la bomba. No cambie el número de valores (secciones) ya que esto resultará en un comportamiento impredecible del modelo. [^](#toc)

### sgs_wallthicks
espesores de pared personalizados para la tubería SGS como se lee en los archivos de salida modificados. Se utiliza si custom_sgs_pipe_sizes es verdadero. No cambie el número de valores (secciones) ya que esto resultará en un comportamiento impredecible del modelo. [^](#toc)

### tanks_in_parallel
true si el ramal del tanque de almacenamiento frío y caliente está en paralelo con el campo solar (caso tradicional), o false si los tanques están en serie con el campo solar (solo aplicable para almacenamiento directo). Predeterminado = true. [^](#toc)

### T_tank_hot_inlet_min
la temperatura mínima de campo htf que puede entrar en el tanque caliente [C]. Por debajo de esta temperatura se abre la válvula de derivación y el campo recircula. Predeterminado = 400 C. [^](#toc)
				
### V_hdr_cold_max
velocidad máxima permitida en el cabezal frío en condiciones de diseño. Este valor puede excederse si también se excede el mínimo, pero solo si esto lo pone menos fuera de rango. [^](#toc)

### V_hdr_cold_min
velocidad mínima permitida en el cabezal frío en condiciones de diseño. Este valor puede excederse si también se excede el máximo, pero solo si esto lo pone menos fuera de rango. [^](#toc)

### V_hdr_hot_max
velocidad máxima permitida en el cabezal caliente en condiciones de diseño. Este valor puede excederse si también se excede el mínimo, pero solo si esto lo pone menos fuera de rango. [^](#toc)

### V_hdr_hot_min
velocidad mínima permitida en el cabezal caliente en condiciones de diseño. Este valor puede excederse si también se excede el máximo, pero solo si esto lo pone menos fuera de rango. [^](#toc)

### V_tes_des
velocidad del punto de diseño para dimensionar los diámetros de la tubería TES [m/s]. Predeterminado = 1.85 m/s. [^](#toc)


## Salidas
### pipe_header_diams
diámetros internos en metros de todas las secciones del cabezal en los cabezales fríos y calientes en un subcampo. El primer diámetro es el que está antes del primer conjunto de bucles en el cabezal frío y el último diámetro es el que sigue al último juego de bucles en el cabezal caliente. [^](#toc)

### pipe_header_expansions
número de expansiones o contracciones en la sección de cabezal dada [^](#toc)

### pipe_header_lengths
longitudes en metros de todas las secciones del cabezal, incluidas las longitudes agregadas de cualquier bucle de expansión. La primera longitud es la que está antes del primer conjunto de bucles en el encabezado en frío y la última longitud es la que sigue al último conjunto de bucles en el encabezado en caliente.  [^](#toc)

### pipe_header_mdot_dsn
caudal másico en kg/s del fluido caloportador en cada sección del colector en las condiciones de diseño. El primer valor está en la sección antes del primer conjunto de bucles en el cabezal frío y el último valor está en la sección después del último conjunto de bucles en el cabezal caliente. El caudal másico de las secciones de colector frío es el mismo que el que entra en la sección, y el caudal másico de las secciones de colector caliente es el mismo que el que sale de la sección. [^](#toc)

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
