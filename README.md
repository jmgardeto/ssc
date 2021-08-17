# SSC (SAM Simulation Core)
[![Build Status](https://travis-ci.com/NREL/ssc.svg?branch=develop)](https://travis-ci.com/NREL/ssc)
[![FOSSA Status](https://app.fossa.io/api/projects/git%2Bgithub.com%2FNREL%2Fssc.svg?type=shield)](https://app.fossa.io/projects/git%2Bgithub.com%2FNREL%2Fssc?ref=badge_shield)

El repositorio de proyectos de código abierto de SSC contiene el código fuente para la tecnología y los modelos financieros contenidos en el Modelo de Asesor de Sistemas (SAM) del Laboratorio Nacional de Energía Renovable. Para obtener más detalles sobre las capacidades de SAM, consulte el sitio web de SAM en [https://sam.nrel.gov/](https://sam.nrel.gov).

Podría pensar en SSC como el hogar de los algoritmos del programa de escritorio SAM. La mayoría de las personas ejecutan el código a través de la interfaz de usuario de escritorio, pero SSC también se puede ejecutar directamente usando el [SAM Sofware Develoment Kit](https://sam.nrel.gov/sdk).

SSC requiere la construcción de otros proyectos de código abierto:

- [Google Test](https://github.com/google/googletest)
- [LK](https://github.com/nrel/lk)
- [wxWidgets](https://www.wxwidgets.org/)
- [WEX](https://github.com/nrel/wex)
- [jsoncpp](https://github.com/open-source-parsers/jsoncpp)

Sin embargo, si elimina SDKtool y TCSconsole de su proyecto SSC, puede crear SSC sin otras dependencias de software. Consulte la [wiki del proyecto SAM](https://github.com/NREL/SAM/wiki) principal para obtener instrucciones de compilación completas y dependencias de software.

SSC incluye directamente el código fuente de otros tres proyectos de código abierto y los construye como parte de su proceso de construcción. Estos proyectos y sus respectivas licencias son:
- [NLopt](https://nlopt.readthedocs.io/en/latest/) - codigo [aqui](https://github.com/NREL/ssc/tree/develop/nlopt), [licencia LGPL](https://nlopt.readthedocs.io/en/latest/NLopt_License_and_Copyright/)
- [lp_solve](http://lpsolve.sourceforge.net/5.5/) - codigo [aqui](https://github.com/NREL/ssc/tree/develop/lpsolve), [licencia LGPL](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html)
- [splinter](https://github.com/bgrimstad/splinter) - codigo [aqui](https://github.com/NREL/ssc/tree/develop/splinter), [licencia MPL](https://github.com/bgrimstad/splinter/blob/master/LICENSE)


Para explorar el código y comprender los algoritmos utilizados en SSC, comience por buscar en el proyecto "SSC" los módulos de cómputo (archivos que comienzan con cmod_) para encontrar el módulo de cómputo para la tecnología o el modelo financiero de interés.

# Contribuyendo

Consulte las pautas de contribución en el [README del proyecto SAM](https://github.com/NREL/SAM/blob/develop/README.md) principal.

# Licencia

SSC tiene licencia con los términos de la cláusula BSD-3, que se encuentran [aquí](https://github.com/NREL/SAM/blob/develop/LICENSE).

# Citando este paquete

System Advisor Model Version 2020.2.29 (2020.2.29). SSC source code. National Renewable Energy Laboratory. Golden, CO. Accessed May 27, 2020. https://github.com/NREL/ssc 
