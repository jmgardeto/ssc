# Instrucciones para mantener la documentación de SSC

Este documento explica cómo usar Doxygen y LaTeX para crear la documentación de la API SSC.

## Visión general


La generación de la documentación de SSC implica una serie de pasos, algunos de los cuales están automatizados y otros son automáticos de por si:

* La documentación de SSC consta de 8 capítulos en el archivo ssc_guide.tex.

* El contenido de los Capítulos 1-7 está en ssc_guide.tex. Fue redactado originalmente por Aron Dobos y revisado por Paul Gilman.

* Doxygen genera contenido para el Capítulo 8 (la referencia del archivo para sscapi.h} en sscapi_8h.tex a partir de comentarios en sscapi.h.

* Un comando include en el archivo ssc_guide.tex inserta el contenido de sscapi_8h.tex en ssc_guide.tex.

* Los ajustes en el archivo de configuración de Doxygen Doxyfile determinan lo que se incluye en la referencia de funciones


Pasos generales para mantenimiento de la documentación
----------------------------------------------------------

1. Según sea necesario, revise el contenido en doc/ssc_guide.tex y los comentarios con formato Doxygen en ssc/sscapi.h.

2. Ejecute ``doxygen`` en Doxyconfig para generar la carpeta de látex con sscapi_8h.tex y otros archivos

3. Ejecute ``pdflatex`` en ssc_guide.tex para generar un archivo PDF; ejecútelo dos veces para generar la tabla de contenido, referencias cruzadas, etc.
