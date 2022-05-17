# MoleWhat: The molecule name game

## Affiliation

Final Project for the Software Engineering 2022-2 class, taught by Dr. José Alfredo Noriega Carmona at the ([National Autonomous University of Mexico](https://www.unam.mx/)), in its  National School of Superior-Level Studies, _Morelia_ Campus (ENES Morelia), as part of its _Bs. in Information Technologies applied to Science_ career plan.

> DEVELOPED BY:
> 
> Mariela Yael Arias Rojo ([marielaAriass](https://github.com/marielaAriass))
> 
> Miriam Guadalupe Valdez Maldonado ([mirluvams](https://github.com/mirluvams))
> 
> Karime Ochoa Jacinto ([Kadkam8a](https://github.com/Kadkam8a))
> 
> Luis Aaron Nieto Cruz ([LuisAaronNietoCruz](https://github.com/LuisAaronNietoCruz))
> 
> Anton Pashkov ([anton-pashkov](https://github.com/anton-pashkov))


## Introduction
The following project intends to provide a functional endpoint in which a potential user can visualize a map containing worldwide flight routes in real time. 

The need for such a project stems as a desire to function as a real-time alert system with which to calculate the range in which a solar flare can cause interference with the equipment on board of airships currently flying at high altitudes, given their vulnerability to EMP events.

As a final project, this aims to be a practical demostration of how a series of computers in parallel can be used to obtain an easily scalable product that could potentially be commercialized and expanded upon once completed, without relying on a single, high performance computer to perform the entire process.

## Objectives
The expected output of this project is a set of four servers, each of which provide an essential part of the project. Their denomination is as follows:

> Data Retrieval Server

In charge of obtaining real-time data from the OpenSky Network API, through the use of Python. Will relay said information in a timely manner to the Storage Server, while also keeping copies of recent data points as required for archive purposes.

> Storage Server

In charge of storing up-to-date data with PostgreSQL, to be sorted and retrieved as needed by the Processing Server.

> Processing Server

In charge of generating map images on demand, based on the request by the Web Server API

> Web Server

In charge of displaying map images to the end user through a modern web interface


## Toolset
The project is to be developed by making use of modern Python 3 libraries, including but not limited to:

* [python 3](https://www.python.org/downloads/) | v 3.10.4
* [Guizero](https://pypi.org/project/guizero/) | v 1.3.0 
* [Sqlite3] | v. 3.36.0
* [drkit](https://www.rdkit.org/)

Additionally, it plans to make use of the following web technologies, on top of the usual development stack:
* [Apache HTTPD](https://httpd.apache.org/)
* [MariaDB](https://mariadb.org/)
* [PHP](https://www.php.net/)
* [Bootstrap](https://getbootstrap.com/)

## Methodology
...

## Usage Instructions & Requirements
...


The contents of this repository are licensed under the GNU General Public License version 3. Visit https://www.gnu.org/licenses/gpl-3.0.html for more information.
