The environmental layers used in the development of the Desert Tortoise habitat model are being distributed in a GIS Raster format in an ESRI ASCII Grid file format. This consists of 4 files. 

e.g.
1. annProx.asc - contains the raster file data
2. annProx.prj - contains the projection information
3. annProx.asc.aux.xml - contains the color table rules
4. annProx.xml - contains the FGDC standard metadata


There are 16 environmental layers included that are described in table 1 of Nussear, K.E., Esque, T.C., Inman, R.D., Gass, L., Thomas, K.A., Wallace, C.S.A., Blainey, J.B., Miller, D.M., and Webb, R.H., 2009, Modeling habitat of the desert tortoise (Gopherus agassizii) in the Mojave and parts of the Sonoran Deserts of California, Nevada, Utah, and Arizona: U.S. Geological Survey Open-File Report 2009-1102, 18 p



These GIS coverages can be opened using a variety of GIS software, including (but not limited to)

1. QGIS - open source: multiple platforms - opens ASC files natively
2. GRASS - open source: multiple platforms - use the 'r.in.gdal', or 'r.in.arc' functions
3. ArcGIS - proprietary: MS Windows only: Use ArcTools's Import ASCII to GRID function
4. ArcView - proprietary: MS Windows only: use the import ASCII Grid function (May need Spatial Analyst)


For those software packages that do not natively support ESRI ASCII Grids, a full description of the data format follows. There are many resources available on the internet that convert ACSCII grid formats to xyz triplets.

Definition of the ESRI ASCII Grid Format
(Copied from the ArcWorkstation 8.3 Help File)
The ASCII file must consist of header information containing a set of keywords, followed by cell values in row-major order. The file format is:
<NCOLS xxx>
<NROWS xxx>
<XLLCENTER xxx | XLLCORNER xxx>
<YLLCENTER xxx | YLLCORNER xxx>
<CELLSIZE xxx>
{NODATA_VALUE xxx}
row 1
row 2
.
.
.
row n
where xxx is a number, and the keyword nodata_value is optional and defaults to -9999. Row 1 of the data is at the top of the grid, row 2 is just under row 1 and so on.


For example, the following is an excerpt from our desert tortoise model output, where the raster file consists of 566 columns and 729 rows, with lower left coordinates (UTM NAD 83) of 314025 E, 3594742 N, a cell size of 1 km2, and a no data value of -9.9999999338158125107e+36. Data values in the rows follow the header, and include values of no-data, and data ranging from 0 to 1 in increments of 0.1

ncols        566
nrows        729
xllcorner    314025.779199999990
yllcorner    3594742.629285299685
cellsize     999.911660000000
NODATA_value -9.9999999338158125107e+36
...
 -9.9999999338158125107e+36 0  0  0  0  0  0  0.2  0.5  0.5  0.9  1  0.7  0.4  0  0  0  0  0  0 
 0  0.1  0.1  0.4  0.1  0.1  0.2  0.1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 ...

The nodata_value is the value in the ASCII file to be assigned to those cells whose true value is unknown. In the grid they will be assigned the keyword NODATA, or NA.
Cell values are delimited by spaces. No carriage returns are necessary at the end of each row in the grid. The number of columns in the header is used to determine when a new row begins.
The number of cell values must be equal to the number of rows times the number of columns, or an error will be returned.