#HOW TO USE
Must have the following packages installed: tinyec, pycryptodome, numdifftools, matplotlib, pandas, numpy.

First open the file `test-reflect_refract.py` and find the lines indicating where various objects representing curves are appended to the list initialObjs. The user may add, remove, or modify these curves as they like so long as they use valid objects found in `reflect_refract.py` and `general_curves.py` (with the exception of Object and ParentLens as these classes are not meant for usage and will fail). Be sure to check the objects for the appropriate variables for each type of object. With the exception of Lens and Linear_Lens, all of the classes can be either "reflection" or "refraction". Once all the objects are added into the list, the user can run the following in the terminal.

```bash
python3 test-reflect_refract.py
```

Additionally, the 3D model can be run as well from the file `test-reflect_refract3d.py` using the same instructions as for the 2D model. Note currently only the ellipsoid objects can be plotted. Once ready, the following command can be run:

```bash
python3 test-reflect_refract3d.py
```