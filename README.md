tessplan is a tool to schedule observations of TESS objects of interest.
It enables the user to determine visible transits of a TOI from any observatory
at any future time, and can determine in a given time window what planets will 
be transiting. It retrieves information from ExoFOP-TESS, so a user can use that
information (such as a TFOP Working Group prioritization) to select possible
planet candidates to observe. 
</p>
To install tessplan with pip:

        pip install tessplan
        

Alternatively you can install the current development version of tessplan:

        git clone https://github.com/benmontet/tessplan
        cd tessplan
        python setup.py install


The [notebooks](../../tree/master/notebooks) directory of this repository hosts a 
Jupyter notebook demonstrating use cases for tessplan.

Citations
---------

If tessplan is useful to your research, and you're feeling especially generous, 
you can cite the astropy-affiliated package 
[astroplan](https://github.com/astropy/astroplan), which handles some of the
under-the-hood computations within tessplan:

```
@ARTICLE{astroplan2018,
   author = {{Morris}, B.~M. and {Tollerud}, E. and {Sip{\H o}cz}, B. and
    {Deil}, C. and {Douglas}, S.~T. and {Berlanga Medina}, J. and
    {Vyhmeister}, K. and {Smith}, T.~R. and {Littlefair}, S. and
    {Price-Whelan}, A.~M. and {Gee}, W.~T. and {Jeschke}, E.},
    title = "{astroplan: An Open Source Observation Planning Package in Python}",
  journal = {\aj},
archivePrefix = "arXiv",
   eprint = {1712.09631},
 primaryClass = "astro-ph.IM",
 keywords = {methods: numerical, methods: observational },
     year = 2018,
    month = mar,
   volume = 155,
      eid = {128},
    pages = {128},
      doi = {10.3847/1538-3881/aaa47e},
   adsurl = {http://adsabs.harvard.edu/abs/2018AJ....155..128M},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
