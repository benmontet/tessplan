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
