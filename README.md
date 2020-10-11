# sdpa-cvxopt
Convert sdp(positive semidefinite programming) problems in SDPA form to cvxopt sdp form (Python)


The SDP library ( http://euler.nmt.edu/~brian/sdplib/) includes a bunch of SDP in  SDPA form (http://euler.nmt.edu/~brian/sdplib/sdplib.pdf). The sdp_cnvrt.py file provides a function called conelp_form, converting the SDPA form to cvxopt.

Required packages:numpy, itertools, cvxopt

An example is provided in example.ipynb.
