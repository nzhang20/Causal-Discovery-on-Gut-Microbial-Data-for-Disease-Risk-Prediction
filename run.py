#!/usr/bin/env python

import sys
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import pylab as py
import itertools
import glob, os
from causallearn.search.ConstraintBased.CDNOD import cdnod
from causallearn.utils.GraphUtils import GraphUtils
from causallearn.utils.PCUtils.BackgroundKnowledge import BackgroundKnowledge
from causallearn.graph.GraphNode import GraphNode
import pydot
from IPython.display import Image, display


from src.eda import *
from src.graph import *


def main(targets):
    if "data" in targets:
        # clean data

    if "eda" in targets:
        # generate plots

    if "graph" in targets:
        data = pd.read_csv("data/clean.csv")

    if "clean" in targets:
        for f in glob.glob("data/*.csv"):
            os.remove(f)
        for f in glob.glob("plots/*.png"):
            os.remove(f)
        for f in glob.glob("graphs/*.png"):
            os.remove(f)

if __name__ == '__main__':
    targets = sys.argv[1:]
    if "all" in targets:
        main(["data", "eda", "graph"])
    else:
        main(targets)