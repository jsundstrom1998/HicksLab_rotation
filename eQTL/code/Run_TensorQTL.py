import sys
import pandas as pd
import session_info
from pyhere import here
from pathlib import Path
import statsmodels.stats.multitest

import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis
print(f'PyTorch {torch.__version__}')
print(f'Pandas {pd.__version__}')

from importlib import metadata

for dist in metadata.distributions():
    print(f"{dist.name}=={dist.version}")